#include "data_io.h"
#include "kmer.h"
#include "sketch.h"
#include "isoform_assignment.h"
#include "sparse_chaining.h"

#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <unordered_map>
#include <nthash/nthash.hpp>
#include <fstream>
#include <sstream>

/**
 * @brief Print the help message describing the program usage.
 *
 * Displays information about program modes (index, quant) and available command-line options.
 *
 * @param program_name The name of the executable.
 */
void print_help(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS] <mode> [arguments]" << std::endl;
    std::cout << "Modes:" << std::endl;
    std::cout << "  index   Build index from reference genome" << std::endl;
    std::cout << "  quant   Quantify using pre-built index and reads" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help              Show this help message and exit" << std::endl;
    std::cout << "  -k, --kmer-length SIZE  Comma separated list of k-mer lengths (default: 81)" << std::endl;
    std::cout << "  -o, --mode MODE         Mode: index or quant (default: quant)" << std::endl;
    std::cout << std::endl;
    std::cout << "Index mode usage:" << std::endl;
    std::cout << "  " << program_name << " index <reference_genome.fasta> <index_output>" << std::endl;
    std::cout << std::endl;
    std::cout << "Quant mode usage:" << std::endl;
    std::cout << "  " << program_name << " quant <index_file> <reads.fastq> <output>" << std::endl;
}

// Global constant for sketch size; adjust as needed.
const float sketch_size = 0.05f;

/**
 * @brief Build an index from a reference genome and save it to a file.
 *
 * This function loads transcripts from a FASTA file, creates a MultiKmerSketch for each transcript
 * by extracting k-mers (using a fractional MinHash approach), builds a mapping from k-mers to transcripts,
 * and finally saves the index to a binary file.
 *
 * @param reference_genome_path Path to the reference genome FASTA file.
 * @param index_output_path Path where the index file will be saved.
 * @param kmer_lengths A vector of k-mer lengths to be used.
 */
void build_and_save_index(const std::string& reference_genome_path,
                          const std::string& index_output_path,
                          std::vector<unsigned>& kmer_lengths)
{
    auto start = std::chrono::high_resolution_clock::now();
    // Load transcripts from FASTA file.
    std::unordered_map<std::string, Transcript> transcripts = load_fasta(reference_genome_path);
    std::unordered_map<std::string, MultiKmerSketch> transcript_sketches;
    
    // For each transcript, check validity and build its k-mer sketches.
    for (const auto& [id, transcript] : transcripts) {
        bool valid = true;
        for (unsigned k : kmer_lengths) {
            if (transcript.sequence.size() < k) {
                valid = false;
                break;
            }
        }
        if (!valid)
            continue;

        MultiKmerSketch mks;
        for (unsigned k : kmer_lengths) {
            mks.sketches[k] = createSketch_FracMinhash_direct(transcript.sequence, k, sketch_size);
        }
        transcript_sketches[id] = std::move(mks);
    }
    
    // Build mapping from k-mer hash values to transcripts.
    auto kmer_to_transcripts = build_kmer_to_transcript_map(transcript_sketches);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Index built in " << elapsed.count() << " seconds." << std::endl;

    // Save the index to file.
    save_index(index_output_path, kmer_lengths, kmer_to_transcripts, transcripts);
}

/**
 * @brief Process a FASTQ file in a single pass to build read sketches.
 *
 * Reads the FASTQ file line by line, validates each read, and creates a MultiKmerSketch for each read
 * using the provided k-mer lengths and sketch size.
 *
 * @param fastq_file Path to the FASTQ file.
 * @param effective_kmer_lengths A vector of k-mer lengths to be used.
 * @param sketch_size The fraction used in sketching.
 * @return An unordered_map mapping read IDs to their corresponding MultiKmerSketch.
 *
 * @throws std::runtime_error if the FASTQ file cannot be opened.
 */
std::unordered_map<std::string, MultiKmerSketch> process_fastq_single_pass(
    const std::string& fastq_file,
    const std::vector<unsigned>& effective_kmer_lengths,
    double sketch_size)
{
    std::unordered_map<std::string, MultiKmerSketch> read_sketches;
    std::ifstream infile(fastq_file);
    if (!infile) {
        throw std::runtime_error("Could not open FASTQ file: " + fastq_file);
    }
    
    std::string line;
    // Read the FASTQ file one record at a time.
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] != '@') {
            continue;  // Skip non-read header lines.
        }
        
        Read read;
        read.id = line.substr(1); // Remove '@'
        std::getline(infile, read.sequence);
        std::getline(infile, line); // Skip "+" line.
        std::getline(infile, read.quality);
        
        // Validate read sequence.
        if (!is_valid_sequence(read.sequence))
            continue;
        
        // Ensure the read is long enough for the largest k-mer.
        unsigned max_k = *std::max_element(effective_kmer_lengths.begin(), effective_kmer_lengths.end());
        if (read.sequence.size() < max_k)
            continue;
        
        MultiKmerSketch mks;
        // Create a sketch for each specified k-mer length.
        for (unsigned k : effective_kmer_lengths) {
            mks.sketches[k] = createSketch_FracMinhash_direct(read.sequence, k, sketch_size);
        }
        
        // Save the constructed sketch for the current read.
        read_sketches[read.id] = std::move(mks);
    }
    
    return read_sketches;
}

/**
 * @brief Quantify transcript abundance using a pre-built index and reads.
 *
 * This function loads a previously built index, processes a FASTQ file to create read sketches,
 * performs sparse chaining to identify homologous transcript segments, estimates isoform abundances
 * using the EM algorithm, assigns reads to isoforms, and finally outputs the results to a CSV file.
 *
 * @param index_path Path to the pre-built index file.
 * @param reads_path Path to the FASTQ file containing reads.
 * @param output_path Path where the output CSV file will be saved.
 * @param kmer_lengths A vector of k-mer lengths specified on the command line.
 */
void quantification(const std::string& index_path,
                    const std::string& reads_path,
                    const std::string& output_path,
                    std::vector<unsigned>& kmer_lengths)
{
    // Load the pre-built index.
    std::vector<unsigned> loaded_kmer_lengths;
    std::unordered_map<unsigned, TranscriptMapping> kmer_to_transcripts;
    std::unordered_map<std::string, Transcript> transcripts;
    load_index(index_path, kmer_lengths, kmer_to_transcripts, transcripts);

    std::cout << "Loading index completed" << std::endl;
    // Use loaded k-mer lengths if available, otherwise use command-line specified lengths.
    std::vector<unsigned> effective_kmer_lengths = loaded_kmer_lengths.empty() ? kmer_lengths : loaded_kmer_lengths;

    // Process FASTQ file to generate read sketches.
    auto read_sketches = process_fastq_single_pass(reads_path, effective_kmer_lengths, sketch_size);
    std::cout << "Loading read completed" << std::endl;

    // Perform sparse chaining and subsequent quantification steps.
    auto homologous_segments = sparse_chain(read_sketches, kmer_to_transcripts, transcripts, effective_kmer_lengths, 0.9);
    std::cout << "Sparse chaining completed" << std::endl;

    auto pi = estimate_isoform_abundance_em(homologous_segments, transcripts, 20, 0.01);
    std::cout << "EM estimation completed" << std::endl;

    auto read_counts = assign_reads_to_isoforms(homologous_segments, pi, transcripts);
    std::cout << "Read assignment completed" << std::endl;

    // Output the results to a CSV file.
    output_to_csv(output_path, read_counts, pi, transcripts);
    std::cout << "Output written to " << output_path << std::endl;
}

/**
 * @brief Main function to control index building and quantification.
 *
 * The program can run in two modes:
 * - "index": Build an index from a reference genome.
 * - "quant": Quantify isoform abundance using a pre-built index and reads.
 *
 * Command-line options allow specifying k-mer lengths and the mode of operation.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Exit status code.
 */
int main(int argc, char* argv[]) {
    // Default mode is "quant" with a default k-mer length of 31.
    std::string mode = "quant";
    std::vector<unsigned> kmer_lengths = { 31 };

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"kmer-length", required_argument, 0, 'k'},
        {"mode", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "hk:o:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_help(argv[0]);
                return 0;
            case 'k': {
                std::string s(optarg);
                kmer_lengths.clear();  // Clear default if user provides values.
                std::istringstream iss(s);
                std::string token;
                while (std::getline(iss, token, ',')) {
                    if (!token.empty()) {
                        kmer_lengths.push_back(std::stoi(token));
                    }
                }
                break;
            }
            case 'o':
                mode = std::string(optarg);
                break;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    // Process mode-specific arguments.
    if (mode == "index") {
        if (optind + 2 > argc) {
            std::cerr << "Usage: " << argv[0] << " index <reference_genome.fasta> <index_output>" << std::endl;
            return 1;
        }
        std::string reference_genome_path = argv[optind];
        std::string index_output_path = argv[optind + 1];
        build_and_save_index(reference_genome_path, index_output_path, kmer_lengths);
    } else if (mode == "quant") {
        if (optind + 3 > argc) {
            std::cerr << "Usage: " << argv[0] << " quant <index_file> <reads.fastq> <output>" << std::endl;
            return 1;
        }
        std::string index_path = argv[optind];
        std::string reads_path = argv[optind + 1];
        std::string output_path = argv[optind + 2];
        quantification(index_path, reads_path, output_path, kmer_lengths);
    } else {
        std::cerr << "Invalid mode. Please choose 'index' or 'quant'." << std::endl;
        return 1;
    }

    return 0;
}
