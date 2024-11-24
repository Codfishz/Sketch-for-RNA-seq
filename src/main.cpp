#include "data_io.h"
#include "kmer.h"
#include "sketch.h"
#include "isoform_assignment.h"
#include "sparse_chaining.h"
//#include "MurmurHash3.h"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <chrono>
#include <unordered_map> 
//
void print_help(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS] <reference_genome.fasta> <reads.fastq>" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help              Show this help message and exit" << std::endl;
    std::cout << "  -t, --threads NUM       Number of threads to use (default: 1)" << std::endl;
    std::cout << "  -m, --memory SIZE       Maximum memory to use (default: 1024 MB)" << std::endl;
    std::cout << "  -k, --kmer-length SIZE  Length of k-mers (default: 31)" << std::endl;
    std::cout << "  -w, --window-size SIZE  Window size for sparse chaining (default: 100)" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  <reference_genome.fasta>  Path to the reference genome in FASTA format" << std::endl;
    std::cout << "  <reads.fastq>             Path to the raw RNA-seq data in FASTQ format" << std::endl;
    std::cout << std::endl;
    std::cout << "Example:" << std::endl;
    std::cout << "  " << program_name << " -t 4 -m 2048 -k 21 reference.fasta reads.fastq" << std::endl;
}

int main(int argc, char* argv[]) {
    int threads = 1;
    int memory = 1024; // 以MB为单位
    int kmer_length = 31;
    int window_size = 100;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"threads", required_argument, 0, 't'},
        {"memory", required_argument, 0, 'm'},
        {"kmer-length", required_argument, 0, 'k'},
        {"window-size", required_argument, 0, 'w'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "ht:m:k:w:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_help(argv[0]);
                return 0;
            case 't':
                threads = std::stoi(optarg);
                break;
            case 'm':
                memory = std::stoi(optarg);
                break;
            case 'k':
                kmer_length = std::stoi(optarg);
                break;
            case 'w':
                window_size = std::stoi(optarg);
                break;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    if (optind + 3 > argc) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome.fasta> <reads.fastq> <gene_annotation.gtf>" << std::endl;
        return 1;
    }

    std::string reference_genome_path = argv[optind];
    std::string reads_path = argv[optind + 1];
    std::string annotation_path = argv[optind + 2];
    std::string output_path = argv[optind + 3];

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // auto annotations_map = parse_gtf(annotation_path);
        // std::unordered_map<std::string, Transcript> transcripts = load_transcripts_with_annotations(reference_genome_path, annotations_map);
        std::unordered_map<std::string, Transcript> transcripts = load_fasta(reference_genome_path);

        int count = 0;

        std::unordered_map<std::string, std::unordered_set<uint32_t>> transcript_sketches;
        for (const auto& [id, transcript] : transcripts) {
            if (transcript.sequence.length() < kmer_length) {
                // std::cerr << "Skipping transcript " << id << " because its length (" << transcript.sequence.length() 
                //         << ") is less than kmer-length (" << kmer_length << ")." << std::endl;
                continue; 
            }

            auto hashed_kmers = extract_and_hash_kmers_murmur(transcript.sequence, kmer_length,67890);
            auto sketch = createSketch_FracMinhash(hashed_kmers, 0.2);
            transcript_sketches[id] = sketch;
        }

        auto kmer_to_transcripts = build_kmer_to_transcript_map(transcript_sketches);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Number of transcripts loaded: " << transcripts.size() << std::endl;
        std::cout << "Time taken to load and process transcripts: " << elapsed.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();

        std::unordered_map<std::string, Read> reads = load_fastq(reads_path);
        
        count = 0;
        std::unordered_map<std::string, std::unordered_set<uint32_t>> read_sketches;
        for (const auto& [id, read] : reads) {
            if (read.sequence.length() < kmer_length) {
                // std::cerr << "Skipping read " << id << " because its length (" << read.sequence.length() 
                //         << ") is less than kmer-length (" << kmer_length << ")." << std::endl;
                continue; 
            }
            auto hashed_kmers = extract_and_hash_kmers_murmur(read.sequence, kmer_length, 67890);
            auto sketch = createSketch_FracMinhash(hashed_kmers, 0.2);
            read_sketches[id] = sketch;
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Number of reads loaded: " << reads.size() << std::endl;
        std::cout << "Time taken to load and process reads: " << elapsed.count() << " seconds." << std::endl;

        size_t num_reads = read_sketches.size();
        std::cout<<"Number of reads sketch:" <<num_reads<<std::endl;


        start = std::chrono::high_resolution_clock::now();

        auto homologous_segments = sparse_chain(read_sketches, kmer_to_transcripts, transcripts, reads, kmer_length, 0.3);
        
        std::cout << "Total number of successful matches: " << homologous_segments.size() << std::endl;
        std::cout<<"finish sparse chain"<<std::endl;

        // count =0;

        // for (const auto& [read_id, transcript_id] : homologous_segments) {
        //     std::cout << read_id << " -> " << transcript_id << std::endl;
        // }

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken to create sparse chain: " << elapsed.count() << " seconds." << std::endl;


        start = std::chrono::high_resolution_clock::now();

        std::unordered_map<std::string, int> read_counts = assign_reads_to_isoforms(homologous_segments, transcripts);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken for read-isoform assignment: " << elapsed.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();

        std::unordered_map<std::string, double> tpm = calculate_tpm(read_counts, transcripts);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken to calculate TPM: " << elapsed.count() << " seconds." << std::endl;

        // for (const auto& [transcript_id, tpm_value] : tpm) {
        //     std::cout << "Transcript " << transcript_id << " has TPM " << tpm_value << std::endl;
        // }

        output_to_csv(output_path, read_counts, tpm, transcripts);
        std::cout <<"output finish" << std::endl;

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
