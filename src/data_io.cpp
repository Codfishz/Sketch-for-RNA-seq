#include "data_io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <unordered_map>

/**
 * @brief Check whether a DNA sequence contains only valid characters (A, T, C, G).
 *
 * This function validates that the input sequence contains only the characters
 * 'A', 'T', 'C', and 'G'. It uses a static boolean array for fast lookup.
 *
 * @param sequence The DNA sequence string to validate.
 * @return true if the sequence is valid, false otherwise.
 */
bool is_valid_sequence(const std::string& sequence) {
    static bool valid_char[256] = { false };

    // Initialize valid characters on first call
    if (!valid_char['A']) {
        valid_char['A'] = true;
        valid_char['T'] = true;
        valid_char['C'] = true;
        valid_char['G'] = true;
    }

    for (char c : sequence) {
        if (!valid_char[static_cast<unsigned char>(c)]) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Load transcripts from a FASTA file.
 *
 * This function reads a FASTA file and constructs a map from transcript ID to Transcript.
 * Each transcript is only added if its sequence is valid (only contains A, T, C, G).
 *
 * @param fasta_file The file path of the FASTA file.
 * @return An unordered_map where the key is the transcript ID and the value is the Transcript struct.
 *
 * @throw std::runtime_error if the file cannot be opened.
 */
std::unordered_map<std::string, Transcript> load_fasta(const std::string& fasta_file) {
    std::unordered_map<std::string, Transcript> transcripts;
    std::ifstream infile(fasta_file, std::ios::in | std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Could not open FASTA file: " + fasta_file);
    }

    std::string line, sequence;
    std::string current_id;
    infile.rdbuf()->pubsetbuf(nullptr, 1 << 20); // Set 1MB buffer for faster reading

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save previous transcript if valid
            if (!current_id.empty() && is_valid_sequence(sequence)) {
                transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
            }
            // Extract transcript ID (assumes ID is terminated by a space)
            current_id = line.substr(1, line.find(' ') - 1);
            sequence.clear();
            sequence.reserve(1000);  // Preallocate buffer to avoid reallocations
        } else {
            sequence += line;
        }
    }
    // Save the last transcript
    if (!current_id.empty()) {
        transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
    }

    return transcripts;
}

/**
 * @brief Load reads from a FASTQ file.
 *
 * This function reads a FASTQ file and loads the reads into an unordered_map.
 * Only reads with valid DNA sequences (only contains A, T, C, G) are kept.
 * It also prints out the total number of reads encountered (including invalid ones).
 *
 * @param fastq_file The file path of the FASTQ file.
 * @return An unordered_map where the key is the read ID and the value is the Read struct.
 *
 * @throw std::runtime_error if the file cannot be opened.
 */
std::unordered_map<std::string, Read> load_fastq(const std::string& fastq_file) {
    int count = 0;
    std::unordered_map<std::string, Read> reads;
    std::ifstream infile(fastq_file);
    if (!infile) {
        throw std::runtime_error("Could not open FASTQ file: " + fastq_file);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] != '@') continue;

        Read read;
        read.id = line.substr(1);
        std::getline(infile, read.sequence);
        std::getline(infile, line); // Skip "+" line
        std::getline(infile, read.quality);
        count++;
        if (is_valid_sequence(read.sequence)) {
            reads[read.id] = read;
        }
    }
    std::cout << "Actual number of reads: " << count << std::endl;
    return reads;
}

/**
 * @brief Output transcript quantification results to a CSV file.
 *
 * This function writes the quantification results into a CSV file with three columns:
 * transcript Name, NumReads, and EM_Abundance. It iterates through the transcripts map
 * and outputs an entry only if corresponding read count and EM abundance data are available.
 *
 * @param filename The output CSV file path.
 * @param read_counts An unordered_map with transcript IDs as keys and read counts as values.
 * @param pi An unordered_map with transcript IDs as keys and EM-based abundance estimates as values.
 * @param transcripts An unordered_map with transcript IDs as keys and Transcript structs as values.
 *
 * @throw std::runtime_error if the output file cannot be opened.
 */
void output_to_csv(const std::string& filename, 
                   const std::unordered_map<std::string, double>& read_counts,
                   const std::unordered_map<std::string, double>& pi,
                   const std::unordered_map<std::string, Transcript>& transcripts) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    outfile << "Name,NumReads,EM_Abundance\n";
    for (const auto& [id, transcript] : transcripts) {
        auto read_count_it = read_counts.find(id);
        auto pi_it = pi.find(id);
        if (read_count_it != read_counts.end() && pi_it != pi.end()) {
            outfile << id << "," << read_count_it->second << "," << pi_it->second << "\n";
        }
    }

    outfile.close();
}

/**
 * @brief Save the full index to a binary file for efficient reloading.
 *
 * This function serializes k-mer lengths, transcripts, and k-mer to transcript mappings into a binary file.
 * The index file can later be loaded to avoid rebuilding the index from scratch.
 *
 * @param index_output_path The file path where the index will be saved.
 * @param kmer_lengths A vector of k-mer lengths used in the index.
 * @param kmer_to_transcripts An unordered_map mapping each k-mer length to a TranscriptMapping.
 * @param transcripts An unordered_map of transcripts.
 */
void save_index(const std::string& index_output_path,
                std::vector<unsigned>& kmer_lengths,
                const std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
                const std::unordered_map<std::string, Transcript>& transcripts) {
    std::ofstream ofs(index_output_path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Unable to open file for writing: " << index_output_path << std::endl;
        return;
    }

    // Save k-mer lengths
    size_t kmerCount = kmer_lengths.size();
    ofs.write(reinterpret_cast<const char*>(&kmerCount), sizeof(kmerCount));
    for (unsigned k : kmer_lengths) {
        ofs.write(reinterpret_cast<const char*>(&k), sizeof(k));
    }

    // Save transcript data
    size_t transcriptCount = transcripts.size();
    ofs.write(reinterpret_cast<const char*>(&transcriptCount), sizeof(transcriptCount));
    for (const auto& [id, transcript] : transcripts) {
        size_t idLen = id.size();
        ofs.write(reinterpret_cast<const char*>(&idLen), sizeof(idLen));
        ofs.write(id.data(), idLen);

        size_t seqLen = transcript.sequence.size();
        ofs.write(reinterpret_cast<const char*>(&seqLen), sizeof(seqLen));
        ofs.write(transcript.sequence.data(), seqLen);

        ofs.write(reinterpret_cast<const char*>(&transcript.length), sizeof(transcript.length));
    }

    // Save k-mer to transcript mappings
    size_t mapCount = kmer_to_transcripts.size();
    ofs.write(reinterpret_cast<const char*>(&mapCount), sizeof(mapCount));
    for (const auto& [kmer_length, mapping] : kmer_to_transcripts) {
        ofs.write(reinterpret_cast<const char*>(&kmer_length), sizeof(kmer_length));
        size_t mappingSize = mapping.size();
        ofs.write(reinterpret_cast<const char*>(&mappingSize), sizeof(mappingSize));
        for (const auto& [kmer, vec] : mapping) {
            ofs.write(reinterpret_cast<const char*>(&kmer), sizeof(kmer));
            size_t vecSize = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));
            for (const auto& pair : vec) {
                const std::string& transcript_id = pair.first;
                size_t tidLen = transcript_id.size();
                ofs.write(reinterpret_cast<const char*>(&tidLen), sizeof(tidLen));
                ofs.write(transcript_id.data(), tidLen);
                // Sketch pointers are not saved (set to nullptr on load)
            }
        }
    }

    ofs.close();
    std::cout << "Index saved to " << index_output_path << std::endl;
}

/**
 * @brief Load index data from a binary file into memory.
 *
 * This function deserializes the binary index file and restores k-mer lengths,
 * transcript data, and k-mer to transcript mappings.
 *
 * @param index_path The file path of the binary index file.
 * @param kmer_lengths A vector to store the loaded k-mer lengths.
 * @param kmer_to_transcripts An unordered_map to store the loaded k-mer to transcript mappings.
 * @param transcripts An unordered_map to store the loaded transcript data.
 */
void load_index(const std::string& index_path,
                std::vector<unsigned>& kmer_lengths,
                std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
                std::unordered_map<std::string, Transcript>& transcripts) {
    std::ifstream ifs(index_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: Unable to open file for reading: " << index_path << std::endl;
        return;
    }

    // Load k-mer lengths
    size_t kmerCount = 0;
    ifs.read(reinterpret_cast<char*>(&kmerCount), sizeof(kmerCount));
    kmer_lengths.resize(kmerCount);
    for (size_t i = 0; i < kmerCount; i++) {
        unsigned k;
        ifs.read(reinterpret_cast<char*>(&k), sizeof(k));
        kmer_lengths[i] = k;
    }

    // Load transcript data
    size_t transcriptCount = 0;
    ifs.read(reinterpret_cast<char*>(&transcriptCount), sizeof(transcriptCount));
    transcripts.clear();
    for (size_t i = 0; i < transcriptCount; i++) {
        size_t idLen = 0;
        ifs.read(reinterpret_cast<char*>(&idLen), sizeof(idLen));
        std::string id(idLen, '\0');
        ifs.read(&id[0], idLen);

        size_t seqLen = 0;
        ifs.read(reinterpret_cast<char*>(&seqLen), sizeof(seqLen));
        std::string sequence(seqLen, '\0');
        ifs.read(&sequence[0], seqLen);

        int length;
        ifs.read(reinterpret_cast<char*>(&length), sizeof(length));

        transcripts[id] = Transcript{id, sequence, length};
    }

    // Load k-mer to transcript mappings
    size_t mapCount = 0;
    ifs.read(reinterpret_cast<char*>(&mapCount), sizeof(mapCount));
    kmer_to_transcripts.clear();
    for (size_t i = 0; i < mapCount; i++) {
        unsigned kmer_length;
        ifs.read(reinterpret_cast<char*>(&kmer_length), sizeof(kmer_length));
        size_t mappingSize = 0;
        ifs.read(reinterpret_cast<char*>(&mappingSize), sizeof(mappingSize));
        TranscriptMapping mapping;
        for (size_t j = 0; j < mappingSize; j++) {
            uint32_t kmer;
            ifs.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
            size_t vecSize = 0;
            ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
            std::vector<std::pair<std::string, const SketchType*>> vec;
            for (size_t h = 0; h < vecSize; h++) {
                size_t tidLen = 0;
                ifs.read(reinterpret_cast<char*>(&tidLen), sizeof(tidLen));
                std::string transcript_id(tidLen, '\0');
                ifs.read(&transcript_id[0], tidLen);
                vec.emplace_back(transcript_id, nullptr);  // Sketch pointer restored as nullptr
            }
            mapping[kmer] = vec;
        }
        kmer_to_transcripts[kmer_length] = mapping;
    }

    ifs.close();
    std::cout << "Index loaded from " << index_path << std::endl;
}
