#ifndef DATA_IO_H
#define DATA_IO_H

/**
 * @file data_io.h
 * @brief Declarations for data input/output functions and basic data structures.
 * @author ziruichen
 * @date 2025/04/07
 */

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

// Type alias for a sketch used in k-mer based algorithms.
using SketchType = std::unordered_set<uint32_t>;

// Type alias for mapping k-mers (represented as uint32_t) to transcript information.
// The vector stores pairs of transcript IDs and pointers to their corresponding sketch data.
using TranscriptMapping = std::unordered_map<uint32_t, std::vector<std::pair<std::string, const SketchType*>>>;

/**
 * @brief Represents a transcript.
 *
 * Contains the transcript ID, sequence, and its length.
 */
struct Transcript {
    std::string id;       ///< Transcript identifier.
    std::string sequence; ///< Nucleotide sequence of the transcript.
    int length;           ///< Length of the transcript sequence.
};

/**
 * @brief Represents a sequencing read.
 *
 * Stores the read identifier, its nucleotide sequence, and quality string.
 */
struct Read {
    std::string id;       ///< Read identifier.
    std::string sequence; ///< Nucleotide sequence of the read.
    std::string quality;  ///< Quality string corresponding to the read.
};

/**
 * @brief Check whether a DNA sequence contains only valid characters (A, T, C, G).
 *
 * @param sequence The DNA sequence string.
 * @return true if the sequence is valid, false otherwise.
 */
bool is_valid_sequence(const std::string& sequence);

/**
 * @brief Load transcripts from a FASTA file.
 *
 * @param fasta_file The file path to the FASTA file.
 * @return An unordered_map mapping transcript IDs to Transcript structures.
 */
std::unordered_map<std::string, Transcript> load_fasta(const std::string& fasta_file);

/**
 * @brief Load sequencing reads from a FASTQ file.
 *
 * Only valid reads (containing only A, T, C, G) are kept.
 *
 * @param fastq_file The file path to the FASTQ file.
 * @return An unordered_map mapping read IDs to Read structures.
 */
std::unordered_map<std::string, Read> load_fastq(const std::string& fastq_file);

/**
 * @brief Output quantification results to a CSV file.
 *
 * Writes the transcript quantification results in a simplified CSV format containing:
 * Name, NumReads, and EM_Abundance.
 *
 * @param filename The output CSV file path.
 * @param read_counts A mapping from transcript IDs to the number of reads.
 * @param pi A mapping from transcript IDs to EM-based abundance estimates.
 * @param transcripts A mapping from transcript IDs to Transcript structures.
 */
void output_to_csv(const std::string& filename, 
    const std::unordered_map<std::string, double>& read_counts,
    const std::unordered_map<std::string, double>& pi,
    const std::unordered_map<std::string, Transcript>& transcripts);

/**
 * @brief Save the full index to a binary file.
 *
 * Serializes k-mer lengths, transcript data, and k-mer-to-transcript mappings for later use.
 *
 * @param index_output_path The file path to save the index.
 * @param kmer_lengths A vector of k-mer lengths used in the index.
 * @param kmer_to_transcripts A mapping from k-mer lengths to their corresponding TranscriptMapping.
 * @param transcripts A mapping from transcript IDs to Transcript structures.
 */
void save_index(const std::string& index_output_path,
    std::vector<unsigned>& kmer_lengths,
    const std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts);

/**
 * @brief Load index data from a binary file.
 *
 * Deserializes the binary file and restores k-mer lengths, transcript data, and mappings.
 *
 * @param index_path The file path to the binary index.
 * @param kmer_lengths A vector to store the loaded k-mer lengths.
 * @param kmer_to_transcripts A mapping to store the loaded TranscriptMapping.
 * @param transcripts A mapping to store the loaded transcript data.
 */
void load_index(const std::string& index_path,
    std::vector<unsigned>& kmer_lengths,
    std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
    std::unordered_map<std::string, Transcript>& transcripts);

#endif // DATA_IO_H
