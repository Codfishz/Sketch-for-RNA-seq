#ifndef SKETCHING_H
#define SKETCHING_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstdint>

/**
 * @brief Alias for a set of k-mer hash values.
 *
 * This type represents the sketch for a sequence as an unordered_set of 32-bit hash values.
 */
using SketchType = std::unordered_set<uint32_t>;

/**
 * @brief Mapping from a k-mer hash value to transcript information.
 *
 * The key is the k-mer hash value (uint32_t). The value is a vector of pairs,
 * each containing a transcript ID (std::string) and a pointer to the corresponding SketchType.
 */
using TranscriptMapping = std::unordered_map<uint32_t, std::vector<std::pair<std::string, const SketchType*>>>;

/**
 * @brief Structure to store multiple k-mer sketches for a transcript.
 *
 * The sketches map stores sketches for different k-mer lengths.
 * The key is the k-mer length and the value is the corresponding SketchType.
 */
struct MultiKmerSketch {
    std::unordered_map<unsigned, SketchType> sketches;
};

/**
 * @brief Create a fractional MinHash sketch for a given DNA sequence.
 *
 * This function extracts k-mers from the input sequence using nthash (or equivalent method) and
 * selects a fraction of the k-mers based on the provided fraction value. Only k-mers with hash values
 * below a threshold (calculated as the maximum uint32_t multiplied by fraction) are included in the sketch.
 *
 * @param sequence The input DNA sequence.
 * @param k The length of the k-mers.
 * @param fraction The fraction threshold for including k-mers (e.g., 0.05 for 5%).
 * @return std::unordered_set<uint32_t> A set containing the hash values of the selected k-mers.
 */
std::unordered_set<uint32_t> createSketch_FracMinhash_direct(const std::string &sequence, int k, double fraction);

/**
 * @brief Build a mapping from k-mers to transcripts based on precomputed sketches.
 *
 * This function iterates over a collection of transcript sketches and constructs a mapping from k-mer hash values
 * to transcript information. The outer key of the returned map is the k-mer length, and for each k-mer in the sketch,
 * the transcript ID and a pointer to the sketch are stored.
 *
 * @param transcript_sketches An unordered_map mapping transcript IDs to their corresponding MultiKmerSketch.
 * @return std::unordered_map<unsigned, TranscriptMapping> A mapping from k-mer lengths to TranscriptMapping.
 */
std::unordered_map<unsigned, TranscriptMapping> build_kmer_to_transcript_map(const std::unordered_map<std::string, MultiKmerSketch>& transcript_sketches);

#endif // SKETCHING_H
