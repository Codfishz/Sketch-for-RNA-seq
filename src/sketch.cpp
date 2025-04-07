#include "sketch.h"
#include <vector>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <stdexcept>
#include <nthash/nthash.hpp>

/**
 * @brief Create a fractional MinHash sketch for a given DNA sequence.
 *
 * This function extracts k-mers from the input sequence using the nthash library and
 * selects a fraction of the k-mers based on the provided fraction value. Only those k-mers
 * whose hash value is less than or equal to the threshold (calculated as the maximum uint32_t value multiplied by fraction)
 * are inserted into the sketch.
 *
 * @param sequence The input DNA sequence.
 * @param k The length of the k-mers.
 * @param fraction The fraction threshold for selecting k-mers (e.g., 0.05 means keeping 5%).
 * @return std::unordered_set<uint32_t> A set containing the hash values of selected k-mers.
 *
 * @throws std::runtime_error if the sequence length is shorter than k.
 */
std::unordered_set<uint32_t> createSketch_FracMinhash_direct(const std::string &sequence, int k, double fraction) {
    const uint32_t H = std::numeric_limits<uint32_t>::max();
    const uint32_t threshold = static_cast<uint32_t>(H * fraction);
    std::unordered_set<uint32_t> sketch;
    // Reserve capacity to reduce rehash overhead; note: sequence.size() - k + 1 is the maximum number of k-mers.
    sketch.reserve(sequence.size() - k + 1);

    nthash::NtHash nth(sequence, 1, k);
    while (nth.roll()) {
        uint32_t hash_value = nth.get_forward_hash();
        if (hash_value <= threshold) {
            sketch.insert(hash_value);
        }
    }
    return sketch;
}

/**
 * @brief Build a mapping from k-mers to transcripts based on precomputed sketches.
 *
 * This function iterates over a collection of transcript sketches (each associated with a transcript ID)
 * and constructs a mapping from k-mer hash values to transcript information. The outer key is the k-mer length.
 * For each k-mer in the sketch, the transcript ID and a pointer to the corresponding sketch are added to the mapping.
 *
 * @param transcript_sketches An unordered_map where the key is the transcript ID and the value is a MultiKmerSketch.
 * @return std::unordered_map<unsigned, TranscriptMapping> A mapping from k-mer lengths to their corresponding TranscriptMapping.
 */
std::unordered_map<unsigned, TranscriptMapping> build_kmer_to_transcript_map(const std::unordered_map<std::string, MultiKmerSketch>& transcript_sketches) {
    // The outer key of the returned map is the k-mer length.
    std::unordered_map<unsigned, TranscriptMapping> kmer_to_transcripts;

    // Iterate over each transcript and its corresponding MultiKmerSketch.
    for (const auto& [transcript_id, mks] : transcript_sketches) {
        // For each sketch corresponding to a specific k-mer length.
        for (const auto& kv : mks.sketches) {
            unsigned k = kv.first;
            // Obtain the sketch (a set of k-mer hash values) for the current k.
            const SketchType& sketch = kv.second;

            // Ensure there is a mapping entry for the current k-mer length.
            TranscriptMapping& mapping = kmer_to_transcripts[k];

            // For each k-mer in the sketch, insert the transcript information.
            for (const auto& kmer : sketch) {
                mapping[kmer].emplace_back(transcript_id, &sketch);
            }
        }
    }

    return kmer_to_transcripts;
}
