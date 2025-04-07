#include "sparse_chaining.h"
#include "kmer.h"
#include "sketch.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>

/**
 * @brief Perform sparse chaining to identify homologous transcript segments for each read.
 *
 * This function compares k-mer sketches of reads with precomputed transcript sketches. For each read,
 * it computes matching counts for each transcript across different k-mer lengths, applies a fractional
 * threshold to determine candidate transcripts, and finally sorts these candidates by their score.
 *
 * @param read_sketches An unordered_map mapping read IDs to their MultiKmerSketch.
 * @param kmer_to_transcripts A mapping for each k-mer length that maps k-mer hash values to a vector
 *        of transcript information (transcript ID and pointer to the corresponding sketch).
 * @param transcripts An unordered_map mapping transcript IDs to Transcript structures.
 * @param kmer_lengths A vector of k-mer lengths used to generate the sketches.
 * @param fraction A fraction value used to compute a threshold for matching counts.
 *                 Only transcripts with counts above fraction*max_count for every k-mer length are retained.
 * @return An unordered_map mapping each read ID to a sorted vector of candidate transcripts.
 *         Each candidate is represented as a pair of (transcript ID, final score).
 */
std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>
sparse_chain(
    const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
    // New mapping: key is k-mer length
    const std::unordered_map<unsigned, 
           std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts,
    const std::vector<unsigned>& kmer_lengths, 
    double fraction)
{
    std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> homologous_segments;

    // Iterate over each read.
    for (const auto& [read_id, multi_sketch] : read_sketches) {
        // For each transcript, count the matches across different k-mer lengths.
        // Use a map: transcript_id -> vector<int> where each index corresponds to a k-mer length.
        std::unordered_map<std::string, std::vector<int>> match_counts;

        // For each k-mer length, count the matches.
        for (size_t i = 0; i < kmer_lengths.size(); i++) {
            unsigned k = kmer_lengths[i];
            // Get the mapping for current k-mer length from the global mapping.
            auto map_it = kmer_to_transcripts.find(k);
            if (map_it == kmer_to_transcripts.end())
                continue;
            const auto& current_mapping = map_it->second;
            // Get the sketch for the current k-mer length from the read's MultiKmerSketch.
            auto sketch_it = multi_sketch.sketches.find(k);
            if (sketch_it == multi_sketch.sketches.end())
                continue;
            const auto& sketch = sketch_it->second;
            // Iterate over each k-mer in the sketch.
            for (const auto& kmer : sketch) {
                auto it = current_mapping.find(kmer);
                if (it != current_mapping.end()) {
                    for (const auto& [transcript_id, transcript_sketch_ptr] : it->second) {
                        // If transcript_id is not in match_counts, initialize a vector of length kmer_lengths.size().
                        if (match_counts.find(transcript_id) == match_counts.end()) {
                            match_counts[transcript_id] = std::vector<int>(kmer_lengths.size(), 0);
                        }
                        match_counts[transcript_id][i]++;  // Increment count for k-mer length index i.
                    }
                }
            }
        }

        // Compute the maximum match count for each k-mer length among transcripts for the current read.
        std::vector<int> max_counts(kmer_lengths.size(), 0);
        for (const auto& [transcript_id, counts_vec] : match_counts) {
            for (size_t i = 0; i < counts_vec.size(); i++) {
                if (counts_vec[i] > max_counts[i])
                    max_counts[i] = counts_vec[i];
            }
        }
        // Compute threshold for each k-mer length.
        std::vector<double> thresholds(kmer_lengths.size(), 0.0);
        for (size_t i = 0; i < max_counts.size(); i++) {
            thresholds[i] = fraction * max_counts[i];
        }

        // Filter transcripts that meet the threshold at every k-mer length.
        std::vector<std::pair<std::string, int>> candidate_transcripts;
        for (const auto& [transcript_id, counts_vec] : match_counts) {
            bool meets_all = true;
            int final_score = 0;
            for (size_t i = 0; i < counts_vec.size(); i++) {
                if (counts_vec[i] < thresholds[i]) {
                    meets_all = false;
                    break;
                }
                // You can choose different methods for the final score; here we sum the counts.
                final_score += counts_vec[i];
            }
            if (meets_all) {
                candidate_transcripts.emplace_back(transcript_id, final_score);
            }
        }

        // Sort candidate transcripts in descending order of their final score.
        std::sort(candidate_transcripts.begin(), candidate_transcripts.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        homologous_segments[read_id] = candidate_transcripts;
    }

    return homologous_segments;
}
