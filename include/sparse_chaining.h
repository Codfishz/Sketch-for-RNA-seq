#ifndef SPARSE_CHAINING_H
#define SPARSE_CHAINING_H

#include <vector>
#include <utility>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "data_io.h"
#include "kmer.h"
#include "sketch.h"

/**
 * @brief Perform sparse chaining for transcript-read matching.
 *
 * This function compares k-mer sketches of reads against precomputed transcript sketches to
 * identify homologous transcript segments for each read. For each read, the function counts
 * matching k-mers for each transcript across multiple k-mer lengths, applies a fractional
 * threshold to filter out low-confidence matches, and returns a sorted list of candidate transcripts
 * along with a final score.
 *
 * @param read_sketches An unordered_map mapping read IDs to their MultiKmerSketch.
 * @param kmer_to_transcripts A mapping for each k-mer length. The outer key is the k-mer length,
 *        and the inner map maps a k-mer hash value (uint32_t) to a vector of pairs. Each pair contains
 *        a transcript ID (std::string) and a pointer to the corresponding sketch (std::unordered_set<uint32_t>*).
 * @param transcripts An unordered_map mapping transcript IDs to Transcript structures.
 * @param kmer_lengths A vector of k-mer lengths used to generate the sketches.
 * @param fraction The fraction used to compute thresholds. For each k-mer length, only transcripts
 *        with a match count above (fraction * max_count) in the current read are retained.
 * @return std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> A mapping from read IDs
 *         to a vector of candidate transcripts. Each candidate is represented as a pair, where the first
 *         element is the transcript ID and the second element is the final score.
 */
std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> sparse_chain(
    const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
    // New mapping: key is k-mer length.
    const std::unordered_map<unsigned, 
           std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts,
    const std::vector<unsigned>& kmer_lengths, 
    double fraction);

#endif // SPARSE_CHAINING_H
