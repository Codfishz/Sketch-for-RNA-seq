#ifndef SPARSE_CHAINING_H
#define SPARSE_CHAINING_H

#include <vector>
#include <utility>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include "data_io.h"
#include "kmer.h"
#include "sketch.h"

// 稀疏链比对函数声明

// std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> sparse_chain(
//     const std::unordered_map<std::string, std::unordered_set<uint32_t>>& read_sketches,
//     const std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>& kmer_to_transcripts,
//     const std::unordered_map<std::string, Transcript>& transcripts,
//     const std::unordered_map<std::string, Read>& reads,
//     int kmer_length, double fraction) ;
std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>
sparse_chain(
    const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
    // 新的 mapping，key 为 kmer length
    const std::unordered_map<unsigned, 
           std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts,
    const std::unordered_map<std::string, Read>& reads,
    const std::vector<unsigned>& kmer_lengths, 
    double fraction);
// std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> sparse_chain(
//     const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
//     const std::pair<
//         std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>,
//         std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>
//     >& kmer_to_transcripts_pair,
//     const std::unordered_map<std::string, Transcript>& transcripts,
//     const std::unordered_map<std::string, Read>& reads,
//     int kmer_length, double fraction);
// std::vector<std::string> find_best_match_orderedminhash(const std::string& read_id,
//                                                         const std::unordered_map<std::string, Read>& reads,
//                                                         const std::vector<std::string>& candidate_transcripts,
//                                                         const std::unordered_map<std::string, Transcript>& transcripts,
//                                                         int kmer_length, double fraction);
// double compare_relative_positions(const std::vector<std::pair<uint32_t, size_t>>& read_minhash,
//                                   const std::vector<std::pair<uint32_t, size_t>>& transcript_minhash);
int edit_distance(const std::string& str1, const std::string& str2);
#endif // SPARSE_CHAINING_H
