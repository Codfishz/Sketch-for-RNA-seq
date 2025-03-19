#ifndef SKETCHING_H
#define SKETCHING_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstdint>

using SketchType = std::unordered_set<uint32_t>;
using TranscriptMapping = std::unordered_map<uint32_t, std::vector<std::pair<std::string, const SketchType*>>>;

struct MultiKmerSketch {
    std::unordered_map<unsigned, SketchType> sketches;
};

std::vector<uint32_t> createSketch_bottomk(const std::vector<uint32_t>& hashed_kmers, int k);

std::vector<uint32_t> createSketch_Minimizer(const std::vector<uint32_t>& hashed_kmers, int window_size);

// std::vector<uint32_t> createSketch_FracMinhash(const std::vector<uint32_t>& hashed_kmers, double fraction);
std::unordered_set<uint32_t> createSketch_FracMinhash(const std::unordered_set<uint32_t>& hashed_kmers, double fraction);
std::unordered_set<uint32_t> createSketch_FracMinhash_vector(const std::vector<uint64_t>& hashed_kmers, double fraction);

// std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>> build_kmer_to_transcript_map(
//     const std::unordered_map<std::string, std::unordered_set<uint32_t>>& transcript_sketches);
// std::pair<
//     std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>,
//     std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>
// >
// build_kmer_to_transcript_map(const std::unordered_map<std::string, MultiKmerSketch>& transcript_sketches);
std::unordered_map<unsigned, TranscriptMapping>
build_kmer_to_transcript_map(const std::unordered_map<std::string, MultiKmerSketch>& transcript_sketches);
// std::vector<uint32_t> createSketch_Minhash(const std::vector<std::string>& kmers, int num_hashes); 
std::vector<std::pair<uint32_t, size_t>> createSketch_OrderedMinhash(const std::vector<std::pair<uint32_t, size_t>>& hashed_kmers, double fraction);
#endif // SKETCHING_H