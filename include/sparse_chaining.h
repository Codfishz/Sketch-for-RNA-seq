#ifndef SPARSE_CHAINING_H
#define SPARSE_CHAINING_H

#include <vector>
#include <utility>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

// 稀疏链比对函数声明

std::unordered_map<std::string, std::pair<std::string, int>> sparse_chain(
    const std::unordered_map<std::string, std::unordered_set<uint32_t>>& read_sketches,
    const std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>& kmer_to_transcripts);
#endif // SPARSE_CHAINING_H
