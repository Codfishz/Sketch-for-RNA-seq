#ifndef KMER_H
#define KMER_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>
#include "MurmurHash3.h"



// std::unordered_map<int, std::unordered_set<uint32_t>> extract_and_hash_kmers_multiple(
//     const std::string &sequence,
//     const std::vector<int> &k_values,
//     uint32_t base,
//     uint32_t mod);

// std::unordered_set<uint32_t> extract_and_hash_kmers(const std::string& sequence, int k);
// bool is_valid_kmer(const std::string& kmer);
// std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions(const std::string& sequence, int kmer_length); 
std::unordered_set<uint32_t> extract_and_hash_kmers_murmur(const std::string& sequence, int k, uint32_t seed);
std::unordered_set<uint32_t> extract_and_hash_kmers_nthash(const std::string &sequence, int k);
// std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions_murmur(const std::string& sequence, int kmer_length, uint32_t seed);
#endif // KMER_H
