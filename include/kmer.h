#ifndef KMER_H
#define KMER_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>


std::unordered_set<uint32_t> extract_and_hash_kmers(const std::string& sequence, int k);
bool is_valid_kmer(const std::string& kmer);

#endif // KMER_H
