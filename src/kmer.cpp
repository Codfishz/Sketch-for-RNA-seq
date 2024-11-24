#include "kmer.h"
#include <stdexcept>
#include <unordered_set>
//#include <iostream>
#include "MurmurHash3.h"


std::unordered_set<uint32_t> extract_and_hash_kmers(const std::string& sequence, int k) {
    if (sequence.size() < k) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    std::unordered_set<uint32_t> hash_set;

    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        uint32_t hash_value = std::hash<std::string>{}(kmer);
        hash_set.insert(hash_value); 
    }

    return hash_set;
}

// Function to extract and hash kmers with positions
std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions(const std::string& sequence, int kmer_length) {
    std::vector<std::pair<uint32_t, size_t>> hashed_kmers;
    std::hash<std::string> hash_fn;
    if (sequence.size() < kmer_length) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    for (size_t i = 0; i <= sequence.length() - kmer_length; ++i) {
        std::string kmer = sequence.substr(i, kmer_length);
        uint32_t hashed_value = hash_fn(kmer);
        hashed_kmers.emplace_back(hashed_value, i);
    }

    return hashed_kmers;
}


std::unordered_set<uint32_t> extract_and_hash_kmers_murmur(const std::string& sequence, int k, uint32_t seed) {
    if (sequence.size() < k) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    std::unordered_set<uint32_t> hash_set;

    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        uint32_t hash_value;
        MurmurHash3_x86_32(kmer.c_str(), kmer.length(), seed, &hash_value);  
        hash_set.insert(hash_value);
    }

    return hash_set;
}

std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions_murmur(const std::string& sequence, int kmer_length, uint32_t seed) {
    std::vector<std::pair<uint32_t, size_t>> hashed_kmers;

    if (sequence.size() < kmer_length) {
        throw std::runtime_error("Sequence length is shorter than kmer_length");
    }

    for (size_t i = 0; i <= sequence.length() - kmer_length; ++i) {
        std::string kmer = sequence.substr(i, kmer_length);
        uint32_t hashed_value;
        MurmurHash3_x86_32(kmer.c_str(), kmer.length(), seed, &hashed_value); 
        hashed_kmers.emplace_back(hashed_value, i);
    }

    return hashed_kmers;
}


// bool is_valid_kmer(const std::string& kmer) {
//     for (char c : kmer) {
//         if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
//             return false;
//         }
//     }
//     return true;
// }
