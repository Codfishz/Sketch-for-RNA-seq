#include "kmer.h"
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include "MurmurHash3.h"
#include <nthash/nthash.hpp>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <string>


inline uint32_t nucleotide_to_int(char c) {
    switch(c) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        default: return 0;
    }
}

// std::unordered_map<int, std::unordered_set<uint32_t>> extract_and_hash_kmers_multiple(
//     const std::string &sequence,
//     const std::vector<int> &k_values,
//     uint32_t base = 101,
//     uint32_t mod  = 1000000007)
// {
//     std::unordered_map<int, std::unordered_set<uint32_t>> result;
//     // 对于每个 k，初始化结果集合；如果序列长度小于 k，则报错
//     for (int k : k_values) {
//         if (sequence.size() < static_cast<size_t>(k)) {
//             throw std::runtime_error("Sequence length is shorter than one of the k-mer lengths");
//         }
//         result[k] = std::unordered_set<uint32_t>();
//     }
    
//     // currentHashes 用于保存每个 k 对应的当前滚动哈希值
//     std::vector<uint32_t> currentHashes(k_values.size(), 0);
    
//     // powers 保存每个 k 的最高幂 base^(k-1) mod mod，用于在滑动窗口中减去最左侧字符的贡献
//     std::vector<uint32_t> powers(k_values.size(), 1);
//     for (size_t j = 0; j < k_values.size(); j++) {
//         int k = k_values[j];
//         uint32_t p = 1;
//         for (int i = 1; i < k; i++) {
//             p = (p * base) % mod;
//         }
//         powers[j] = p;
//     }
    
//     // 一次遍历序列，更新每个 k 值的滚动哈希
//     for (size_t i = 0; i < sequence.size(); i++) {
//         // 使用 nucleotide_to_int 映射当前字符
//         uint32_t currentVal = nucleotide_to_int(sequence[i]);
//         // 对于每个 k 值更新
//         for (size_t j = 0; j < k_values.size(); j++) {
//             int k = k_values[j];
//             // 当窗口已满时，先减去窗口最左侧字符的贡献
//             if (i >= static_cast<size_t>(k)) {
//                 uint32_t outgoingVal = nucleotide_to_int(sequence[i - k]);
//                 currentHashes[j] = (mod + currentHashes[j] - (outgoingVal * powers[j]) % mod) % mod;
//             }
//             // 更新当前哈希：乘以基数加上当前字符的值，再取 mod
//             currentHashes[j] = (currentHashes[j] * base + currentVal) % mod;
            
//             // 当窗口长度达到 k 时，将当前哈希值保存到结果集合中
//             if (i >= static_cast<size_t>(k - 1)) {
//                 result[k].insert(currentHashes[j]);
//             }
//         }
//     }
    
//     return result;
// }

// std::unordered_set<uint32_t> extract_and_hash_kmers(const std::string& sequence, int k) {
//     if (sequence.size() < k) {
//         throw std::runtime_error("Sequence length is shorter than k-mer length");
//     }

//     std::unordered_set<uint32_t> hash_set;

//     for (size_t i = 0; i <= sequence.size() - k; ++i) {
//         std::string kmer = sequence.substr(i, k);
//         uint32_t hash_value = std::hash<std::string>{}(kmer);
//         hash_set.insert(hash_value); 
//     }

//     return hash_set;
// }

// // Function to extract and hash kmers with positions
// std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions(const std::string& sequence, int kmer_length) {
//     std::vector<std::pair<uint32_t, size_t>> hashed_kmers;
//     std::hash<std::string> hash_fn;
//     if (sequence.size() < kmer_length) {
//         throw std::runtime_error("Sequence length is shorter than k-mer length");
//     }

//     for (size_t i = 0; i <= sequence.length() - kmer_length; ++i) {
//         std::string kmer = sequence.substr(i, kmer_length);
//         uint32_t hashed_value = hash_fn(kmer);
//         hashed_kmers.emplace_back(hashed_value, i);
//     }

//     return hashed_kmers;
// }


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

std::unordered_set<uint32_t> extract_and_hash_kmers_nthash(const std::string &sequence, int k) {
    if (sequence.size() < static_cast<size_t>(k)) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    std::unordered_set<uint32_t> hash_set;
    nthash::NtHash nth(sequence, 1, k);

    while (nth.roll()) {
        uint32_t hash_value = nth.get_forward_hash();
        hash_set.insert(hash_value);
    }
    
    return hash_set;
}

// std::vector<std::pair<uint32_t, size_t>> extract_and_hash_kmers_with_positions_murmur(const std::string& sequence, int kmer_length, uint32_t seed) {
//     std::vector<std::pair<uint32_t, size_t>> hashed_kmers;

//     if (sequence.size() < kmer_length) {
//         throw std::runtime_error("Sequence length is shorter than kmer_length");
//     }

//     for (size_t i = 0; i <= sequence.length() - kmer_length; ++i) {
//         std::string kmer = sequence.substr(i, kmer_length);
//         uint32_t hashed_value;
//         MurmurHash3_x86_32(kmer.c_str(), kmer.length(), seed, &hashed_value); 
//         hashed_kmers.emplace_back(hashed_value, i);
//     }

//     return hashed_kmers;
// }


// bool is_valid_kmer(const std::string& kmer) {
//     for (char c : kmer) {
//         if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
//             return false;
//         }
//     }
//     return true;
// }
