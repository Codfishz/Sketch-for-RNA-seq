#include "kmer.h"
#include <stdexcept>
#include <unordered_set>
//#include <iostream>
//#include "MurmurHash3.h"


// 提取k-mers并直接将哈希值插入哈希表的函数
std::unordered_set<uint32_t> extract_and_hash_kmers(const std::string& sequence, int k) {
    if (sequence.size() < k) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    std::unordered_set<uint32_t> hash_set;

    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        uint32_t hash_value = std::hash<std::string>{}(kmer);
        hash_set.insert(hash_value);  // 直接插入哈希表
    }

    return hash_set;
}



//检查k-mer是否有效的函数（例如，是否包含非碱基字符）
// bool is_valid_kmer(const std::string& kmer) {
//     for (char c : kmer) {
//         if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
//             return false;
//         }
//     }
//     return true;
// }
