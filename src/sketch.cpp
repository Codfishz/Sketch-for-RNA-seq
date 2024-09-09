#include "sketch.h"
#include <vector>
#include <algorithm> 
#include <limits>   
#include <unordered_set>
#include <unordered_map>
#include <stdexcept>

// Bottom-k Sketch
std::vector<uint32_t> createSketch_bottomk(const std::vector<uint32_t>& hashed_kmers, int k) {

    std::vector<uint32_t> min_k_hashes;
    min_k_hashes.reserve(k);

    for (const auto& hash : hashed_kmers) {
        if (min_k_hashes.size() < k) {
            min_k_hashes.push_back(hash);
        } else if (hash < *std::max_element(min_k_hashes.begin(), min_k_hashes.end())) {
            auto max_it = std::max_element(min_k_hashes.begin(), min_k_hashes.end());
            *max_it = hash;
        }
    }

    std::make_heap(min_k_hashes.begin(), min_k_hashes.end());
    std::sort_heap(min_k_hashes.begin(), min_k_hashes.end());
    return min_k_hashes;
}
// std::vector<uint32_t> createSketch_Minhash(const std::vector<std::string>& kmers, int num_hashes) {
//     // 使用先前的函数生成每个k-mer的多个哈希值
//     std::vector<std::vector<uint32_t>> hashes = mutil_hash_kmers(kmers, num_hashes);

//     // 存储每个哈希函数的最小值
//     std::vector<uint32_t> min_hashes(num_hashes, std::numeric_limits<uint32_t>::max());

//     // 遍历所有哈希值，更新每个哈希函数的最小值
//     for (const auto& kmer_hashes : hashes) {
//         for (int i = 0; i < num_hashes; ++i) {
//             if (kmer_hashes[i] < min_hashes[i]) {
//                 min_hashes[i] = kmer_hashes[i];
//             }
//         }
//     }

//     return min_hashes;
// }
std::vector<uint32_t> createSketch_Minimizer(const std::vector<uint32_t>& hashed_kmers, int window_size) {
    if (window_size <= 0) {
        throw std::invalid_argument("Window size must be greater than 0");
    }
    std::vector<uint32_t> minimizers;

    // 遍历所有可能的窗口
    for (size_t i = 0; i <= hashed_kmers.size() - window_size; ++i) {
        // 找到窗口中的最小哈希值
        uint32_t min_hash = std::numeric_limits<uint32_t>::max();
        for (size_t j = i; j < i + window_size; ++j) {
            if (hashed_kmers[j] < min_hash) {
                min_hash = hashed_kmers[j];
            }
        }
        // 只添加新minimizer（避免重复）
        if (minimizers.empty() || minimizers.back() != min_hash) {
            minimizers.push_back(min_hash);
        }
    }

    return minimizers;
}

// 创建FracMinhash Sketch的函数
// std::vector<uint32_t> createSketch_FracMinhash(const std::vector<uint32_t>& hashed_kmers, double fraction) {
//     if (fraction <= 0.0 || fraction > 1.0) {
//         throw std::invalid_argument("Fraction must be between 0 and 1");
//     }

//     // 确定H的值为uint32_t的最大值
//     const uint32_t H = std::numeric_limits<uint32_t>::max();
//     const uint32_t threshold = static_cast<uint32_t>(H * fraction);

//     std::vector<uint32_t> fracminhash;

//     // 选择小于等于H/s的哈希值
//     for (const auto& hash : hashed_kmers) {
//         if (hash <= threshold) {
//             fracminhash.push_back(hash);
//         }
//     }

//     return fracminhash;
// }
std::unordered_set<uint32_t> createSketch_FracMinhash(const std::unordered_set<uint32_t>& hashed_kmers, double fraction) {
    if (fraction <= 0.0 || fraction > 1.0) {
        throw std::invalid_argument("Fraction must be between 0 and 1");
    }

    const uint32_t H = std::numeric_limits<uint32_t>::max();
    const uint32_t threshold = static_cast<uint32_t>(H * fraction);

    std::unordered_set<uint32_t> fracminhash;

    for (const auto& hash : hashed_kmers) {
        if (hash <= threshold) {
            fracminhash.insert(hash);  // 插入符合条件的哈希值
        }
    }

    return fracminhash;
}

std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>> build_kmer_to_transcript_map(
    const std::unordered_map<std::string, std::unordered_set<uint32_t>>& transcript_sketches) {

    std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>> kmer_to_transcripts;

    for (const auto& [transcript_id, sketch] : transcript_sketches) {
        for (const auto& kmer : sketch) {
            // 存储 transcript_id 和指向 sketch 的指针
            kmer_to_transcripts[kmer].emplace_back(transcript_id, &sketch);
        }
    }

    return kmer_to_transcripts;
}