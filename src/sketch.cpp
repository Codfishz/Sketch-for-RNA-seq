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
//     
//     std::vector<std::vector<uint32_t>> hashes = mutil_hash_kmers(kmers, num_hashes);

//    
//     std::vector<uint32_t> min_hashes(num_hashes, std::numeric_limits<uint32_t>::max());

//   
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

    for (size_t i = 0; i <= hashed_kmers.size() - window_size; ++i) {
        uint32_t min_hash = std::numeric_limits<uint32_t>::max();
        for (size_t j = i; j < i + window_size; ++j) {
            if (hashed_kmers[j] < min_hash) {
                min_hash = hashed_kmers[j];
            }
        }
        if (minimizers.empty() || minimizers.back() != min_hash) {
            minimizers.push_back(min_hash);
        }
    }

    return minimizers;
}

// 
// std::vector<uint32_t> createSketch_FracMinhash(const std::vector<uint32_t>& hashed_kmers, double fraction) {
//     if (fraction <= 0.0 || fraction > 1.0) {
//         throw std::invalid_argument("Fraction must be between 0 and 1");
//     }

//    
//     const uint32_t H = std::numeric_limits<uint32_t>::max();
//     const uint32_t threshold = static_cast<uint32_t>(H * fraction);

//     std::vector<uint32_t> fracminhash;

//    
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
            fracminhash.insert(hash); 
        }
    }

    return fracminhash;
}

std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>> build_kmer_to_transcript_map(
    const std::unordered_map<std::string, std::unordered_set<uint32_t>>& transcript_sketches) {

    std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>> kmer_to_transcripts;

    for (const auto& [transcript_id, sketch] : transcript_sketches) {
        for (const auto& kmer : sketch) {
            kmer_to_transcripts[kmer].emplace_back(transcript_id, &sketch);
        }
    }

    return kmer_to_transcripts;
}

// Function to create an ordered MinHash sketch
std::vector<std::pair<uint32_t, size_t>> createSketch_OrderedMinhash(const std::vector<std::pair<uint32_t, size_t>>& hashed_kmers, double fraction) {
    if (fraction <= 0.0 || fraction > 1.0) {
        throw std::invalid_argument("Fraction must be between 0 and 1");
    }

    const size_t total_hashes = hashed_kmers.size();
    const size_t num_to_select = static_cast<size_t>(total_hashes * fraction);

    // Copy and sort by hash value to select minimum hashes in the original order
    std::vector<std::pair<uint32_t, size_t>> sorted_hashed_kmers = hashed_kmers;
    std::partial_sort(sorted_hashed_kmers.begin(), sorted_hashed_kmers.begin() + num_to_select, sorted_hashed_kmers.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });

    // Collect the selected hashes along with their original position
    std::vector<std::pair<uint32_t, size_t>> ordered_minhash(sorted_hashed_kmers.begin(), sorted_hashed_kmers.begin() + num_to_select);

    return ordered_minhash;
}