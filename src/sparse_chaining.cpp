#include "sparse_chaining.h"
#include "kmer.h"
#include "sketch.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <fstream>
// Function to perform sparse chaining
std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>
sparse_chain(
    const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
    // 新的 mapping，key 为 kmer length
    const std::unordered_map<unsigned, 
           std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts,
    const std::unordered_map<std::string, Read>& reads,
    const std::vector<unsigned>& kmer_lengths, 
    double fraction)
{
    std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> homologous_segments;

    // 遍历每个 read
    for (const auto& [read_id, multi_sketch] : read_sketches) {
        // 为每个 transcript，在不同 kmer length 上分别统计匹配数。
        // 使用 map: transcript_id -> vector<int>，其中 vector 的下标对应 kmer_lengths 中的顺序
        std::unordered_map<std::string, std::vector<int>> match_counts;
        // 初始化 match_counts 时不必预先分配，后面遇到时再填充。

        // 针对每个 kmer length，依次统计匹配数
        for (size_t i = 0; i < kmer_lengths.size(); i++) {
            unsigned k = kmer_lengths[i];
            // 从全局 mapping 中获取当前 k 的 mapping
            auto map_it = kmer_to_transcripts.find(k);
            if (map_it == kmer_to_transcripts.end())
                continue;
            const auto& current_mapping = map_it->second;
            // 从 read 的 multi_sketch 中获取对应 k 的 sketch
            auto sketch_it = multi_sketch.sketches.find(k);
            if (sketch_it == multi_sketch.sketches.end())
                continue;
            const auto& sketch = sketch_it->second;
            // 遍历当前 sketch 中的每个 kmer
            for (const auto& kmer : sketch) {
                auto it = current_mapping.find(kmer);
                if (it != current_mapping.end()) {
                    for (const auto& [transcript_id, transcript_sketch_ptr] : it->second) {
                        // 若 transcript_id 不在 match_counts 中，初始化一个与 kmer_lengths.size() 等长的 vector
                        if (match_counts.find(transcript_id) == match_counts.end()) {
                            match_counts[transcript_id] = std::vector<int>(kmer_lengths.size(), 0);
                        }
                        match_counts[transcript_id][i]++;  // 在第 i 个 kmer length 上计数
                    }
                }
            }
        }

        // 对于每个 transcript，分别计算每个 kmer length 的最大匹配数（仅在当前 read 内部）
        std::vector<int> max_counts(kmer_lengths.size(), 0);
        for (const auto& [transcript_id, counts_vec] : match_counts) {
            for (size_t i = 0; i < counts_vec.size(); i++) {
                if (counts_vec[i] > max_counts[i])
                    max_counts[i] = counts_vec[i];
            }
        }
        // 计算每个 kmer length 的阈值
        std::vector<double> thresholds(kmer_lengths.size(), 0.0);
        for (size_t i = 0; i < max_counts.size(); i++) {
            thresholds[i] = fraction * max_counts[i];
        }

        // 筛选出在所有 kmer length 上均达到阈值的 transcripts
        std::vector<std::pair<std::string, int>> candidate_transcripts;
        for (const auto& [transcript_id, counts_vec] : match_counts) {
            bool meets_all = true;
            int final_score = 0;
            for (size_t i = 0; i < counts_vec.size(); i++) {
                if (counts_vec[i] < thresholds[i]) {
                    meets_all = false;
                    break;
                }
                // 这里可以根据需要选择累加、取最小值或其它方式作为最终得分
                final_score += counts_vec[i];
            }
            if (meets_all) {
                candidate_transcripts.emplace_back(transcript_id, final_score);
            }
        }

        // 按得分降序排序候选 transcript
        std::sort(candidate_transcripts.begin(), candidate_transcripts.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        homologous_segments[read_id] = candidate_transcripts;
    }

    return homologous_segments;
}

// std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> sparse_chain(
//     const std::unordered_map<std::string, MultiKmerSketch>& read_sketches,
//     const std::pair<
//         std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>,
//         std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>
//     >& kmer_to_transcripts_pair,
//     const std::unordered_map<std::string, Transcript>& transcripts,
//     const std::unordered_map<std::string, Read>& reads,
//     int kmer_length, double fraction)
// {
//     // 分别取出两个 mapping：一个对应 sketch1，一个对应 sketch2
//     const auto& kmer_to_transcripts1 = kmer_to_transcripts_pair.first;
//     const auto& kmer_to_transcripts2 = kmer_to_transcripts_pair.second;

//     std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> homologous_segments;

//     // 针对每个 read
//     for (const auto& [read_id, multi_sketch] : read_sketches) {
//         std::unordered_map<std::string, int> match_counts1; // kmer1 匹配计数
//         std::unordered_map<std::string, int> match_counts2; // kmer2 匹配计数

//         // 对 read 的 sketch1 中的 kmer 进行匹配
//         for (const auto& kmer : multi_sketch.sketch1) {
//             auto it = kmer_to_transcripts1.find(kmer);
//             if (it != kmer_to_transcripts1.end()) {
//                 for (const auto& [transcript_id, transcript_sketch_ptr] : it->second) {
//                     match_counts1[transcript_id]++;
//                 }
//             }
//         }

//         // 对 read 的 sketch2 中的 kmer 进行匹配
//         for (const auto& kmer : multi_sketch.sketch2) {
//             auto it = kmer_to_transcripts2.find(kmer);
//             if (it != kmer_to_transcripts2.end()) {
//                 for (const auto& [transcript_id, transcript_sketch_ptr] : it->second) {
//                     match_counts2[transcript_id]++;
//                 }
//             }
//         }

//         // 分别计算两侧的最大匹配数
//         int max_count1 = 0;
//         for (const auto& kv : match_counts1) {
//             if (kv.second > max_count1) {
//                 max_count1 = kv.second;
//             }
//         }
//         int max_count2 = 0;
//         for (const auto& kv : match_counts2) {
//             if (kv.second > max_count2) {
//                 max_count2 = kv.second;
//             }
//         }

//         // 计算两侧的阈值
//         double threshold1 = 0.9 * max_count1;
//         double threshold2 = 0.9 * max_count2;

//         std::vector<std::pair<std::string, int>> candidate_transcripts;
//         // 只考虑同时在两侧都有匹配的 transcript
//         for (const auto& [transcript_id, count1] : match_counts1) {
//             auto it2 = match_counts2.find(transcript_id);
//             if (it2 != match_counts2.end()) {
//                 int count2 = it2->second;
//                 // 要求两侧匹配数均达到各自最大值的 90%
//                 if (count1 >= threshold1 && count2 >= threshold2) {
//                     // 可以用 min(count1, count2) 或其它方式作为最终得分
//                     int final_score = std::min(count1, count2);
//                     candidate_transcripts.emplace_back(transcript_id, final_score);
//                 }
//             }
//         }

//         // 对候选 transcript 根据得分进行降序排序
//         std::sort(candidate_transcripts.begin(), candidate_transcripts.end(),
//                   [](const auto& a, const auto& b) { return a.second > b.second; });

//         homologous_segments[read_id] = candidate_transcripts;
//     }

//     return homologous_segments;
// }

// std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> sparse_chain(
//     const std::unordered_map<std::string, std::unordered_set<uint32_t>>& read_sketches,
//     const std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>& kmer_to_transcripts,
//     const std::unordered_map<std::string, Transcript>& transcripts,
//     const std::unordered_map<std::string, Read>& reads,
//     int kmer_length, double fraction) {

//     std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> homologous_segments;

//     for (const auto& [read_id, read_sketch] : read_sketches) {
//         std::unordered_map<std::string, int> transcript_match_counts;

//         for (const auto& kmer : read_sketch) {
//             if (kmer_to_transcripts.find(kmer) != kmer_to_transcripts.end()) {
//                 for (const auto& [transcript_id, transcript_sketch] : kmer_to_transcripts.at(kmer)) {
//                     transcript_match_counts[transcript_id]++;
//                 }
//             }
//         }

//         // 将 match_count 结果存储为候选列表
//         std::vector<std::pair<std::string, int>> candidate_transcripts;
//         for (const auto& [transcript_id, match_count] : transcript_match_counts) {
//             if (transcripts.find(transcript_id) != transcripts.end()) {
//                 candidate_transcripts.emplace_back(transcript_id, match_count);
//             }
//         }

//         homologous_segments[read_id] = candidate_transcripts;
//     }

//     return homologous_segments;
// }



// Function to find the best matching transcripts using Ordered MinHash for candidates
// std::vector<std::string> find_best_match_orderedminhash(const std::string& read_id,
//                                                         const std::unordered_map<std::string, Read>& reads,
//                                                         const std::vector<std::string>& candidate_transcripts,
//                                                         const std::unordered_map<std::string, Transcript>& transcripts,
//                                                         int kmer_length, double fraction) {
//     // Get the read sequence
//     auto read_it = reads.find(read_id);
//     if (read_it == reads.end()) {
//         std::cerr << "Error: Read ID not found: " << read_id << std::endl;
//         return {};
//     }
//     const auto& read_sequence = read_it->second.sequence;
//     if (read_sequence.length() < kmer_length) {
//         std::cerr << "Error: Read length is shorter than kmer length" << std::endl;
//         return {};
//     }
//     // Create Ordered MinHash sketch for the read
//     auto read_hashed_kmers = extract_and_hash_kmers_with_positions_murmur(read_sequence, kmer_length, 12345);
//     if (read_hashed_kmers.empty()) {
//         std::cerr << "Error: No hashed kmers generated for read ID: " << read_id << std::endl;
//         return {};
//     }
    
//     auto read_ordered_sketch = createSketch_OrderedMinhash(read_hashed_kmers, fraction);

//     std::string target_transcript_id = "ENST00000692674.1|ENSG00000289007.2|-|-|ENST00000692674|ENSG00000289007|1520|lncRNA|";
//     auto target_it = transcripts.find(target_transcript_id);
//     if (target_it == transcripts.end()) {
//         std::cerr << "Error: Target transcript ID not found: " << target_transcript_id << std::endl;
//         return {};
//     }
//     //auto target_transcript_seq = target_it->second.sequence;
//     double target_score;

//     std::vector<std::string> best_transcript_ids;
//     double best_similarity_score = 0.0;

//     for (const auto& candidate_transcript_id : candidate_transcripts) {
//         auto transcript_it = transcripts.find(candidate_transcript_id);
//         if (transcript_it == transcripts.end()) {
//             std::cerr << "Warning: Transcript ID not found: " << candidate_transcript_id << std::endl;
//             continue;
//         }
//         const auto& transcript_sequence = transcript_it->second.sequence;
//         if (transcript_sequence.length() < kmer_length) {
//             std::cerr << "Warning: Transcript length is shorter than kmer length for transcript ID: " << candidate_transcript_id << std::endl;
//             continue;
//         }
//         auto transcript_hashed_kmers = extract_and_hash_kmers_with_positions_murmur(transcript_sequence, kmer_length, 12345);
//         if (transcript_hashed_kmers.empty()) {
//             std::cerr << "Warning: No hashed kmers generated for transcript ID: " << candidate_transcript_id << std::endl;
//             continue;
//         }
//         auto transcript_ordered_sketch = createSketch_OrderedMinhash(transcript_hashed_kmers, fraction);
//         double similarity_score = compare_relative_positions(read_ordered_sketch, transcript_ordered_sketch);
//         // if (candidate_transcript_id == target_transcript_id) {
//         //     target_score = similarity_score;
//         // }
//         if (similarity_score > best_similarity_score) {
//             best_similarity_score = similarity_score;
//             best_transcript_ids.clear();
//             best_transcript_ids.push_back(candidate_transcript_id);
//         } else if (similarity_score == best_similarity_score) {
//             best_transcript_ids.push_back(candidate_transcript_id);
//         }
        
//     }

//     return best_transcript_ids;
// }

// // Function to compare positions of common hashes between read and transcript
// double compare_relative_positions(const std::vector<std::pair<uint32_t, size_t>>& read_minhash,
//                                   const std::vector<std::pair<uint32_t, size_t>>& transcript_minhash) {
//     // Create a map from hash to position for both read and transcript
//     std::unordered_map<uint32_t, size_t> read_positions;
//     for (const auto& [hash, position] : read_minhash) {
//         read_positions[hash] = position;
//     }

//     std::unordered_map<uint32_t, size_t> transcript_positions;
//     for (const auto& [hash, position] : transcript_minhash) {
//         transcript_positions[hash] = position;
//     }

//     // Extract common hashes
//     std::vector<uint32_t> common_hashes;
//     for (const auto& [hash, position] : read_positions) {
//         if (transcript_positions.find(hash) != transcript_positions.end()) {
//             common_hashes.push_back(hash);
//         }
//     }

//     size_t total_common = common_hashes.size();
//     if (total_common < 2) {
//         return 0.0;  // Not enough common hashes to compare relative positions
//     }

//     // Compare relative positions
//     size_t similar_order_count = 0;
//     for (size_t i = 0; i < total_common - 1; ++i) {
//         uint32_t hash1 = common_hashes[i];
//         uint32_t hash2 = common_hashes[i + 1];
//         if ((read_positions[hash1] < read_positions[hash2] && transcript_positions[hash1] < transcript_positions[hash2]) ||
//             (read_positions[hash1] > read_positions[hash2] && transcript_positions[hash1] > transcript_positions[hash2])) {
//             similar_order_count++;
//         }
//     }

//     // Return a score that represents the proportion of common hashes with similar relative order
//     return static_cast<double>(similar_order_count) / (total_common - 1);
// }


int edit_distance(const std::string& str1, const std::string& str2) {
    size_t len1 = str1.length();
    size_t len2 = str2.length();

    // Create a DP table to store results of subproblems
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));

    // Fill the DP table
    for (size_t i = 0; i <= len1; ++i) {
        for (size_t j = 0; j <= len2; ++j) {
            if (i == 0) {
                dp[i][j] = j;  // Minimum operations = j (all insertions)
            } else if (j == 0) {
                dp[i][j] = i;  // Minimum operations = i (all deletions)
            } else if (str1[i - 1] == str2[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];  // No operation needed if characters match
            } else {
                dp[i][j] = 1 + std::min({dp[i - 1][j],     // Deletion
                                         dp[i][j - 1],     // Insertion
                                         dp[i - 1][j - 1]  // Replacement
                                         });
            }
        }
    }
    return dp[len1][len2];
}