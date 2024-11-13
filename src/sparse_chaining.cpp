#include "sparse_chaining.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <string>
#include <utility>

// 稀疏链比对函数实现

std::unordered_map<std::string, std::string> sparse_chain(
    const std::unordered_map<std::string, std::unordered_set<uint32_t>>& read_sketches,
    const std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, std::string> homologous_segments;

    for (const auto& [read_id, read_sketch] : read_sketches) {
        std::unordered_map<std::string, int> transcript_match_counts;

        // 遍历 read 的每个 k-mer
        for (const auto& kmer : read_sketch) {
            // 查找包含该 k-mer 的所有转录本 sketch
            if (kmer_to_transcripts.find(kmer) != kmer_to_transcripts.end()) {
                for (const auto& [transcript_id, transcript_sketch] : kmer_to_transcripts.at(kmer)) {
                    transcript_match_counts[transcript_id]++;
                }
            }
        }

        // 找到与当前 read 匹配度最高的转录本
        std::string best_transcript_id;
        int max_match_count = 0;
        int same_max_match = 0;
        for (const auto& [transcript_id, match_count] : transcript_match_counts) {
            // std::cout<<(transcript_id)<<std::endl;
            // if (transcript_id=="ENST00000692674.1|ENSG00000289007.2|-|-|ENST00000692674|ENSG00000289007|1520|lncRNA|"){
            //     std::cout<<"correct match:"<<read_id<<":"<<match_count<<std::endl;
            // }
            if (match_count == max_match_count){
                same_max_match+=1;
            }
            if (match_count > max_match_count) {
                same_max_match = 0;
                max_match_count = match_count;
                best_transcript_id = transcript_id;
            }
        }

        if (!best_transcript_id.empty()) {
            // std::cout << same_max_match <<std::endl;
            //std::cout << read_id << " -> " << best_transcript_id << " "<< max_match_count <<std::endl;
            homologous_segments[read_id] = best_transcript_id;
        }
    }

    return homologous_segments;
}
