#include "isoform_assignment.h"
#include <iostream> 
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

// 从稀疏链比对结果中分配 reads 到 isoforms

std::unordered_map<std::string, int> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::string>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, int> read_counts;

    for (const auto& [read_id, transcript_id] : homologous_segments) {
        // 检查该转录本是否存在于 transcripts 中
        if (transcripts.find(transcript_id) != transcripts.end()) {
            // 累加每个转录本的 reads 数量
            read_counts[transcript_id] += 1;
        }
    }

    return read_counts;
}

// 计算TPM
std::unordered_map<std::string, double> calculate_tpm(
    const std::unordered_map<std::string, int>& read_counts,
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, double> tpm;
    double sum_rpk = 0.0;

    // 计算每个转录本的RPK（Reads Per Kilobase）
    for (const auto& [transcript_id, count] : read_counts) {
        if (transcripts.find(transcript_id) != transcripts.end()) {
            const auto& transcript = transcripts.at(transcript_id);
            double rpk = count / (transcript.length / 1000.0);
            tpm[transcript_id] = rpk;
            sum_rpk += rpk;
        }
    }

    // 标准化RPK以获得TPM
    for (auto& [transcript_id, rpk] : tpm) {
        tpm[transcript_id] = rpk / sum_rpk * 1e6;
    }

    return tpm;
}