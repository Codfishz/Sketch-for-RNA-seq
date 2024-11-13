#ifndef ISOFORM_ASSIGNMENT_H
#define ISOFORM_ASSIGNMENT_H

#include <vector>
#include <string>
#include <unordered_map>
#include "data_io.h"

// 从稀疏链比对结果中分配 reads 到 isoforms 的函数声明


std::unordered_map<std::string, int> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::string>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts);


std::unordered_map<std::string, double> calculate_tpm(
    const std::unordered_map<std::string, int>& read_counts,
    const std::unordered_map<std::string, Transcript>& transcripts);

#endif // ISOFORM_ASSIGNMENT_H
