#ifndef ISOFORM_ASSIGNMENT_H
#define ISOFORM_ASSIGNMENT_H

#include <vector>
#include <string>
#include <unordered_map>
#include "data_io.h"

// 从稀疏链比对结果中分配 reads 到 isoforms 的函数声明


std::unordered_map<std::string, double> estimate_isoform_abundance_em(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts,
    int max_iterations, double convergence_threshold);

// void update_pi_with_newton(
//     std::unordered_map<std::string, double>& pi,
//     const std::unordered_map<std::string, double>& posterior_sums,
//     double regularization_lambda,
//     double step_size);

std::unordered_map<std::string, double> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, double>& pi,
    const std::unordered_map<std::string, Transcript>& transcripts);


std::unordered_map<std::string, double> calculate_tpm(
    const std::unordered_map<std::string, double>& read_counts, // 加权后的 read_counts
    const std::unordered_map<std::string, Transcript>& transcripts);

#endif // ISOFORM_ASSIGNMENT_H
