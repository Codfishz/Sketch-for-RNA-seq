#include "isoform_assignment.h"
#include <iostream> 
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <chrono>

// std::unordered_map<std::string, double> estimate_isoform_abundance_em(
//     const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
//     const std::unordered_map<std::string, Transcript>& transcripts,
//     int max_iterations = 100, double convergence_threshold = 0.05) {

//     // Initialize pi with a uniform distribution
//     std::unordered_map<std::string, double> pi;
//     for (const auto& [transcript_id, _] : transcripts) {
//         pi[transcript_id] = 1.0 / transcripts.size();
//     }

//     float previous_total_change=1;
//     // EM algorithm iterations
//    for (int iteration = 0; iteration < max_iterations; ++iteration) {
//         // E step: Calculate posterior sums
//         std::unordered_map<std::string, double> posterior_sums;
//         for (const auto& [transcript_id, _] : transcripts) {
//             posterior_sums[transcript_id] = 0.0;
//         }

//         double total_change = 0.0;

//         for (const auto& [read_id, candidates] : homologous_segments) {
//             double denominator = 0.0;
//             std::vector<double> numerators(candidates.size());

//             for (size_t i = 0; i < candidates.size(); ++i) {
//                 const auto& [transcript_id, match_count] = candidates[i];
//                 double value = pi[transcript_id] * static_cast<double>(match_count);
//                 numerators[i] = value;
//                 denominator += value;
//             }

//             if (denominator > 1e-10) { // 防止分母为零
//                 for (size_t i = 0; i < candidates.size(); ++i) {
//                     const auto& [transcript_id, _] = candidates[i];
//                     posterior_sums[transcript_id] += numerators[i] / denominator;
//                 }
//             }
//         }

//         // M step: Update pi using Newton's method
//         double regularization_lambda = 1e-6;
//         double step_size = 0.5; // 限制步长
//         update_pi_with_newton(pi, posterior_sums, regularization_lambda, step_size);

//         // Calculate total change
//         double relative_change = 0.0;
//         double pi_sum = 0.0;
//         for (const auto& [transcript_id, new_pi] : pi) {
//             double expected_value = posterior_sums[transcript_id] / homologous_segments.size();
//             total_change += std::abs(new_pi - expected_value);
//             pi_sum += new_pi;
//         }
        
//         relative_change = total_change / pi_sum;
//         if (iteration > 0 && total_change > previous_total_change) {
//             step_size *= 0.5; // 当变化增大时，减小步长
//         } else {
//             step_size = std::min(step_size * 1.1, 1.0); // 否则逐渐增大步长
//         }
//         previous_total_change= total_change;

//         // Check for convergence
//         std::cout << "Iteration " << iteration << ", Relative Change: " << relative_change << std::endl;
//         if (relative_change < convergence_threshold) {
//             std::cout << "Converged after iteration: " << iteration << std::endl;
//             break;
//         }
//     }

//     return pi;
// }


// void update_pi_with_newton(
//     std::unordered_map<std::string, double>& pi,
//     const std::unordered_map<std::string, double>& posterior_sums,
//     double regularization_lambda,
//     double step_size = 0.1) { // 加入步长限制

//     for (auto& [transcript_id, old_pi] : pi) {
//         auto it = posterior_sums.find(transcript_id);
//         if (it == posterior_sums.end()) continue;

//         double posterior_sum = it->second;

//         // 计算梯度和 Hessian
//         double gradient = posterior_sum / old_pi;
//         double hessian = -posterior_sum / (old_pi * old_pi) + regularization_lambda;

//         // 防止 Hessian 为零或负值导致不稳定
//         if (hessian >= -1e-8 && hessian <= 1e-8) {
//             hessian = (hessian < 0) ? -1e-8 : 1e-8;
//         }

//         // 牛顿更新公式
//         double delta = -gradient / hessian;

//         // 加入步长限制
//         double new_pi = old_pi + step_size * delta;

//         // 确保 \(\pi_i\) 非负
//         pi[transcript_id] = std::max(new_pi, 1e-10);
//     }
// }
std::unordered_map<std::string, double> estimate_isoform_abundance_em(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts,
    int max_iterations = 100, double convergence_threshold = 0.05) {

    // 在函数内部：
    auto start_total = std::chrono::high_resolution_clock::now();

    // Uniform distribution 初始化
    auto start_uniform = std::chrono::high_resolution_clock::now();
    // 原代码
    std::unordered_map<std::string, double> pi;
    for (const auto& [transcript_id, _] : transcripts) {
        pi[transcript_id] = 1.0 / transcripts.size();
    }
    auto end_uniform = std::chrono::high_resolution_clock::now();
    std::cout << "Uniform distribution time: " 
            << std::chrono::duration<double>(end_uniform - start_uniform).count() 
            << " seconds" << std::endl;

    // EM 主循环
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        auto start_iteration = std::chrono::high_resolution_clock::now();

        // E-step
        auto start_e_step = std::chrono::high_resolution_clock::now();
        std::unordered_map<std::string, double> posterior_sums;
        double total_change = 0.0;

        for (const auto& [read_id, candidates] : homologous_segments) {
            double denominator = 0.0;
            std::vector<double> numerators;
            for (const auto& [transcript_id, match_count] : candidates) {
                double p_r_ti = static_cast<double>(match_count);
                double value = pi[transcript_id] * p_r_ti;
                numerators.push_back(value);
                denominator += value;
            }

            if (denominator > 0.0) {
                for (size_t i = 0; i < candidates.size(); ++i) {
                    const std::string& transcript_id = candidates[i].first;
                    double posterior = numerators[i] / denominator;
                    posterior_sums[transcript_id] += posterior;
                }
            }
        }
        auto end_e_step = std::chrono::high_resolution_clock::now();
        std::cout << "E-step time: " 
                << std::chrono::duration<double>(end_e_step - start_e_step).count() 
                << " seconds" << std::endl;

        // M-step
        float pseudocount = 0.01;
        auto start_m_step = std::chrono::high_resolution_clock::now();
        size_t R = homologous_segments.size();
        for (auto& [transcript_id, old_pi] : pi) {
            double new_pi = posterior_sums[transcript_id] + pseudocount / R + pseudocount;
            total_change += std::abs(new_pi - old_pi);
            pi[transcript_id] = new_pi;
        }
        auto end_m_step = std::chrono::high_resolution_clock::now();
        std::cout << "M-step time: " 
                << std::chrono::duration<double>(end_m_step - start_m_step).count() 
                << " seconds" << std::endl;

        auto end_iteration = std::chrono::high_resolution_clock::now();
        std::cout << "Iteration " << iteration << " time: " 
                << std::chrono::duration<double>(end_iteration - start_iteration).count() 
                << " seconds" << std::endl;

        if (total_change < convergence_threshold) {
            break;
        }
    }

    auto end_total = std::chrono::high_resolution_clock::now();
    std::cout << "Total time: " 
            << std::chrono::duration<double>(end_total - start_total).count() 
            << " seconds" << std::endl;
   
    return pi;
}



std::unordered_map<std::string, double> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, double>& pi,
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, double> weighted_read_counts;

    for (const auto& [read_id, candidates] : homologous_segments) {
        double total_probability = 0.0;

        // 计算分母：所有候选 transcript 的概率加权和
        for (const auto& [transcript_id, match_count] : candidates) {
            if (pi.find(transcript_id) != pi.end()) {
                total_probability += pi.at(transcript_id) * match_count;
            }
        }

        // 计算每个 transcript 的加权分配
        for (const auto& [transcript_id, match_count] : candidates) {
            if (pi.find(transcript_id) != pi.end() && total_probability > 0.0) {
                double probability = (pi.at(transcript_id) * match_count) / total_probability;
                weighted_read_counts[transcript_id] += probability;
            }
        }
    }

    return weighted_read_counts;
}



// 计算TPM
std::unordered_map<std::string, double> calculate_tpm(
    const std::unordered_map<std::string, double>& read_counts, // 加权后的 read_counts
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, double> tpm;
    double sum_rpk = 0.0;

    // 计算每个转录本的 RPK（Reads Per Kilobase）
    for (const auto& [transcript_id, count] : read_counts) {
        if (transcripts.find(transcript_id) != transcripts.end()) {
            const auto& transcript = transcripts.at(transcript_id);

            // RPK = read_count / (transcript_length / 1000)
            double rpk = count / (transcript.length / 1000.0);
            tpm[transcript_id] = rpk;
            sum_rpk += rpk;
        }
    }

    // 归一化 RPK 到 TPM（Transcripts Per Million）
    for (auto& [transcript_id, rpk] : tpm) {
        tpm[transcript_id] = rpk / sum_rpk * 1e6; // 将 RPK 标准化到百万级别
    }

    return tpm;
}

