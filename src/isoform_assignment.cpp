#include "isoform_assignment.h"
#include <iostream> 
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <chrono>

std::unordered_map<std::string, double> estimate_isoform_abundance_em(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts,
    int max_iterations = 100, double convergence_threshold = 0.05) {

    // Inside the function:

    // Initialize with a uniform distribution
    std::unordered_map<std::string, double> pi;
    for (const auto& [transcript_id, _] : transcripts) {
        pi[transcript_id] = 1.0 / transcripts.size();
    }

    // EM main loop
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::unordered_map<std::string, double> posterior_sums;  // Accumulate posterior sums
        double total_change = 0.0;

        const double probability_threshold = 1e-6;  // Low probability pruning
        const double epsilon = 1e-10;  // Avoid division by zero

        for (const auto& [read_id, candidates] : homologous_segments) {
            double denominator = 0.0;
            std::vector<double> numerators(candidates.size()); // Pre-allocate vector for improved access speed

            // Calculate denominator (direct computation, no caching)
            size_t idx = 0;
            for (const auto& [transcript_id, match_count] : candidates) {
                double p_r_ti = static_cast<double>(match_count);
                double value = pi[transcript_id] * p_r_ti;
                numerators[idx++] = value;  // Store into vector to avoid extra lookup
                denominator += value;
            }

            if (denominator > epsilon) {
                double inv_denominator = 1.0 / denominator;
                idx = 0;
                for (const auto& [transcript_id, match_count] : candidates) {
                    double posterior = numerators[idx++] * inv_denominator;
                    posterior_sums[transcript_id] += posterior;
                }
            }
        }

        // M-step
        float pseudocount = 0.01;
        size_t R = homologous_segments.size();
        for (auto& [transcript_id, old_pi] : pi) {
            double new_pi = posterior_sums[transcript_id] + pseudocount / R + pseudocount;
            total_change += std::abs(new_pi - old_pi);
            pi[transcript_id] = new_pi;
        }

        if (total_change < convergence_threshold) {
            break;
        }
    }
   
    return pi;
}

std::unordered_map<std::string, double> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, double>& pi,
    const std::unordered_map<std::string, Transcript>& transcripts) {

    std::unordered_map<std::string, double> weighted_read_counts;

    for (const auto& [read_id, candidates] : homologous_segments) {
        double total_probability = 0.0;

        // Calculate the denominator: weighted sum of probabilities for all candidate transcripts
        for (const auto& [transcript_id, match_count] : candidates) {
            if (pi.find(transcript_id) != pi.end()) {
                total_probability += pi.at(transcript_id) * match_count;
            }
        }

        // Calculate the weighted assignment for each transcript
        for (const auto& [transcript_id, match_count] : candidates) {
            if (pi.find(transcript_id) != pi.end() && total_probability > 0.0) {
                double probability = (pi.at(transcript_id) * match_count) / total_probability;
                weighted_read_counts[transcript_id] += probability;
            }
        }
    }

    return weighted_read_counts;
}
