#ifndef ISOFORM_ASSIGNMENT_H
#define ISOFORM_ASSIGNMENT_H

#include <vector>
#include <string>
#include <unordered_map>
#include "data_io.h"

/**
 * @brief Estimate isoform abundances using the EM algorithm.
 *
 * This function takes the sparse chaining results (homologous segments) and uses the Expectation-Maximization (EM)
 * algorithm to estimate the relative abundance (pi) of each isoform. The process iterates up to a maximum number of
 * iterations or until the change between iterations falls below the convergence threshold.
 *
 * @param homologous_segments An unordered_map where each key is a read ID and the value is a vector of pairs.
 *        Each pair contains a transcript ID and an integer score indicating the strength of the match for that transcript.
 * @param transcripts An unordered_map mapping transcript IDs to Transcript structures.
 * @param max_iterations The maximum number of iterations for the EM algorithm.
 * @param convergence_threshold The threshold below which the algorithm is considered to have converged.
 * @return std::unordered_map<std::string, double> A map where the key is a transcript ID and the value is the estimated abundance (pi) for that isoform.
 */
std::unordered_map<std::string, double> estimate_isoform_abundance_em(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, Transcript>& transcripts,
    int max_iterations, double convergence_threshold);

/**
 * @brief Assign reads to isoforms based on estimated abundances.
 *
 * This function uses the homologous segments from sparse chaining and the isoform abundance estimates (pi) to assign
 * reads to their corresponding isoforms. It returns a mapping from transcript IDs to the number of reads assigned.
 *
 * @param homologous_segments An unordered_map where each key is a read ID and the value is a vector of pairs.
 *        Each pair contains a transcript ID and an integer score from the sparse chaining process.
 * @param pi An unordered_map mapping transcript IDs to their estimated abundances from the EM algorithm.
 * @param transcripts An unordered_map mapping transcript IDs to Transcript structures.
 * @return std::unordered_map<std::string, double> A map where the key is a transcript ID and the value is the number of reads assigned to that isoform.
 */
std::unordered_map<std::string, double> assign_reads_to_isoforms(
    const std::unordered_map<std::string, std::vector<std::pair<std::string, int>>>& homologous_segments,
    const std::unordered_map<std::string, double>& pi,
    const std::unordered_map<std::string, Transcript>& transcripts);

#endif // ISOFORM_ASSIGNMENT_H
