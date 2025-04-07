#ifndef KMER_H
#define KMER_H

#include <vector>
#include <string>
#include <cstdint>
#include <unordered_set>

/**
 * @brief Extracts and hashes k-mers from a DNA sequence using nthash.
 *
 * This function extracts all possible k-mers from the given DNA sequence and computes their
 * hash values using the nthash algorithm. The resulting hash values are stored in an unordered_set.
 *
 * @param sequence The input DNA sequence.
 * @param k The length of the k-mers to be extracted.
 * @return std::unordered_set<uint32_t> A set containing the hash values of all extracted k-mers.
 *
 * @throws std::runtime_error if the sequence length is shorter than the specified k-mer length.
 */
std::unordered_set<uint32_t> extract_and_hash_kmers_nthash(const std::string &sequence, int k);

#endif // KMER_H
