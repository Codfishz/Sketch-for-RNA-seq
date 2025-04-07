#include "kmer.h"
#include <stdexcept>
#include <unordered_set>
#include <nthash/nthash.hpp>
#include <string>

/**
 * @brief Extract and hash k-mers from a given DNA sequence using nthash.
 *
 * This function extracts all k-mers from the input sequence and computes their hash values
 * using the nthash library. It returns an unordered_set of uint32_t hash values.
 *
 * @param sequence The input DNA sequence from which k-mers are extracted.
 * @param k The length of each k-mer.
 * @return std::unordered_set<uint32_t> A set containing hash values of the extracted k-mers.
 *
 * @throw std::runtime_error if the sequence length is shorter than the specified k-mer length.
 */
std::unordered_set<uint32_t> extract_and_hash_kmers_nthash(const std::string &sequence, int k) {
    if (sequence.size() < static_cast<size_t>(k)) {
        throw std::runtime_error("Sequence length is shorter than k-mer length");
    }

    std::unordered_set<uint32_t> hash_set;
    // Initialize nthash with the sequence, seed value 1, and k-mer length k
    nthash::NtHash nth(sequence, 1, k);

    // Iterate through the sequence to compute forward hashes for each k-mer
    while (nth.roll()) {
        uint32_t hash_value = nth.get_forward_hash();
        hash_set.insert(hash_value);
    }
    
    return hash_set;
}
