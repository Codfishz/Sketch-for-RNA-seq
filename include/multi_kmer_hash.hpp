#ifndef MULTI_KMER_HASH_HPP
#define MULTI_KMER_HASH_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

class MultiKmerHash {
    std::string seq;
    std::vector<unsigned> ks;
    std::unordered_map<unsigned, uint64_t> hashes;
    size_t pos;

    void init();

public:
    MultiKmerHash(const std::string& sequence, const std::vector<unsigned>& k_lengths);
    
    bool roll();

    const std::unordered_map<unsigned, uint64_t>& current_hashes() const;

    size_t current_pos() const;
};

// 在这里添加函数声明
std::unordered_map<unsigned, std::vector<uint64_t>>
extract_hashes_with_multikmer(const std::string& sequence,
                              const std::vector<unsigned>& ks);

#endif // MULTI_KMER_HASH_HPP
