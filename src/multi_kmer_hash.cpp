#include "multi_kmer_hash.hpp"
#include <algorithm>
#include <chrono>
static const uint8_t CONVERT_TAB[256] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3
};

inline uint64_t srol(uint64_t val, int shift = 1) {
    return (val << shift) | (val >> (64 - shift));
}

static const uint64_t SEED_TAB[4] = {
    0x3c8bfbb395c60474ULL, 0x3193c18562a02b4cULL,
    0x20323ed082572324ULL, 0x295549f54be24456ULL
};

uint64_t base_hash(const char* seq, unsigned k) {
    uint64_t h = 0;
    for (unsigned i = 0; i < k; i++) {
        h = srol(h) ^ SEED_TAB[CONVERT_TAB[(unsigned char)seq[i]]];
    }
    return h;
}

uint64_t next_hash(uint64_t prev_hash, unsigned k, char char_out, char char_in) {
    return srol(prev_hash) 
           ^ SEED_TAB[CONVERT_TAB[(unsigned char)char_in]]
           ^ srol(SEED_TAB[CONVERT_TAB[(unsigned char)char_out]], k);
}

MultiKmerHash::MultiKmerHash(const std::string& sequence, const std::vector<unsigned>& k_lengths)
    : seq(sequence), ks(k_lengths), pos(0) 
{
    init();
}

void MultiKmerHash::init() {
    hashes.clear();
    for (auto k : ks) {
        if (seq.size() >= k)
            hashes[k] = base_hash(seq.data(), k);
    }
}

bool MultiKmerHash::roll() {
    if (++pos + *std::max_element(ks.begin(), ks.end()) > seq.size())
        return false;

    for (auto k : ks) {
        if (pos + k <= seq.size()) {
            hashes[k] = next_hash(
                hashes[k], k,
                seq[pos - 1],
                seq[pos + k - 1]
            );
        }
    }
    return true;
}

const std::unordered_map<unsigned, uint64_t>& MultiKmerHash::current_hashes() const {
    return hashes;
}

size_t MultiKmerHash::current_pos() const {
    return pos;
}


std::unordered_map<unsigned, std::vector<uint64_t>>
extract_hashes_with_multikmer(const std::string& sequence,
                              const std::vector<unsigned>& ks) {
    MultiKmerHash mkh(sequence, ks);
    std::unordered_map<unsigned, std::vector<uint64_t>> hash_results;
    for (unsigned k : ks)
        hash_results[k] = {};

    do {
        auto current_hashes = mkh.current_hashes();
        for (auto& [k, hash] : current_hashes) {
            hash_results[k].push_back(hash);
        }
    } while (mkh.roll());

    return hash_results;
}