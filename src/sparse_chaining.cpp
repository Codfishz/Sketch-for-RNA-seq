#include "sparse_chaining.h"
#include "kmer.h"
#include "sketch.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <string>
#include <utility>

// Function to perform sparse chaining
std::unordered_map<std::string, std::string> sparse_chain(
    const std::unordered_map<std::string, std::unordered_set<uint32_t>>& read_sketches,
    const std::unordered_map<uint32_t, std::vector<std::pair<std::string, const std::unordered_set<uint32_t>*>>>& kmer_to_transcripts,
    const std::unordered_map<std::string, Transcript>& transcripts,
    const std::unordered_map<std::string, Read>& reads,
    int kmer_length, double fraction) {

    std::unordered_map<std::string, std::string> homologous_segments;
    for (const auto& [read_id, read_sketch] : read_sketches) {
        std::unordered_map<std::string, int> transcript_match_counts;


        for (const auto& kmer : read_sketch) {
            if (kmer_to_transcripts.find(kmer) != kmer_to_transcripts.end()) {
                for (const auto& [transcript_id, transcript_sketch] : kmer_to_transcripts.at(kmer)) {
                    transcript_match_counts[transcript_id]++;
                }
            }
        }

        std::vector<std::string> candidate_transcripts;
        int max_match_count = 0;
        for (const auto& [transcript_id, match_count] : transcript_match_counts) {
            if (match_count > max_match_count) {
                max_match_count = match_count;
                candidate_transcripts.clear();
                candidate_transcripts.push_back(transcript_id);
            } else if (match_count == max_match_count) {
                candidate_transcripts.push_back(transcript_id);
            }
        }
        // std::string target_transcript_id = "ENST00000692674.1|ENSG00000289007.2|-|-|ENST00000692674|ENSG00000289007|1520|lncRNA|";
        // if (std::find(candidate_transcripts.begin(), candidate_transcripts.end(), target_transcript_id) != candidate_transcripts.end()) {
        //     std::cout << target_transcript_id << "In candidate" << std::endl;
        // }

        if (!candidate_transcripts.empty()) {
            std::string best_transcript_id = find_best_match_orderedminhash(read_id, reads, candidate_transcripts, transcripts, kmer_length, fraction);
            if (!best_transcript_id.empty()) {
                homologous_segments[read_id] = best_transcript_id;
            }
        }
    }

    return homologous_segments;
}

// Function to find the best matching transcript using Ordered MinHash for candidates
std::string find_best_match_orderedminhash(const std::string& read_id,
                                           const std::unordered_map<std::string, Read>& reads,
                                           const std::vector<std::string>& candidate_transcripts,
                                           const std::unordered_map<std::string, Transcript>& transcripts,
                                           int kmer_length, double fraction) {
    // Get the read sequence
    auto read_it = reads.find(read_id);
    if (read_it == reads.end()) {
        std::cerr << "Error: Read ID not found: " << read_id << std::endl;
        return "";
    }
    const auto& read_sequence = read_it->second.sequence;
    if (read_sequence.length() < kmer_length) {
        std::cerr << "Error: Read length is shorter than kmer length" << std::endl;
        return "";
    }

    // Create Ordered MinHash sketch for the read
    auto read_hashed_kmers = extract_and_hash_kmers_with_positions_murmur(read_sequence, kmer_length,12345);
    if (read_hashed_kmers.empty()) {
        std::cerr << "Error: No hashed kmers generated for read ID: " << read_id << std::endl;
        return "";
    }
    auto read_ordered_sketch = createSketch_OrderedMinhash(read_hashed_kmers, fraction);

    std::string best_transcript_id;
    double best_similarity_score = 0.0;


    std::string target_transcript_id = "ENST00000692674.1|ENSG00000289007.2|-|-|ENST00000692674|ENSG00000289007|1520|lncRNA|";
    double target_score;
    // int count = 0;
    // Iterate through candidate transcripts

    for (const auto& candidate_transcript_id : candidate_transcripts) {
        auto transcript_it = transcripts.find(candidate_transcript_id);
        if (transcript_it == transcripts.end()) {
            std::cerr << "Warning: Transcript ID not found: " << candidate_transcript_id << std::endl;
            continue;
        }
        const auto& transcript_sequence = transcript_it->second.sequence;
        if (transcript_sequence.length() < kmer_length) {
            std::cerr << "Warning: Transcript length is shorter than kmer length for transcript ID: " << candidate_transcript_id << std::endl;
            continue;
        }

        auto transcript_hashed_kmers = extract_and_hash_kmers_with_positions_murmur(transcript_sequence, kmer_length,12345);
        if (transcript_hashed_kmers.empty()) {
            std::cerr << "Warning: No hashed kmers generated for transcript ID: " << candidate_transcript_id << std::endl;
            continue;
        }
        auto transcript_ordered_sketch = createSketch_OrderedMinhash(transcript_hashed_kmers, fraction);


        double similarity_score = compare_relative_positions(read_ordered_sketch, transcript_ordered_sketch);
        // if (similarity_score==1){
        //     count++;
        // }
        if (candidate_transcript_id == target_transcript_id){
            target_score = similarity_score;
        }
        if (similarity_score > best_similarity_score) {
            best_similarity_score = similarity_score;
            best_transcript_id = candidate_transcript_id;
        } 
    }
    // if (best_transcript_id!=target_transcript_id){
    //     std::cout<<"best ID"<<best_transcript_id<<" with score "<<best_similarity_score<<std::endl;
    //     std::cout<<"target ID"<<target_transcript_id<<" with score "<<target_score<<std::endl;
    //     // std::cout<<"Num of candidate:" <<candidate_transcripts.size();
    //     // std::cout<<" count "<<count<<std::endl;
    // } 


    return best_transcript_id;
}

// Function to compare positions of common hashes between read and transcript
double compare_relative_positions(const std::vector<std::pair<uint32_t, size_t>>& read_minhash,
                                  const std::vector<std::pair<uint32_t, size_t>>& transcript_minhash) {
    // Create a map from hash to position for both read and transcript
    std::unordered_map<uint32_t, size_t> read_positions;
    for (const auto& [hash, position] : read_minhash) {
        read_positions[hash] = position;
    }

    std::unordered_map<uint32_t, size_t> transcript_positions;
    for (const auto& [hash, position] : transcript_minhash) {
        transcript_positions[hash] = position;
    }

    // Extract common hashes
    std::vector<uint32_t> common_hashes;
    for (const auto& [hash, position] : read_positions) {
        if (transcript_positions.find(hash) != transcript_positions.end()) {
            common_hashes.push_back(hash);
        }
    }

    // Compare relative positions
    size_t similar_order_count = 0;
    size_t total_common = common_hashes.size();
    for (size_t i = 0; i < total_common - 1; ++i) {
        uint32_t hash1 = common_hashes[i];
        uint32_t hash2 = common_hashes[i + 1];
        if ((read_positions[hash1] < read_positions[hash2] && transcript_positions[hash1] < transcript_positions[hash2]) ||
            (read_positions[hash1] > read_positions[hash2] && transcript_positions[hash1] > transcript_positions[hash2])) {
            similar_order_count++;
        }
    }

    // Return a score that represents the proportion of common hashes with similar relative order
    return total_common > 1 ? static_cast<double>(similar_order_count) / (total_common - 1) : 0.0;
}