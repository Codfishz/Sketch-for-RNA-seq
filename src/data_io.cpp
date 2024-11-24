#include "data_io.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <unordered_set>

// 读取FASTQ文件
#include <unordered_map>

bool is_valid_sequence(const std::string& sequence) {
    static bool valid_char[256] = { false };

    // 初始化合法字符，只在第一次调用时执行
    if (!valid_char['A']) {
        valid_char['A'] = true;
        valid_char['T'] = true;
        valid_char['C'] = true;
        valid_char['G'] = true;
    }

    for (char c : sequence) {
        if (!valid_char[static_cast<unsigned char>(c)]) {
            return false;
        }
    }
    return true;
}

// Read FASTA
std::unordered_map<std::string, Transcript> load_fasta(const std::string& fasta_file) {
    std::unordered_map<std::string, Transcript> transcripts;
    std::ifstream infile(fasta_file, std::ios::in | std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Could not open FASTA file: " + fasta_file);
    }

    std::string line, sequence;
    std::string current_id;

    infile.rdbuf()->pubsetbuf(nullptr, 1 << 20); // 

    while (std::getline(infile, line)) {
        if (line.empty()) continue; 

        if (line[0] == '>') {
            if (!current_id.empty() && is_valid_sequence(sequence)) {
                transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
            }
            current_id = line.substr(1, line.find(' ') - 1);  
            sequence.clear();
            sequence.reserve(1000);  
        } else {
            sequence += line;
        }
    }
    if (!current_id.empty()) {
        transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
    }

    return transcripts;
}

std::unordered_map<std::string, Read> load_fastq(const std::string& fastq_file) {
    int count = 0;
    std::unordered_map<std::string, Read> reads;
    std::ifstream infile(fastq_file);
    if (!infile) {
        throw std::runtime_error("Could not open FASTQ file: " + fastq_file);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] != '@') continue; 
        Read read;
        read.id = line.substr(1);  
        std::getline(infile, read.sequence);
        std::getline(infile, line); 
        std::getline(infile, read.quality);
        count++;
        if (is_valid_sequence(read.sequence)) {
            reads[read.id] = read;  
        } 
    }
    std::cout<<"Actual number of reads: "<<count<<std::endl;
    return reads;
}

void output_to_csv(const std::string& filename, 
                   const std::unordered_map<std::string, int>& read_counts,
                   const std::unordered_map<std::string, double>& tpms,
                   const std::unordered_map<std::string, Transcript>& transcripts) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    outfile << "Name,Length,EffectiveLength,TPM,NumReads\n";

    for (const auto& [id, transcript] : transcripts) {
        auto read_count_it = read_counts.find(id);
        auto tpm_it = tpms.find(id);

        if (read_count_it != read_counts.end() && tpm_it != tpms.end()) {
            int read_count = read_count_it->second;
            double tpm = tpm_it->second;
            outfile << id << "," << transcript.length << "," 
                    << transcript.length << "," // 
                    << tpm << "," << read_count << "\n";
        } 
    }

    outfile.close();
}

// std::vector<Read> quality_control(const std::vector<Read>& reads) {
//     std::vector<Read> filtered_reads;
//     const int min_quality = 20; 
//     const int min_length = 50;  

//     for (const auto& read : reads) {
//         bool pass = true;
//         for (char q : read.quality) {
//             if ((q - 33) < min_quality) {
//                 pass = false;
//                 break;
//             }
//         }
//         if (pass && read.sequence.size() >= min_length) {
//             filtered_reads.push_back(read);
//         }
//     }
//     return filtered_reads;
// }
