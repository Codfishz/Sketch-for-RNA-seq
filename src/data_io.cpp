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

// 读取FASTA文件
std::unordered_map<std::string, Transcript> load_fasta(const std::string& fasta_file) {
    std::unordered_map<std::string, Transcript> transcripts;
    std::ifstream infile(fasta_file, std::ios::in | std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Could not open FASTA file: " + fasta_file);
    }

    std::string line, sequence;
    std::string current_id;

    // 设置输入缓冲区大小（可选）
    infile.rdbuf()->pubsetbuf(nullptr, 1 << 20); // 设置 1MB 的缓冲区

    while (std::getline(infile, line)) {
        if (line.empty()) continue; // 跳过空行

        if (line[0] == '>') {
            if (!current_id.empty() && is_valid_sequence(sequence)) {
                transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
            }
            current_id = line.substr(1, line.find(' ') - 1);  // 提取 ID，去掉 '>'
            sequence.clear();
            sequence.reserve(1000);  // 预留空间，减少 realloc 次数
        } else {
            sequence += line;
        }
    }
    // 处理最后一个转录本
    if (!current_id.empty()) {
        transcripts.emplace(std::move(current_id), Transcript{current_id, std::move(sequence), static_cast<int>(sequence.size())});
    }

    return transcripts;
}

// 读取FASTQ文件并返回一个hash table
std::unordered_map<std::string, Read> load_fastq(const std::string& fastq_file) {
    int count = 0;
    std::unordered_map<std::string, Read> reads;
    std::ifstream infile(fastq_file);
    if (!infile) {
        throw std::runtime_error("Could not open FASTQ file: " + fastq_file);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] != '@') continue; // 跳过无效行
        Read read;
        read.id = line.substr(1);  // 提取 read ID
        std::getline(infile, read.sequence);
        std::getline(infile, line); // 跳过 '+'
        std::getline(infile, read.quality);
        count++;
        if (is_valid_sequence(read.sequence)) {
            reads[read.id] = read;  // 将合法的 read 存入 hash table
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

    // 写入文件头
    outfile << "Name,Length,EffectiveLength,TPM,NumReads\n";

    // 写入每个转录本的统计信息
    for (const auto& [id, transcript] : transcripts) {
        auto read_count_it = read_counts.find(id);
        auto tpm_it = tpms.find(id);

        if (read_count_it != read_counts.end() && tpm_it != tpms.end()) {
            int read_count = read_count_it->second;
            double tpm = tpm_it->second;
            outfile << id << "," << transcript.length << "," 
                    << transcript.length << "," // 假设 EffectiveLength 和 Length 相同
                    << tpm << "," << read_count << "\n";
        } 
    }

    outfile.close();
}

// // 质量控制函数：过滤掉低质量的reads
// std::vector<Read> quality_control(const std::vector<Read>& reads) {
//     std::vector<Read> filtered_reads;
//     const int min_quality = 20; // 设置最小质量分数
//     const int min_length = 50;  // 设置最小read长度

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
