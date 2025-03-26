#include "data_io.h"
#include "kmer.h"
#include "sketch.h"
#include "isoform_assignment.h"
#include "sparse_chaining.h"
#include "multi_kmer_hash.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <chrono>
#include <unordered_map>
#include <nthash/nthash.hpp>
#include <fstream>
#include <sstream>


// 高效保存索引：以二进制格式写入 kmer_lengths、transcripts 和 kmer_to_transcripts
void save_index(const std::string& index_output_path,
                std::vector<unsigned>& kmer_lengths,
                const std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
                const std::unordered_map<std::string, Transcript>& transcripts)
{
    std::ofstream ofs(index_output_path, std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Unable to open file for writing: " << index_output_path << std::endl;
        return;
    }
    
    // 保存 kmer_lengths：先写入数量，再写入每个 unsigned 数值
    size_t kmerCount = kmer_lengths.size();
    ofs.write(reinterpret_cast<const char*>(&kmerCount), sizeof(kmerCount));
    for (unsigned k : kmer_lengths) {
        ofs.write(reinterpret_cast<const char*>(&k), sizeof(k));
    }
    
    // 保存 transcripts：写入 transcript 数量，每个 transcript 写入 id、sequence 和 length
    size_t transcriptCount = transcripts.size();
    ofs.write(reinterpret_cast<const char*>(&transcriptCount), sizeof(transcriptCount));
    for (const auto& [id, transcript] : transcripts) {
        // 写入 id：先写入 id 长度，再写入字符数据
        size_t idLen = id.size();
        ofs.write(reinterpret_cast<const char*>(&idLen), sizeof(idLen));
        ofs.write(id.data(), idLen);
        
        // 写入 sequence 同理
        size_t seqLen = transcript.sequence.size();
        ofs.write(reinterpret_cast<const char*>(&seqLen), sizeof(seqLen));
        ofs.write(transcript.sequence.data(), seqLen);
        
        // 写入 transcript.length
        ofs.write(reinterpret_cast<const char*>(&transcript.length), sizeof(transcript.length));
    }
    
    // 保存 kmer_to_transcripts：
    // 写入 kmer_to_transcripts 的数量（即有多少个不同的 k-mer 长度）
    size_t mapCount = kmer_to_transcripts.size();
    ofs.write(reinterpret_cast<const char*>(&mapCount), sizeof(mapCount));
    for (const auto& [kmer_length, mapping] : kmer_to_transcripts) {
        // 写入 kmer_length（unsigned）
        ofs.write(reinterpret_cast<const char*>(&kmer_length), sizeof(kmer_length));
        // 写入当前 k-mer 长度对应的 TranscriptMapping 数量
        size_t mappingSize = mapping.size();
        ofs.write(reinterpret_cast<const char*>(&mappingSize), sizeof(mappingSize));
        for (const auto& [kmer, vec] : mapping) {
            // 写入 kmer (uint32_t)
            ofs.write(reinterpret_cast<const char*>(&kmer), sizeof(kmer));
            // 写入 vector 的大小
            size_t vecSize = vec.size();
            ofs.write(reinterpret_cast<const char*>(&vecSize), sizeof(vecSize));
            for (const auto& pair : vec) {
                const std::string& transcript_id = pair.first;
                size_t tidLen = transcript_id.size();
                ofs.write(reinterpret_cast<const char*>(&tidLen), sizeof(tidLen));
                ofs.write(transcript_id.data(), tidLen);
                // pair.second 为指针，这里不保存实际内容，读取时设为 nullptr
            }
        }
    }
    
    ofs.close();
    std::cout << "Index saved to " << index_output_path << std::endl;
}


// 高效加载索引：从二进制文件中恢复 kmer_lengths、transcripts 和 kmer_to_transcripts
void load_index(const std::string& index_path,
                std::vector<unsigned>& kmer_lengths,
                std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
                std::unordered_map<std::string, Transcript>& transcripts)
{
    std::ifstream ifs(index_path, std::ios::binary);
    if (!ifs) {
        std::cerr << "Error: Unable to open file for reading: " << index_path << std::endl;
        return;
    }
    
    // 读取 kmer_lengths
    size_t kmerCount = 0;
    ifs.read(reinterpret_cast<char*>(&kmerCount), sizeof(kmerCount));
    kmer_lengths.resize(kmerCount);
    for (size_t i = 0; i < kmerCount; i++) {
        unsigned k;
        ifs.read(reinterpret_cast<char*>(&k), sizeof(k));
        kmer_lengths[i] = k;
    }
    
    // 读取 transcripts
    size_t transcriptCount = 0;
    ifs.read(reinterpret_cast<char*>(&transcriptCount), sizeof(transcriptCount));
    transcripts.clear();
    for (size_t i = 0; i < transcriptCount; i++) {
        // 读取 id
        size_t idLen = 0;
        ifs.read(reinterpret_cast<char*>(&idLen), sizeof(idLen));
        std::string id(idLen, '\0');
        ifs.read(&id[0], idLen);
        
        // 读取 sequence
        size_t seqLen = 0;
        ifs.read(reinterpret_cast<char*>(&seqLen), sizeof(seqLen));
        std::string sequence(seqLen, '\0');
        ifs.read(&sequence[0], seqLen);
        
        int length;
        ifs.read(reinterpret_cast<char*>(&length), sizeof(length));
        
        transcripts[id] = Transcript{id, sequence, length};
    }
    
    // 读取 kmer_to_transcripts
    size_t mapCount = 0;
    ifs.read(reinterpret_cast<char*>(&mapCount), sizeof(mapCount));
    kmer_to_transcripts.clear();
    for (size_t i = 0; i < mapCount; i++) {
        unsigned kmer_length;
        ifs.read(reinterpret_cast<char*>(&kmer_length), sizeof(kmer_length));
        size_t mappingSize = 0;
        ifs.read(reinterpret_cast<char*>(&mappingSize), sizeof(mappingSize));
        TranscriptMapping mapping;
        for (size_t j = 0; j < mappingSize; j++) {
            uint32_t kmer;
            ifs.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
            size_t vecSize = 0;
            ifs.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
            std::vector<std::pair<std::string, const SketchType*>> vec;
            for (size_t h = 0; h < vecSize; h++) {
                size_t tidLen = 0;
                ifs.read(reinterpret_cast<char*>(&tidLen), sizeof(tidLen));
                std::string transcript_id(tidLen, '\0');
                ifs.read(&transcript_id[0], tidLen);
                vec.emplace_back(transcript_id, nullptr);
            }
            mapping[kmer] = vec;
        }
        kmer_to_transcripts[kmer_length] = mapping;
    }
    
    ifs.close();
    std::cout << "Index loaded from " << index_path << std::endl;
}


// 帮助信息打印函数
void print_help(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS] <mode> [arguments]" << std::endl;
    std::cout << "Modes:" << std::endl;
    std::cout << "  index   Build index from reference genome" << std::endl;
    std::cout << "  quant   Quantify using pre-built index and reads" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help              Show this help message and exit" << std::endl;
    std::cout << "  -k, --kmer-length SIZE  Comma separated list of k-mer lengths (default: 81)" << std::endl;
    std::cout << "  -w, --window-size SIZE  Window size for sparse chaining (default: 100)" << std::endl;
    std::cout << "  -o, --mode MODE         Mode: index or quant (default: quant)" << std::endl;
    std::cout << std::endl;
    std::cout << "Index mode usage:" << std::endl;
    std::cout << "  " << program_name << " index <reference_genome.fasta> <index_output>" << std::endl;
    std::cout << std::endl;
    std::cout << "Quant mode usage:" << std::endl;
    std::cout << "  " << program_name << " quant <index_file> <reads.fastq> <output>" << std::endl;
}

// 全局常量 sketch_size（可以根据需要调整）
const float sketch_size = 0.05;

// index 模式：构建索引并保存
void build_and_save_index(const std::string& reference_genome_path,
                          const std::string& index_output_path,
                          std::vector<unsigned>& kmer_lengths)
{
    auto start = std::chrono::high_resolution_clock::now();
    // 加载参考转录本
    std::unordered_map<std::string, Transcript> transcripts = load_fasta(reference_genome_path);
    std::unordered_map<std::string, MultiKmerSketch> transcript_sketches;
    for (const auto& [id, transcript] : transcripts) {
        bool valid = true;
        for (unsigned k : kmer_lengths) {
            if (transcript.sequence.size() < k) {
                valid = false;
                break;
            }
        }
        if (!valid)
            continue;

        MultiKmerSketch mks;
        for (unsigned k : kmer_lengths) {
            auto hashed_kmers = extract_and_hash_kmers_nthash(transcript.sequence, k);
            mks.sketches[k] = createSketch_FracMinhash(hashed_kmers, sketch_size);
        }
        transcript_sketches[id] = std::move(mks);
    }
    // 构建 k-mer 到转录本的映射
    auto kmer_to_transcripts = build_kmer_to_transcript_map(transcript_sketches);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Index built in " << elapsed.count() << " seconds." << std::endl;

    // 保存索引到文件
    save_index(index_output_path, kmer_lengths, kmer_to_transcripts, transcripts);
}

// quant 模式：加载 index 进行 reads 量化
void quantification(const std::string& index_path,
                    const std::string& reads_path,
                    const std::string& output_path,
                    std::vector<unsigned>& kmer_lengths)
{
    auto start = std::chrono::high_resolution_clock::now();
    // 加载预构建的 index
    std::vector<unsigned> loaded_kmer_lengths;
    std::unordered_map<unsigned, TranscriptMapping> kmer_to_transcripts;
    std::unordered_map<std::string, Transcript> transcripts;
    load_index(index_path,kmer_lengths,kmer_to_transcripts,transcripts);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Loading index completed in " << elapsed.count() << " seconds." << std::endl;
    // 若加载的 kmer 长度为空，则使用命令行指定的
    std::vector<unsigned> effective_kmer_lengths = loaded_kmer_lengths.empty() ? kmer_lengths : loaded_kmer_lengths;

    // 加载 reads
    auto reads = load_fastq(reads_path);
    std::unordered_map<std::string, MultiKmerSketch> read_sketches;
    for (const auto& [id, read] : reads) {
        bool valid = true;
        for (unsigned k : effective_kmer_lengths) {
            if (read.sequence.size() < k) {
                valid = false;
                break;
            }
        }
        if (!valid)
            continue;

        MultiKmerSketch mks;
        for (unsigned k : effective_kmer_lengths) {
            auto hashed_kmers = extract_and_hash_kmers_nthash(read.sequence, k);
            mks.sketches[k] = createSketch_FracMinhash(hashed_kmers, sketch_size);
        }
        read_sketches[id] = std::move(mks);
    }

    // 执行 sparse chaining 和后续量化步骤
    start = std::chrono::high_resolution_clock::now();
    auto homologous_segments = sparse_chain(read_sketches, kmer_to_transcripts, transcripts, reads, effective_kmer_lengths, 0.75);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Sparse chaining completed in " << elapsed.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    auto pi = estimate_isoform_abundance_em(homologous_segments, transcripts, 20, 0.01);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "EM estimation completed in " << elapsed.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    auto read_counts = assign_reads_to_isoforms(homologous_segments, pi, transcripts);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Read assignment completed in " << elapsed.count() << " seconds." << std::endl;

    start = std::chrono::high_resolution_clock::now();
    auto tpm = calculate_tpm(read_counts, transcripts);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "TPM calculation completed in " << elapsed.count() << " seconds." << std::endl;

    // 输出结果
    output_to_csv(output_path, read_counts, tpm, pi, transcripts);
    std::cout << "Output written to " << output_path << std::endl;
}

int main(int argc, char* argv[]) {
    // 默认模式为 quant，默认 k-mer 长度为 31
    std::string mode = "quant";
    std::vector<unsigned> kmer_lengths = { 67 };

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"kmer-length", required_argument, 0, 'k'},
        {"window-size", required_argument, 0, 'w'},
        {"mode", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "hk:w:o:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_help(argv[0]);
                return 0;
            case 'k': {
                std::string s(optarg);
                kmer_lengths.clear();  // 若用户提供则清除默认值
                std::istringstream iss(s);
                std::string token;
                while (std::getline(iss, token, ',')) {
                    if (!token.empty()) {
                        kmer_lengths.push_back(std::stoi(token));
                    }
                }
                break;
            }
            case 'o':
                mode = std::string(optarg);
                break;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    // 根据不同模式解析位置参数
    if (mode == "index") {
        if (optind + 2 > argc) {
            std::cerr << "Usage: " << argv[0] << " index <reference_genome.fasta> <index_output>" << std::endl;
            return 1;
        }
        std::string reference_genome_path = argv[optind];
        std::string index_output_path = argv[optind + 1];
        build_and_save_index(reference_genome_path, index_output_path, kmer_lengths);
    } else if (mode == "quant") {
        if (optind + 3 > argc) {
            std::cerr << "Usage: " << argv[0] << " quant <index_file> <reads.fastq> <output>" << std::endl;
            return 1;
        }
        std::string index_path = argv[optind];
        std::string reads_path = argv[optind + 1];
        std::string output_path = argv[optind + 2];
        quantification(index_path, reads_path, output_path, kmer_lengths);
    } else {
        std::cerr << "Invalid mode. Please choose 'index' or 'quant'." << std::endl;
        return 1;
    }

    return 0;
}
