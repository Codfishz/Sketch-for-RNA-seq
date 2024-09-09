#include "data_io.h"
#include "kmer.h"
#include "sketch.h"
#include "isoform_assignment.h"
#include "sparse_chaining.h"
//#include "MurmurHash3.h"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <chrono>
#include <unordered_map> 
//
//帮助函数
void print_help(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS] <reference_genome.fasta> <reads.fastq>" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help              Show this help message and exit" << std::endl;
    std::cout << "  -t, --threads NUM       Number of threads to use (default: 1)" << std::endl;
    std::cout << "  -m, --memory SIZE       Maximum memory to use (default: 1024 MB)" << std::endl;
    std::cout << "  -k, --kmer-length SIZE  Length of k-mers (default: 31)" << std::endl;
    std::cout << "  -w, --window-size SIZE  Window size for sparse chaining (default: 100)" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  <reference_genome.fasta>  Path to the reference genome in FASTA format" << std::endl;
    std::cout << "  <reads.fastq>             Path to the raw RNA-seq data in FASTQ format" << std::endl;
    std::cout << std::endl;
    std::cout << "Example:" << std::endl;
    std::cout << "  " << program_name << " -t 4 -m 2048 -k 21 reference.fasta reads.fastq" << std::endl;
}

int main(int argc, char* argv[]) {
    // 默认参数值
    int threads = 1;
    int memory = 1024; // 以MB为单位
    int kmer_length = 31;
    int window_size = 100;

    // 处理命令行选项
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"threads", required_argument, 0, 't'},
        {"memory", required_argument, 0, 'm'},
        {"kmer-length", required_argument, 0, 'k'},
        {"window-size", required_argument, 0, 'w'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "ht:m:k:w:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                print_help(argv[0]);
                return 0;
            case 't':
                threads = std::stoi(optarg);
                break;
            case 'm':
                memory = std::stoi(optarg);
                break;
            case 'k':
                kmer_length = std::stoi(optarg);
                break;
            case 'w':
                window_size = std::stoi(optarg);
                break;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    if (optind + 3 > argc) {
        std::cerr << "Usage: " << argv[0] << " <reference_genome.fasta> <reads.fastq> <gene_annotation.gtf>" << std::endl;
        return 1;
    }

    std::string reference_genome_path = argv[optind];
    std::string reads_path = argv[optind + 1];
    std::string annotation_path = argv[optind + 2];
    std::string output_path = argv[optind + 3];

    try {
        auto start = std::chrono::high_resolution_clock::now();

        // 解析 GTF 文件，加载转录本和注释信息
        // auto annotations_map = parse_gtf(annotation_path);
        // std::unordered_map<std::string, Transcript> transcripts = load_transcripts_with_annotations(reference_genome_path, annotations_map);
        std::unordered_map<std::string, Transcript> transcripts = load_fasta(reference_genome_path);

        int count = 0;

        // 创建转录本的 sketches
        std::unordered_map<std::string, std::unordered_set<uint32_t>> transcript_sketches;
        for (const auto& [id, transcript] : transcripts) {
            if (transcript.sequence.length() < kmer_length) {
                // std::cerr << "Skipping transcript " << id << " because its length (" << transcript.sequence.length() 
                //         << ") is less than kmer-length (" << kmer_length << ")." << std::endl;
                continue;  // 跳过该转录本
            }

            auto hashed_kmers = extract_and_hash_kmers(transcript.sequence, kmer_length);
            auto sketch = createSketch_FracMinhash(hashed_kmers, 0.2);
            transcript_sketches[id] = sketch;
        }

        auto kmer_to_transcripts = build_kmer_to_transcript_map(transcript_sketches);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        // 输出转录本的数量和运行时间
        std::cout << "Number of transcripts loaded: " << transcripts.size() << std::endl;
        std::cout << "Time taken to load and process transcripts: " << elapsed.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // 加载 reads
        std::unordered_map<std::string, Read> reads = load_fastq(reads_path);
        
        count = 0;
        // 创建 reads 的 sketches
        std::unordered_map<std::string, std::unordered_set<uint32_t>> read_sketches;
        for (const auto& [id, read] : reads) {
            if (read.sequence.length() < kmer_length) {
                // std::cerr << "Skipping read " << id << " because its length (" << read.sequence.length() 
                //         << ") is less than kmer-length (" << kmer_length << ")." << std::endl;
                continue;  // 跳过该转录本
            }
            auto hashed_kmers = extract_and_hash_kmers(read.sequence, kmer_length);
            auto sketch = createSketch_FracMinhash(hashed_kmers, 0.2);
            read_sketches[id] = sketch;
        }

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;

        // 输出 reads 的数量和运行时间
        std::cout << "Number of reads loaded: " << reads.size() << std::endl;
        std::cout << "Time taken to load and process reads: " << elapsed.count() << " seconds." << std::endl;
   
        //构建sparse chain
        // start = std::chrono::high_resolution_clock::now();
        // std::vector<std::vector<std::pair<int, int>>> all_homologous_segments;
        // for (size_t i = 0; i < read_sketches.size(); ++i) {        
        //     for (size_t j = 0; j < transcript_sketches.size(); ++j) {
        //         if(j%5000==0){
        //             std::cout<<j<<std::endl;
        //         }   
        //         auto homologous_segments = sparse_chain(read_sketches[i], transcript_sketches[j]);
        //         all_homologous_segments.push_back(homologous_segments);
        //         // std::cout << "Read " << i << " and Transcript " << j << " homologous segments:" << std::endl;
        //         // for (const auto& segment : homologous_segments) {
        //         //     std::cout << "Reference position: " << segment.first << ", Read position: " << segment.second << std::endl;
        //         // }
        //     }
        // }
        // end = std::chrono::high_resolution_clock::now();
        // elapsed = end - start;
        // std::cout << "Time taken for sparse chaining: " << elapsed.count() << " seconds." << std::endl;
        // std::cout << "All homologous segments found:" << std::endl;
        // for (size_t i = 0; i < all_homologous_segments.size(); ++i) {
        //     std::cout << "For Read " << i << ":" << std::endl;
        //     for (const auto& segment : all_homologous_segments[i]) {
        //         std::cout << "Reference position: " << segment.first << ", Read position: " << segment.second << std::endl;
        //     }
        // }
        // end = std::chrono::high_resolution_clock::now();
        // elapsed = end - start;
        // std::cout << "Time taken for sparse chaining: " << elapsed.count() << " seconds." << std::endl;
        start = std::chrono::high_resolution_clock::now();

        auto homologous_segments = sparse_chain(read_sketches, kmer_to_transcripts);
        std::cout << "Total number of successful matches: " << homologous_segments.size() << std::endl;
        std::cout<<"finish sparse chain"<<std::endl;
        // 输出匹配结果
        // count =0;
        // for (const auto& [read_id, match] : homologous_segments) {
        //     const std::string& transcript_id = match.first;
        //     int match_count = match.second;

        //     std::cout << "Read ID: " << read_id << " best matches with Transcript ID: " 
        //             << transcript_id << " with " << match_count << " matching k-mers." << std::endl;
        // }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken to create sparse chain: " << elapsed.count() << " seconds." << std::endl;


        // 分配reads给 isoforms 和计算TPM
        start = std::chrono::high_resolution_clock::now();

        // 假设 `transcripts` 已经包含了加载好的 `Transcript` 结构体
        std::unordered_map<std::string, int> read_counts = assign_reads_to_isoforms(homologous_segments, transcripts);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken for read-isoform assignment: " << elapsed.count() << " seconds." << std::endl;

        start = std::chrono::high_resolution_clock::now();

        std::unordered_map<std::string, double> tpm = calculate_tpm(read_counts, transcripts);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "Time taken to calculate TPM: " << elapsed.count() << " seconds." << std::endl;

        //输出TPM结果
        // for (const auto& [transcript_id, tpm_value] : tpm) {
        //     std::cout << "Transcript " << transcript_id << " has TPM " << tpm_value << std::endl;
        // }

        output_to_csv(output_path, read_counts, tpm, transcripts);

    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
