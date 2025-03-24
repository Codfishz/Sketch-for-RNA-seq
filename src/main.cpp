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

void save_index(const std::string& index_output_path,
	std::vector<unsigned>& kmer_lengths,
	const std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
	const std::unordered_map<std::string, Transcript>& transcripts)
{
	std::ofstream ofs(index_output_path);
	if (!ofs) {
		std::cerr << "Error: Unable to open file for writing: " << index_output_path << std::endl;
		return;
	}

	// 保存 kmer_lengths
	ofs << kmer_lengths.size() << "\n";
	for (auto k : kmer_lengths) {
		ofs << k << "\n";
	}

	// 保存 transcripts（每行保存：id <tab> sequence <tab> length）
	ofs << transcripts.size() << "\n";
	for (const auto& [id, transcript] : transcripts) {
		ofs << id << "\t" << transcript.sequence << "\t" << transcript.length << "\n";
	}

	// 保存 kmer_to_transcripts
	// 第一行写入 kmer_to_transcripts 中有多少个不同的 k-mer 长度（即 map 的 size）
	ofs << kmer_to_transcripts.size() << "\n";
	for (const auto& [kmer_length, mapping] : kmer_to_transcripts) {
		ofs << kmer_length << "\n";
		// 写入当前 k-mer 长度对应的 TranscriptMapping 中有多少个 k-mer（map 的 size）
		ofs << mapping.size() << "\n";
		for (const auto& [kmer, vec] : mapping) {
			ofs << kmer << "\n";
			// 写入 vector 中有多少个 transcript 记录
			ofs << vec.size() << "\n";
			for (const auto& pair : vec) {
				// 保存 transcript id；忽略 pair.second（指针部分）
				ofs << pair.first << "\n";
			}
		}
	}

	ofs.close();
	std::cout << "Index saved to " << index_output_path << std::endl;
}



void load_index(const std::string& index_path,
	std::vector<unsigned>& kmer_lengths,
	std::unordered_map<unsigned, TranscriptMapping>& kmer_to_transcripts,
	std::unordered_map<std::string, Transcript>& transcripts)
{
	std::ifstream ifs(index_path);
	if (!ifs) {
		std::cerr << "Error: Unable to open file for reading: " << index_path << std::endl;
		return;
	}

	std::string line;
	size_t count = 0;

	// 读取 kmer_lengths
	if (std::getline(ifs, line)) {
		std::istringstream iss(line);
		iss >> count;
		kmer_lengths.clear();
		for (size_t i = 0; i < count; i++) {
			if (std::getline(ifs, line)) {
				unsigned k;
				std::istringstream kss(line);
				kss >> k;
				kmer_lengths.push_back(k);
			}
		}
	}

	// 读取 transcripts
	if (std::getline(ifs, line)) {
		std::istringstream iss(line);
		iss >> count;
		transcripts.clear();
		for (size_t i = 0; i < count; i++) {
			if (std::getline(ifs, line)) {
				std::istringstream lineStream(line);
				std::string id, sequence;
				int length;
				// 使用制表符分隔 id 和 sequence，再读入 length
				if (std::getline(lineStream, id, '\t') &&
					std::getline(lineStream, sequence, '\t') &&
					(lineStream >> length)) {
					transcripts[id] = Transcript{id, sequence, length};
				}
			}
		}
	}

	// 读取 kmer_to_transcripts
	if (std::getline(ifs, line)) {
		std::istringstream iss(line);
		iss >> count;
		kmer_to_transcripts.clear();
		for (size_t i = 0; i < count; i++) {
			unsigned kmer_length = 0;
			if (std::getline(ifs, line)) {
				std::istringstream iss_k(line);
				iss_k >> kmer_length;
			}
			size_t mappingSize = 0;
			if (std::getline(ifs, line)) {
				std::istringstream iss_map(line);
				iss_map >> mappingSize;
			}
			TranscriptMapping mapping;
			for (size_t j = 0; j < mappingSize; j++) {
				uint32_t kmer = 0;
				if (std::getline(ifs, line)) {
					std::istringstream iss_kmer(line);
					iss_kmer >> kmer;
				}
				size_t vecSize = 0;
				if (std::getline(ifs, line)) {
					std::istringstream iss_vec(line);
					iss_vec >> vecSize;
				}
				std::vector<std::pair<std::string, const SketchType*>> vec;
				for (size_t h = 0; h < vecSize; h++) {
					std::string transcript_id;
					if (std::getline(ifs, transcript_id)) {
						// 由于我们没有保存 sketch 指针，故将其设为 nullptr
						vec.emplace_back(transcript_id, nullptr);
					}
				}
				mapping[kmer] = vec;
			}
			kmer_to_transcripts[kmer_length] = mapping;
		}
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
	// 加载预构建的 index
	std::vector<unsigned> loaded_kmer_lengths;
	std::unordered_map<unsigned, TranscriptMapping> kmer_to_transcripts;
	std::unordered_map<std::string, Transcript> transcripts;
	load_index(index_path,kmer_lengths,kmer_to_transcripts,transcripts);
	
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
	auto start = std::chrono::high_resolution_clock::now();
	auto homologous_segments = sparse_chain(read_sketches, kmer_to_transcripts, transcripts, reads, effective_kmer_lengths, 0.8);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
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
	std::vector<unsigned> kmer_lengths = { 31 };
	
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
