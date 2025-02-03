#ifndef DATA_IO_H
#define DATA_IO_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map> 

// 结构体定义
struct Transcript {
    std::string id;
    std::string sequence;
    int length;
};

// struct Transcript {
//     std::string id;
//     std::string sequence;
//     std::string gene_id;
//     int start;
//     int end;
//     int length;
// };

// 基因注释信息结构
struct GeneAnnotation {
    std::string gene_id;
    std::string transcript_id;
    int start;
    int end;
    int length;
};

struct Read {
    std::string id;
    std::string sequence;
    std::string quality;
};

// 函数声明

std::unordered_map<std::string, Transcript> load_fasta(const std::string& fasta_file); 
std::unordered_map<std::string, Read> load_fastq(const std::string& fastq_file);
void output_to_csv(const std::string& filename, 
                   const std::unordered_map<std::string, double>& read_counts,
                   const std::unordered_map<std::string, double>& tpms,
                   const std::unordered_map<std::string, double>& pi,
                   const std::unordered_map<std::string, Transcript>& transcripts)  ;

// std::vector<Read> quality_control(const std::vector<Read>& reads);

#endif // DATA_IO_H
