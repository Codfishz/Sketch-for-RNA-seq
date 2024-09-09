# RNA-seq Isoform Quantification using Sketching and Sparse Chaining
This repository contains an implementation of an RNA-seq Isoform Quantification tool, leveraging sketching techniques and sparse chaining to efficiently quantify gene isoforms from RNA-seq data. The goal of this project is to optimize the alignment and quantification process for RNA-seq reads using a combination of advanced hashing methods (such as FracMinhash) and sparse chaining algorithms.

Key Features
Efficient Isoform Quantification: Utilizes sketching techniques to reduce memory usage and accelerate the quantification process.
Sparse Chaining for Alignment: Implements sparse chaining to find homologous segments between RNA-seq reads and transcriptome reference sequences, reducing computational overhead.
FASTQ and FASTA Support: Reads RNA-seq data from FASTQ format and uses reference transcriptomes in FASTA format for efficient alignment.
Transcript Annotation Handling: Supports the use of GTF/GFF files for gene annotation to aid in the isoform quantification process.
TPM Calculation: Supports calculating TPM (Transcripts Per Million) to estimate gene expression levels.
Workflow Overview
Data Input: The user provides a reference genome in FASTA format and raw RNA-seq data in FASTQ format.
Quality Control: The RNA-seq reads undergo quality control (QC) before further processing.
Sketch Generation: K-mers are extracted from both the reference genome and the reads, and are hashed to generate sketches using techniques like FracMinhash.
Sparse Chaining: The tool uses sparse chaining to align the RNA-seq reads to the reference transcriptome by finding homologous segments.
Isoform Assignment: Reads are assigned to isoforms using gene annotation from GTF files, with read-isoform correction to ensure accurate assignment.
Isoform Expression Estimation: Isoform expression levels are calculated using the EM algorithm, with output in TPM.
Installation
Clone the repository:

bash
复制代码
git clone https://github.com/your-username/rna-seq-isoform-quantification.git
cd rna-seq-isoform-quantification
Requirements
C++ (C++17 or newer)
CMake
g++/clang++
SeqAn (for efficient sequence handling)
MurmurHash (for hashing k-mers)
Build Instructions
bash
复制代码
mkdir build && cd build
cmake ..
make
Usage
Run the tool with the following command:

bash
复制代码
./isoform_quantify -t 4 -m 2048 -k 31 -w 100 reference.fasta reads.fastq annotations.gtf output.csv
Options:

-t: Number of threads to use.
-m: Maximum memory to use (in MB).
-k: Length of k-mers to extract.
-w: Window size for sparse chaining.
Example Output
The tool outputs an output.csv file with isoform expression estimates, including columns like:

Transcript ID
TPM (Transcripts Per Million)
Read Counts
