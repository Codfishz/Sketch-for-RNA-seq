# RNA-seq Isoform Quantification using Sketching and Sparse Chaining

**RNA-seq Isoform Quantification** is a high-effi command-line tool designed to quantify gene isoforms from RNA-seq data. This tool leverages advanced sketching techniques, including fractional MinHash (FracMinhash), together with sparse chaining algorithms to improve the accuracy and speed of alignment and quantification tasks.

## Key Features
- **Efficient Isoform Quantification**  
  Uses innovative sketching methods to reduce memory usage and accelerate RNA-seq analysis.
- **Sparse Chaining for Alignment**  
  Implements sparse chaining to identify homologous segments between RNA-seq reads and the reference transcriptome with reduced computational overhead.
- **Support for FASTQ and FASTA Formats**  
  Reads RNA-seq data from FASTQ files and aligns against transcriptomes provided in FASTA format.
- **TPM Calculation**  
  Calculates Gene Abundance to standardize gene expression estimates across samples.

## Workflow Overview
1. **Data Input**  
   The user supplies a reference transcriptome (FASTA format) and raw RNA-seq reads (FASTQ format).  
2. **Sketch Generation**  
   K-mers are extracted from both the reference and the reads and are hashed using FracMinhash to create compact sketches.
4. **Sparse Chaining**  
   Homologous segments between RNA-seq reads and transcripts are identified using a sparse chaining algorithm, drastically reducing the search space.
5. **Isoform Assignment**  
   Reads are assigned to transcripts, with additional correction steps for improved accuracy.
6. **Isoform Expression Estimation**  
   The EM algorithm is employed to estimate isoform abundances, with the final output given in TPM.

## Installation

### Prerequisites
- **C++ Compiler**: A C++17 compliant compiler (e.g., g++).
- **nthash Library**: Used for fast ks-mer hashing. Please install it from the [nthash GitHub repository](https://github.com/BIMSBbioinfo/nthash).
- **Bash**: For running the provided build and testing scripts.
- **Unix-like Environment**: Linux/macOS (adapt instructions if targeting Windows).

### Building the Project
Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/your-username/rna-seq-isoform-quantification.git
cd rna-seq-isoform-quantification

## Usage

The tool supports two primary modes of operation: **index** and **quant**.

### Command-Line Options

- `-h, --help`  
  Display the help message and exit.

- `-k, --kmer-length SIZE`  
  Provide a comma-separated list of k-mer lengths to be used (default: 31).

- `-o, --mode MODE`  
  Specify the mode of operation. Available options are:
  - `index`: Build an index from a reference transcriptome.
  - `quant`: Quantify isoform expression using a pre-built index and RNA-seq reads (default).

### Examples

#### 1. Indexing Mode

Build an index from your reference transcriptome (FASTA file) and save the output to a specified index file:

```bash
./build/test -o index /path/to/reference_genome.fasta /path/to/index_output

#### 2. Quantification Mode

Use a pre-built index along with RNA-seq reads (FASTQ file) to quantify isoform expression. The results will be saved in a CSV file:

```bash
./build/test -o quant /path/to/index_file /path/to/reads.fastq /path/to/output.csv


## Directory Structure

The repository is organized as follows:

```
rna-seq-isoform-quantification/
├── build/               # Build output (compiled executable, analysis reports, etc.)
├── include/             # Header files for the project
├── src/                 # Source code files
├── Test_Data/           # Sample data for testing (FASTA, FASTQ, GTF/GFF, etc.)
├── scripts/             # Bash scripts for building, running tests, and accuracy evaluation
├── LICENSE              # License file (e.g., MIT License)
├── README.md            # This README file
└── CODE_OF_CONDUCT.md   # Code of Conduct for contributors
```

This layout helps you quickly understand the project's structure and find the resources and tools you need.
```

You can copy this snippet directly into your README file under a section titled **"Directory Structure"**. This will provide a clear and professional overview of your repository's organization.