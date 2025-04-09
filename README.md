```markdown
# RNA-seq Isoform Quantification using Sketching and Sparse Chaining

**RNA-seq Isoform Quantification** is a high-efficiency command-line tool designed to quantify gene isoforms from RNA-seq data. This tool leverages advanced sketching techniques, including fractional MinHash (FracMinhash), together with sparse chaining algorithms to improve the accuracy and speed of alignment and quantification tasks.

---

## Key Features
- **Efficient Isoform Quantification**  
  Uses innovative sketching methods to reduce memory usage and accelerate RNA-seq analysis.
- **Sparse Chaining for Alignment**  
  Implements sparse chaining to identify homologous segments between RNA-seq reads and the reference transcriptome with reduced computational overhead.
- **Support for FASTQ and FASTA Formats**  
  Reads RNA-seq data from FASTQ files and aligns against transcriptomes provided in FASTA format.
- **TPM Calculation**  
  Calculates gene abundance to standardize gene expression estimates across samples.

---

## Workflow Overview
1. **Data Input**  
   The user supplies a reference transcriptome (FASTA format) and raw RNA-seq reads (FASTQ format).
2. **Sketch Generation**  
   K-mers are extracted from both the reference and the reads and are hashed using FracMinhash to create compact sketches.
3. **Sparse Chaining**  
   Homologous segments between RNA-seq reads and transcripts are identified using a sparse chaining algorithm, drastically reducing the search space.
4. **Isoform Assignment**  
   Reads are assigned to transcripts, with additional correction steps for improved accuracy.
5. **Isoform Expression Estimation**  
   The EM algorithm is employed to estimate isoform abundances, with the final output given in TPM.

---

## Installation

### Prerequisites
- **C++ Compiler**: A C++17 compliant compiler (e.g., g++).
- **nthash Library**: Used for fast k-mer hashing. Please install it from the [nthash GitHub repository](https://github.com/BIMSBbioinfo/nthash).
- **Bash**: For running the provided build and testing scripts.
- **Unix-like Environment**: Linux/macOS (adapt instructions if targeting Windows).

### Building the Project
Clone the repository and navigate to the project directory:
```bash
git clone https://github.com/your-username/rna-seq-isoform-quantification.git
cd rna-seq-isoform-quantification
```
Build the project using the provided `build.sh` script:
```bash
chmod +x build.sh
./build.sh
```
This script will compile all source files in the `src` directory, link against the nthash library, and output the executable (by default named `test`) into the `build/` directory.

---

## Usage

The tool supports two primary modes of operation: **index** and **quant**.

### Command-Line Options

- **`-h, --help`**  
  Display the help message and exit.

- **`-k, --kmer-length SIZE`**  
  Provide a comma-separated list of k-mer lengths to be used (default: 31).

- **`-o, --mode MODE`**  
  Specify the mode of operation. Available options are:
  - `index`: Build an index from a reference transcriptome.
  - `quant`: Quantify isoform expression using a pre-built index and RNA-seq reads (default).

### Examples

#### 1. Indexing Mode

Build an index from your reference transcriptome (FASTA file) and save the output to a specified index file:
```bash
./build/test -o index /path/to/reference_genome.fasta /path/to/index_output
```

#### 2. Quantification Mode

Use a pre-built index along with RNA-seq reads (FASTQ file) to quantify isoform expression. The results will be saved in a CSV file:
```bash
./build/test -o quant /path/to/index_file /path/to/reads.fastq /path/to/output.csv
```

> **Note:**  
> If you built the project with profiling enabled (using `-g -pg`), after running the tool you can generate a profiling report using:
> ```bash
> gprof ./build/test gmon.out > build/analysis.txt
> ```

---

## Directory Structure

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

This layout helps you quickly understand the project's structure and locate the necessary resources and tools.

---



## Acknowledgements

- **nthash Library**: Thanks to the developers of nthash for providing a fast k-mer hashing solution.

_Last updated on YYYY-MM-DD._
```
