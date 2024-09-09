# RNA-seq Isoform Quantification using Sketching and Sparse Chaining

This repository contains an implementation of an RNA-seq Isoform Quantification tool, leveraging sketching techniques and sparse chaining to efficiently quantify gene isoforms from RNA-seq data. The goal of this project is to optimize the alignment and quantification process for RNA-seq reads using a combination of advanced hashing methods (such as FracMinhash) and sparse chaining algorithms.

## Key Features
- **Efficient Isoform Quantification**: Utilizes sketching techniques to reduce memory usage and accelerate the quantification process.
- **Sparse Chaining for Alignment**: Implements sparse chaining to find homologous segments between RNA-seq reads and transcriptome reference sequences, reducing computational overhead.
- **FASTQ and FASTA Support**: Reads RNA-seq data from FASTQ format and uses reference transcriptomes in FASTA format for efficient alignment.
- **Transcript Annotation Handling**: Supports the use of GTF/GFF files for gene annotation to aid in the isoform quantification process.
- **TPM Calculation**: Supports calculating TPM (Transcripts Per Million) to estimate gene expression levels.

## Workflow Overview
1. **Data Input**: The user provides a reference genome in FASTA format and raw RNA-seq data in FASTQ format.
2. **Quality Control**: The RNA-seq reads undergo quality control (QC) before further processing.
3. **Sketch Generation**: K-mers are extracted from both the reference genome and the reads, and are hashed to generate sketches using techniques like FracMinhash.
4. **Sparse Chaining**: The tool uses sparse chaining to align the RNA-seq reads to the reference transcriptome by finding homologous segments.
5. **Isoform Assignment**: Reads are assigned to isoforms using gene annotation from GTF files, with read-isoform correction to ensure accurate assignment.
6. **Isoform Expression Estimation**: Isoform expression levels are calculated using the EM algorithm, with output in TPM.

## Installation
Clone the repository:
```bash
git clone https://github.com/your-username/rna-seq-isoform-quantification.git
cd rna-seq-isoform-quantification
