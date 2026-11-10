# Bacterial DNA Sequence Analysis

## Project Overview

This project focuses on the analysis of bacterial 16S ribosomal RNA (16S rRNA) DNA sequences using R programming.

Six bacterial 16S rRNA sequences were downloaded from the NCBI database and analyzed to explore their nucleotide composition, GC content, and sequence alignment.

## Dataset

- Source: NCBI (National Center for Biotechnology Information)
- Data type: Bacterial 16S rRNA sequences
- Number of sequences: 6
- File format: FASTA (`.fasta`)

## Analysis Performed

### 1. Sequence Import
- Imported bacterial DNA sequences into R.
- Examined the sequence data.

### 2. Nucleotide Composition Analysis
- Calculated the number of Adenine (A), Thymine (T), Guanine (G), and Cytosine (C) nucleotides.
- Compared nucleotide composition across the six bacterial sequences.

### 3. GC Content Analysis
- Calculated the GC content percentage for each sequence.
- Compared GC content among the bacterial sequences.

### 4. Data Visualization
Generated visualizations to present the analysis results:
- Nucleotide composition comparison.
- GC content comparison.

### 5. Multiple Sequence Alignment
- Performed multiple sequence alignment using the DECIPHER package in R.
- Saved the aligned bacterial DNA sequences in FASTA format.

## Technologies and Tools

- R Programming Language
- RStudio
- DECIPHER package
- NCBI database

## Project Structure

```text
Bacterial_DNA_Analysis/
│
├── code/
│   └── 01_bacterial_sequence_analysis.R
│
├── data/
│   └── bacterial_16S_sequences.fasta
│
├── results/
│   ├── nucleotide_analysis.csv
│   ├── nucleotide_composition.png
│   ├── gc_content_comparison.png
│   └── aligned_sequences.fasta
│
├── report/
│
└── README.md
```

## Project Status

Completed the initial sequence import, nucleotide composition analysis, GC content analysis, visualization, and multiple sequence alignment stages.