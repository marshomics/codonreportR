# CodonReporter

CodonReporter is a comprehensive bioinformatics tool designed to analyze codon usage patterns, identify ribosomal protein genes, and compute metrics like Codon Adaptation Index (CAI) and Relative Synonymous Codon Usage (RSCU). It integrates various scripts and programs to provide a streamlined workflow for genomic data analysis.

---

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [`usage` Command](#usage-command)
  - [`ribosomal` Command](#ribosomal-command)
  - [`cai` Command](#cai-command)
  - [`rscu` Command](#rscu-command)
- [Outputs](#outputs)
- [Examples](#examples)
- [Dependencies](#dependencies)
- [License](#license)

---

## Features

- **Codon Usage Analysis**: Compute codon counts and normalize them across multiple FASTA files.
- **Ribosomal Protein Analysis**: Identify ribosomal protein genes using HMMER's `hmmsearch` and combine results.
- **Codon Adaptation Index (CAI) Calculation**: Compute CAI values for genes using ribosomal proteins as a reference.
- **Relative Synonymous Codon Usage (RSCU) Calculation**: Compute RSCU values for genes and calculate median values for both ribosomal and non-ribosomal genes.
- **Flexible Workflow**: Modular commands allow for customized analysis pipelines.
- **Parallel Processing**: Utilizes multiprocessing to handle large datasets efficiently.

---

## Prerequisites

- **Operating System**: Unix/Linux-based system recommended.
- **Python**: Version 3.6 or higher.
- **R**: Version 3.5 or higher.

---

## Installation

1. **Clone the Repository**:

   ```
   git clone https://github.com/yourusername/codonreporter.git
   cd codonreporter
   ```

2.	**Set Up the Environment**:
	•	Ensure that python3 and Rscript are available in your PATH.
	•	Install required Python packages:

    ```
    pip install pandas
    ```

	•	Install required R packages by running R and executing:

    ```R
    install.packages(c("optparse", "data.table"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("coRdon", "Biostrings", "IRanges"))
    ```

3.	**Set Execution Permissions**:

    ```
    chmod +x codonreporter
    chmod +x scripts/*.R
    chmod +x scripts/*.py
    chmod +x programs/Compute_RSCU_gene.pyz
    ```

## Usage

The codonreporter script is the main entry point and supports multiple sub-commands:

    ```
    ./codonreporter [command] [options]
    ```

### `usage` Command

Description: Performs codon usage analysis on a set of FASTA files.

Arguments:

	•	-m, --metadata_file: Path to the metadata file containing filenames to process.
	•	-f, --fasta_dir: Directory containing input FASTA files.
	•	-o, --output_dir: Output directory for codon counts.
	•	-c, --combined_output: Path to the combined codon usage output file.

### `ribosomal` Command

Description: Identifies ribosomal protein genes using hmmsearch, combines outputs, filters results, and splits sequences.

Arguments:

	•	-i, --input_files: Input FASTA files to process (accepts multiple files).
	•	-o, --output_dir: Output directory for hmmsearch results.
	•	-hmm, --hmm_file: Path to the .hmm file.
	•	--cpu: Number of CPU cores to use (default: 1).
	•	-c, --combined_output: Path for the combined hmmsearch output file.
	•	-f, --filtered_output: Path for the filtered combined output file.
	•	--bitscore: Bitscore threshold for filtering (default: 25.0).

### `cai` Command

Description: Calculates the Codon Adaptation Index (CAI) for genes using ribosomal proteins as a reference.

Arguments:

	•	-i, --input_dir: Input directory containing FASTA files.
	•	-hmm, --hmm_table_path: Path to the filtered HMM table.
	•	-m, --metadata_table_path: Path to the metadata table.
	•	-o, --output_dir: Output directory for CAI results.

### `rscu` Command

Description: Computes RSCU values and calculates median values for both ribosomal and non-ribosomal protein genes.

Arguments:

	•	-ri, --ribosomal_input_dir: Input directory with ribosomal protein FASTA files.
	•	-ni, --non_ribosomal_input_dir: Input directory with non-ribosomal protein FASTA files.
	•	-ro, --ribosomal_output_dir: Output directory for ribosomal RSCU results.
	•	-no, --non_ribosomal_output_dir: Output directory for non-ribosomal RSCU results.
	•	-rm, --ribosomal_median_output: Output file for ribosomal median values.
	•	-nm, --non_ribosomal_median_output: Output file for non-ribosomal median values.

## Outputs

### `usage` Command Outputs

	•	Individual Codon Counts: For each input FASTA file, generates codon count files.
	•	Combined Codon Usage File: A combined file aggregating codon counts from all processed files.

### `ribosomal` Command Outputs

	•	HMMsearch Results: Individual hmmsearch output files for each input FASTA file.
	•	Combined HMMsearch Output: A single file combining all hmmsearch results.
	•	Filtered Output: A file containing filtered results based on bitscore threshold.
	•	Split Sequences: FASTA files split into ribosomal and non-ribosomal sequences.

### `cai` Command Outputs

	•	CAI Results: CAI values calculated for genes, outputted for each input FASTA file.

### `rscu` Command Outputs

	•	RSCU CSV Files: RSCU values computed for each gene in the input FASTA files.
	•	Median Values: Tab-separated files containing median RSCU values for both ribosomal and non-ribosomal genes.

## Dependencies

	Python Packages:
	•	pandas
	•	argparse

Install using:

    ```
    pip install pandas argparse
    ```

	R Packages:
	•	optparse
	•	data.table
	•	coRdon
	•	Biostrings
	•	IRanges

    Install using R:

    ```
    install.packages(c("optparse", "data.table"))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("coRdon", "Biostrings", "IRanges"))
    ```

	External Tools:
	•	HMMER3: hmmsearch must be installed and accessible in your PATH.

    ```
    sudo apt-get install hmmer
    ```

    Compute_RSCU_gene.pyz: Place in programs/ directory.
