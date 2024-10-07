#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import shutil
import pandas as pd
import concurrent.futures
from Bio import SeqIO

def run_hmmsearch(input_files, output_dir, hmm_file, cpu_cores):
    # Check if hmmsearch is available
    if not shutil.which('hmmsearch'):
        print("Error: 'hmmsearch' command not found. Please install HMMER3.", file=sys.stderr)
        sys.exit(1)

    for f in input_files:
        if not os.path.isfile(f):
            print(f"Input file does not exist: {f}", file=sys.stderr)
            continue

        # Construct the output file path
        assembly_name = os.path.basename(f)
        output_file = os.path.join(output_dir, f"{assembly_name}_hmmsearch.txt")

        # Build the hmmsearch command
        cmd = [
            'hmmsearch',
            '--cpu', str(cpu_cores),
            '--tblout', output_file,
            hmm_file,
            f
        ]

        # Run the hmmsearch command
        try:
            print(f"Running hmmsearch on {f}...")
            subprocess.run(cmd, check=True)
            print(f"Output saved to {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running hmmsearch on {f}: {e}", file=sys.stderr)

def combine_outputs(output_dir):
    # Define the header for the output file
    header = [
        "target_name",
        "accession1",
        "query_name",
        "accession2",
        "full_sequence_evalue",
        "full_sequence_bitscore",
        "full_sequence_bias",
        "best_1_domain_evalue",
        "best_1_domain_bitscore",
        "best_1_domain_bias",
        "exp",
        "reg",
        "clu",
        "ov",
        "env",
        "dom",
        "rep",
        "inc",
        "description",
        "genome"
    ]

    def process_file(file_path):
        data = []
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if not line.startswith("#"):
                        if line.strip() and "[No hits detected that satisfy reporting thresholds]" not in line:
                            genome = os.path.basename(file_path).replace("_hmmsearch.txt", "")
                            data.append(line.strip().split() + [genome])
        except Exception as e:
            print(f"Error processing file {file_path}: {e}", file=sys.stderr)
        return data

    def process_directory(files, root):
        result = []
        for file in files:
            if file.endswith("_hmmsearch.txt"):
                file_path = os.path.join(root, file)
                result.extend(process_file(file_path))
        return result

    # Initialize a list to store the combined data
    combined_data = []

    # Collect all directories and files to be processed
    tasks = []
    for root, dirs, files in os.walk(output_dir):
        if any(file.endswith("_hmmsearch.txt") for file in files):
            tasks.append((files, root))

    # Use ThreadPoolExecutor to parallelize file processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = []
        for files, root in tasks:
            futures.append(executor.submit(process_directory, files, root))

        for future in concurrent.futures.as_completed(futures):
            combined_data.extend(future.result())

    # Check if any data was collected
    if not combined_data:
        print("No data collected. Please check the output directory and files.", file=sys.stderr)
        sys.exit(1)

    # Convert the combined data to a DataFrame
    df = pd.DataFrame(combined_data, columns=header)
    return df

def filter_combined_output(df, bitscore_threshold):
    # Filter the DataFrame to retain only rows with "full_sequence_bitscore" >= bitscore_threshold
    df['full_sequence_bitscore'] = pd.to_numeric(df['full_sequence_bitscore'], errors='coerce')
    filtered_df = df[df['full_sequence_bitscore'] >= bitscore_threshold]
    return filtered_df

def split_fasta_files(input_dir, ribosomal_table_path, metadata_table_path,
                      output_dir_rp, output_dir_non_rp):
    # Create output directories if they don't exist
    os.makedirs(output_dir_rp, exist_ok=True)
    os.makedirs(output_dir_non_rp, exist_ok=True)

    # Load the ribosomal table
    ribosomal_table = pd.read_csv(ribosomal_table_path, sep='\t')

    # Extract the target_name column as a set for faster lookup
    target_names = set(ribosomal_table['target_name'])

    # Load the metadata table
    metadata_table = pd.read_csv(metadata_table_path, sep='\t')

    # Extract the Bin Id column as a set for faster lookup
    bin_ids = set(metadata_table['Bin Id'])

    # Process each FASTA file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna'):
            # Strip the file extension to get the base filename
            base_filename = os.path.splitext(filename)[0]

            # Check if the base filename is in the set of Bin Ids
            if base_filename in bin_ids:
                # Prepare the output file paths
                rp_output_path = os.path.join(output_dir_rp, f'{base_filename}_rp.fna')
                non_rp_output_path = os.path.join(output_dir_non_rp, f'{base_filename}_non_rp.fna')

                # Skip processing if both output files already exist
                if os.path.exists(rp_output_path) and os.path.exists(non_rp_output_path):
                    print(f'Skipping {filename}: Output files already exist')
                    continue

                input_file_path = os.path.join(input_dir, filename)
                rp_sequences = []
                non_rp_sequences = []

                # Read the sequences from the multi-fasta file
                for record in SeqIO.parse(input_file_path, 'fasta'):
                    if record.id in target_names:
                        rp_sequences.append(record)
                    else:
                        non_rp_sequences.append(record)

                # Write the ribosomal protein sequences to the new file
                if rp_sequences:
                    with open(rp_output_path, 'w') as rp_output_file:
                        SeqIO.write(rp_sequences, rp_output_file, 'fasta')

                # Write the non-ribosomal protein sequences to the new file
                if non_rp_sequences:
                    with open(non_rp_output_path, 'w') as non_rp_output_file:
                        SeqIO.write(non_rp_sequences, non_rp_output_file, 'fasta')

                print(f'Processed {filename}: {len(rp_sequences)} ribosomal, {len(non_rp_sequences)} non-ribosomal sequences')

    print('Processing complete.')

def main():
    parser = argparse.ArgumentParser(description='Run hmmsearch, combine outputs, filter, and split FASTA files.')
    parser.add_argument('-i', '--input_files', required=True, nargs='+',
                        help='Input FASTA files to process')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for hmmsearch results')
    parser.add_argument('-hmm', '--hmm_file', required=True,
                        help='Path to the .hmm file')
    parser.add_argument('--cpu', type=int, default=1,
                        help='Number of CPU cores to use (default: 1)')
    parser.add_argument('-c', '--combined_output', required=True,
                        help='Path for the combined hmmsearch output file')
    parser.add_argument('-f', '--filtered_output', required=True,
                        help='Path for the filtered combined output file')
    parser.add_argument('--bitscore', type=float, default=25.0,
                        help='Bitscore threshold for filtering (default: 25.0)')
    parser.add_argument('--split_input_dir', required=True,
                        help='Input directory containing FASTA files for splitting')
    parser.add_argument('--ribosomal_table', required=True,
                        help='Path to the ribosomal table (filtered combined output)')
    parser.add_argument('--metadata_table', required=True,
                        help='Path to the metadata table')
    parser.add_argument('--rp_output_dir', required=True,
                        help='Output directory for ribosomal protein sequences')
    parser.add_argument('--non_rp_output_dir', required=True,
                        help='Output directory for non-ribosomal protein sequences')

    args = parser.parse_args()

    input_files = args.input_files
    output_dir = args.output_dir
    hmm_file = args.hmm_file
    cpu_cores = args.cpu
    combined_output_file = args.combined_output
    filtered_output_file = args.filtered_output
    bitscore_threshold = args.bitscore
    split_input_dir = args.split_input_dir
    ribosomal_table_path = args.ribosomal_table
    metadata_table_path = args.metadata_table
    output_dir_rp = args.rp_output_dir
    output_dir_non_rp = args.non_rp_output_dir

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Run hmmsearch on input files
    run_hmmsearch(input_files, output_dir, hmm_file, cpu_cores)

    # Combine hmmsearch outputs
    print("Combining hmmsearch output files...")
    combined_df = combine_outputs(output_dir)

    # Save the combined DataFrame to a tab-separated file
    combined_df.to_csv(combined_output_file, sep='\t', index=False)
    print(f"Combined output saved to {combined_output_file}")

    # Filter the combined DataFrame
    print(f"Filtering combined output with bitscore >= {bitscore_threshold}...")
    filtered_df = filter_combined_output(combined_df, bitscore_threshold)

    # Save the filtered DataFrame to a tab-separated file
    filtered_df.to_csv(filtered_output_file, sep='\t', index=False)
    print(f"Filtered output saved to {filtered_output_file}")

    # Split FASTA files into ribosomal and non-ribosomal sequences
    print("Splitting FASTA files into ribosomal and non-ribosomal sequences...")
    split_fasta_files(split_input_dir, filtered_output_file, metadata_table_path,
                      output_dir_rp, output_dir_non_rp)

if __name__ == '__main__':
    main()