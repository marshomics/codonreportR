#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os

def main():
    parser = argparse.ArgumentParser(description='Run hmmsearch on input files.')
    parser.add_argument('-i', '--input_files', required=True, nargs='+',
                        help='Input FASTA files to process')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for hmmsearch results')
    parser.add_argument('-hmm', '--hmm_file', required=True,
                        help='Path to the .hmm file')
    parser.add_argument('--cpu', type=int, default=1,
                        help='Number of CPU cores to use (default: 1)')

    args = parser.parse_args()

    input_files = args.input_files
    output_dir = args.output_dir
    hmm_file = args.hmm_file
    cpu_cores = args.cpu

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

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

if __name__ == '__main__':
    import shutil
    main()