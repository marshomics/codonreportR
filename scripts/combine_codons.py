#!/usr/bin/env python3

import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='Combine codon count files into a single file.')
    parser.add_argument('-i', '--input_directory', required=True, help='Input directory containing codon count files')
    parser.add_argument('-o', '--output_file', required=True, help='Output file path for the combined codon usage')

    args = parser.parse_args()

    input_directory = args.input_directory
    output_file = args.output_file

    # Open the output file in write mode
    with open(output_file, 'w') as outfile:
        header_written = False

        # Loop through all the files in the input directory
        for filename in os.listdir(input_directory):
            # Only consider *.txt files that do not end with "_normalized.txt"
            if filename.endswith(".txt") and not filename.endswith("_normalized.txt"):
                file_path = os.path.join(input_directory, filename)

                with open(file_path, 'r') as infile:
                    lines = infile.readlines()

                    # Write the header only once
                    if not header_written and lines:
                        outfile.write(lines[0].strip() + "\tassembly\n")
                        header_written = True

                    # Write the rest of the file, adding the assembly column
                    assembly_name = filename[:-4]
                    for line in lines[1:]:
                        outfile.write(line.strip() + f"\t{assembly_name}\n")

    print(f"Combined file saved to {output_file}")

if __name__ == '__main__':
    main()