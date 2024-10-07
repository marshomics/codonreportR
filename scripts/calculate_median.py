#!/usr/bin/env python3

import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import argparse

def process_file(filename, input_directory, columns_of_interest):
    file_path = os.path.join(input_directory, filename)
    
    # Read only the necessary columns from the CSV file
    df = pd.read_csv(file_path, usecols=columns_of_interest)
    
    # Calculate the median values for the columns of interest
    medians = df.median()
    
    # Return the filename (without .csv extension) and the median values as a list
    return [filename[:-4]] + medians.tolist()

def save_results(results, output_file, columns_of_interest):
    # Convert the list of results to a DataFrame
    median_df = pd.DataFrame(results, columns=["Filename"] + columns_of_interest)
    
    # Save the DataFrame to a tab-separated file
    median_df.to_csv(output_file, sep='\t', index=False, mode='a', header=not os.path.exists(output_file))

def main():
    parser = argparse.ArgumentParser(description='Calculate median values for each gene.')
    parser.add_argument('-i', '--input_directory', required=True, help='Input directory containing RSCU CSV files')
    parser.add_argument('-o', '--output_file', required=True, help='Output file path for the median values')
    args = parser.parse_args()

    input_directory = args.input_directory
    output_file = args.output_file

    # Define the columns of interest
    columns_of_interest = [
        "Phe-UUU", "Phe-UUC", "Ser4-UCU", "Ser4-UCC", "Ser4-UCA", "Ser4-UCG", "Ser2-AGU", "Ser2-AGC",
        "Leu4-CUU", "Leu4-CUC", "Leu4-CUA", "Leu4-CUG", "Leu2-UUA", "Leu2-UUG", "Tyr-UAU", "Tyr-UAC",
        "Cys-UGU", "Cys-UGC", "Arg4-CGU", "Arg4-CGC", "Arg4-CGA", "Arg4-CGG", "Arg2-AGA", "Arg2-AGG",
        "Pro-CCU", "Pro-CCC", "Pro-CCA", "Pro-CCG", "His-CAU", "His-CAC", "Gln-CAA", "Gln-CAG",
        "Ile-AUU", "Ile-AUC", "Ile-AUA", "Thr-ACU", "Thr-ACC", "Thr-ACA", "Thr-ACG",
        "Asn-AAU", "Asn-AAC", "Lys-AAA", "Lys-AAG", "Val-GUU", "Val-GUC", "Val-GUA", "Val-GUG",
        "Ala-GCU", "Ala-GCC", "Ala-GCA", "Ala-GCG", "Asp-GAU", "Asp-GAC", "Glu-GAA", "Glu-GAG",
        "Gly-GGU", "Gly-GGC", "Gly-GGA", "Gly-GGG"
    ]

    # Use ProcessPoolExecutor to process files in parallel
    with ProcessPoolExecutor() as executor:
        # List all CSV files in the directory
        csv_files = [f for f in os.listdir(input_directory) if f.endswith(".csv")]
        
        # Chunk size for processing
        chunk_size = 1000
        results = []
        
        # Process files in parallel and batch results to save
        for i in range(0, len(csv_files), chunk_size):
            chunk = csv_files[i:i + chunk_size]
            futures = [executor.submit(process_file, filename, input_directory, columns_of_interest) for filename in chunk]
            for future in futures:
                results.append(future.result())
            
            # Save results periodically to reduce memory usage
            save_results(results, output_file, columns_of_interest)
            results.clear()  # Clear results to free up memory

    print(f"Median values saved to {output_file}")

if __name__ == "__main__":
    main()