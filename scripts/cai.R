#!/usr/bin/env Rscript

library(optparse)
library(coRdon)
library(Biostrings)
library(IRanges)
library(data.table)

# Set up command-line options
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", help = "Input directory containing FASTA files"),
  make_option(c("-hmm", "--hmm_table_path"), type = "character", help = "Path to the HMM table"),
  make_option(c("-m", "--metadata_table_path"), type = "character", help = "Path to the metadata table"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory for CAI results")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$input_dir) || is.null(opt$hmm_table_path) || is.null(opt$metadata_table_path) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call. = FALSE)
}

# Assign command-line arguments to variables
input_dir <- opt$input_dir
hmm_table_path <- opt$hmm_table_path
metadata_table_path <- opt$metadata_table_path
output_dir <- opt$output_dir

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the tables
hmm_table <- fread(hmm_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata_table <- fread(metadata_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Get the list of FASTA files in the input directory
fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE)

# Extract the Bin Ids to process from the metadata table
bin_ids_to_process <- metadata_table$`Bin Id`

# Function to generate the output file name
get_output_file <- function(fasta_file) {
  output_file <- gsub("\\.(fasta|fa|fna)$", ".txt", basename(fasta_file))
  return(file.path(output_dir, output_file))
}

# Function to process each file
process_files <- function(fasta_file) {
  output_file <- get_output_file(fasta_file)
  
  # Skip processing if the output file already exists
  if (file.exists(output_file)) {
    message(paste("Output file already exists, skipping:", output_file))
    return()
  }
  
  # Extract the genome name from the file name
  genome_name <- gsub("\\.(fasta|fa|fna)$", "", basename(fasta_file))
  
  # Check if the genome name is in the list of Bin Ids to process
  if (!(genome_name %in% bin_ids_to_process)) {
    message(paste("Genome name not in Bin Ids to process, skipping:", genome_name))
    return()
  }
  
  # Get the rows matching the genome name
  matching_rows <- hmm_table[hmm_table$genome == genome_name, ]
  
  # Check if there are matching rows
  if (nrow(matching_rows) == 0) {
    message(paste("No matching rows found for:", genome_name))
    return()
  }
  
  # Create the rp list object
  rp <- list(rp = matching_rows$target_name)
  
  tryCatch({
    # Read the FASTA file
    fasta <- readSet(file = fasta_file)
    
    # Generate codon table
    codons <- codonTable(fasta)
    
    # Add annotations
    codons@KO <- codons@ID
    
    # Perform CAI analysis
    cai <- CAI(codons, filtering = "none", subsets = rp, id_or_name2 = "11")
    cai_df <- as.data.frame(cai)
    
    # Extract names and widths, then create data frames
    names_df <- data.frame(gene = names(fasta))
    width_df <- data.frame(width = width(fasta))
    
    # Combine CAI results with names and widths
    cai_ann <- cbind(cai_df, names_df, width_df)
    filtered_cai_ann <- subset(cai_ann, width > 240)
    
    # Write the combined data frame to a file
    write.table(filtered_cai_ann, file = output_file, sep = "\t", row.names = FALSE)
    message(paste("Processed and saved CAI results for:", genome_name))
  }, error = function(e) {
    message(paste("Error occurred while processing file:", fasta_file))
    message("Error message:", e$message)
  })
}

# Process each FASTA file
lapply(fasta_files, process_files)