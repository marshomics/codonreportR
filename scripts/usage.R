#!/usr/bin/env Rscript

# Load required libraries
library(coRdon)
library(Biostrings)
library(optparse)

# Set up command-line options
option_list <- list(
  make_option(c("-m", "--metadata_file"), type="character", default=NULL, 
              help="Path to the metadata file", metavar="FILE"),
  make_option(c("-f", "--fasta_dir"), type="character", default=NULL, 
              help="Directory containing the FASTA files", metavar="DIR"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for codon counts", metavar="DIR")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$metadata_file) || is.null(opt$fasta_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

# Read the table of filenames to process
metadata <- read.table(opt$metadata_file, sep = "\t", header = TRUE)

# Directory with the FASTA files
fasta_dir <- opt$fasta_dir

# Output directory for codon counts
output_dir <- opt$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Numeric columns to be divided by the length column
numeric_columns <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
                     "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
                     "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
                     "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
                     "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
                     "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
                     "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
                     "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")

# Loop through each filename in the metadata table
for (filename in metadata$filename) {
  fasta_file <- file.path(fasta_dir, paste0(filename, ".fasta"))
  
  # Check if the FASTA file exists
  if (file.exists(fasta_file)) {
    # Read the FASTA file
    fasta_seqs <- readDNAStringSet(fasta_file)
    
    # Extract headers and sequence lengths
    fasta_headers <- names(fasta_seqs)
    fasta_width <- width(fasta_seqs)
    
    # Read the genome set
    genome <- readSet(file = fasta_file)
    
    # Compute the codon usage table
    genome_codon <- codonTable(genome)
    genome_counts <- codonCounts(genome_codon)
    genome_counts_df <- as.data.frame(genome_counts)
    rownames(genome_counts_df) <- fasta_headers
    genome_counts_df$length <- fasta_width
    
    # Write the original codon counts to a file
    output_file <- file.path(output_dir, paste0(filename, ".txt"))
    write.table(genome_counts_df, file = output_file, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    
    # Normalize specified numeric columns by the length column
    normalized_df <- genome_counts_df
    for (col in numeric_columns) {
      normalized_df[[col]] <- normalized_df[[col]] / normalized_df$length
    }
    
    # Write the normalized codon counts to a second output file
    normalized_output_file <- file.path(output_dir, paste0(filename, "_normalized.txt"))
    write.table(normalized_df, file = normalized_output_file, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
    
  } else {
    cat("File does not exist:", fasta_file, "\n")
  }
}