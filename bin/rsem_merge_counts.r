#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")

# R function to merge RSEM counts from multiple samples with annotation integration
option_list <- list(
    make_option(c("-g", "--genes_dir"), type="character", default=NULL, help="Directory containing gene count files", metavar="character"),
    make_option(c("-t", "--transcripts_dir"), type="character", default=NULL, help="Directory containing transcript count files", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file prefix", metavar="character"),
    make_option(c("-a", "--annotation_matrix"), type="character", default=NULL, help="Annotation matrix file path", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Load annotation matrix if provided
annotation_data <- NULL
if (!is.null(opt$annotation_matrix) && file.exists(opt$annotation_matrix)) {
    cat("Loading annotation matrix from:", opt$annotation_matrix, "\n")
    full_annotation <- read.table(opt$annotation_matrix, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cat("Loaded annotation for", nrow(full_annotation), "genes\n")
    
    # Filter annotation to exclude specific boolean columns:
    # Exclude: has_transcript, has_exon, has_intron, has_5utr, has_3utr
    # Keep: gene_id, chromosome, strand, plus any additional metadata columns
    
    # Define columns to exclude
    excluded_cols <- c("has_transcript", "has_exon", "has_intron", "has_5utr", "has_3utr")
    
    # Get all available columns and exclude the unwanted ones
    available_cols <- colnames(full_annotation)
    selected_cols <- setdiff(available_cols, excluded_cols)
    
    annotation_data <- full_annotation[, selected_cols, drop = FALSE]
    
    cat("Selected annotation columns:", paste(colnames(annotation_data), collapse = ", "), "\n")
    cat("Filtered annotation to", ncol(annotation_data), "columns (", ncol(annotation_data) - 1, "annotation columns)\n")
}

# Process gene counts
cat("Processing gene count files...\n")
gene_files <- list.files(opt$genes_dir, pattern = "\\.genes\\.results$", full.names = TRUE)

if (length(gene_files) == 0) {
    stop("No gene count files found in:", opt$genes_dir)
}

cat("Found", length(gene_files), "gene count files\n")

# Read first file to get gene IDs and gene names
first_gene_file <- read.table(gene_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_info <- first_gene_file[, c("gene_id", "transcript_id.s.")]
colnames(gene_info) <- c("gene_id", "gene_name")

# Initialize count and TPM matrices
gene_counts <- data.frame(gene_id = gene_info$gene_id, gene_name = gene_info$gene_name)
gene_tpm <- data.frame(gene_id = gene_info$gene_id, gene_name = gene_info$gene_name)

# Process each gene file
for (gene_file in gene_files) {
    sample_name <- gsub("\\.genes\\.results$", "", basename(gene_file))
    cat("Processing gene file:", sample_name, "\n")
    
    gene_data <- read.table(gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Extract counts and TPM
    gene_counts[[sample_name]] <- gene_data$expected_count
    gene_tpm[[sample_name]] <- gene_data$TPM
}

# Process transcript counts
cat("Processing transcript count files...\n")
transcript_files <- list.files(opt$transcripts_dir, pattern = "\\.isoforms\\.results$", full.names = TRUE)

if (length(transcript_files) == 0) {
    stop("No transcript count files found in:", opt$transcripts_dir)
}

cat("Found", length(transcript_files), "transcript count files\n")

# Read first file to get transcript IDs and gene IDs
first_transcript_file <- read.table(transcript_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
transcript_info <- first_transcript_file[, c("transcript_id", "gene_id")]

# Initialize transcript count and TPM matrices
transcript_counts <- data.frame(transcript_id = transcript_info$transcript_id, gene_id = transcript_info$gene_id)
transcript_tpm <- data.frame(transcript_id = transcript_info$transcript_id, gene_id = transcript_info$gene_id)

# Process each transcript file
for (transcript_file in transcript_files) {
    sample_name <- gsub("\\.isoforms\\.results$", "", basename(transcript_file))
    cat("Processing transcript file:", sample_name, "\n")
    
    transcript_data <- read.table(transcript_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Extract counts and TPM
    transcript_counts[[sample_name]] <- transcript_data$expected_count
    transcript_tpm[[sample_name]] <- transcript_data$TPM
}

# Create output directory
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Merge gene counts with annotation data if available
final_gene_counts <- gene_counts
final_gene_tpm <- gene_tpm

if (!is.null(annotation_data)) {
    # Merge gene counts with annotation information
    final_gene_counts <- merge(annotation_data, gene_counts, 
                              by = "gene_id", all.y = TRUE, all.x = FALSE)
    final_gene_tpm <- merge(annotation_data, gene_tpm, 
                           by = "gene_id", all.y = TRUE, all.x = FALSE)
    cat("Merged gene counts and TPM with annotation data\n")
}

# Save gene results
write.table(final_gene_counts, 
           file = paste0(opt$output, ".merged.gene_counts.tsv"), 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

write.table(final_gene_tpm, 
           file = paste0(opt$output, ".merged.gene_tpm.tsv"), 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

# Save transcript results (no annotation merge for transcripts as annotation is gene-level)
write.table(transcript_counts, 
           file = paste0(opt$output, ".merged.transcript_counts.tsv"), 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

write.table(transcript_tpm, 
           file = paste0(opt$output, ".merged.transcript_tpm.tsv"), 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

cat("RSEM merge completed successfully!\n")
if (!is.null(annotation_data)) {
    cat("Gene-level files include annotation with", ncol(annotation_data) - 1, "annotation columns\n")
    cat("Annotation columns included:", paste(colnames(annotation_data)[-1], collapse = ", "), "\n")
    cat("Note: Excluded has_transcript, has_exon, has_intron, has_5utr, has_3utr columns. All other annotation and metadata columns are included\n")
}

cat("Output files:\n")
cat("- Gene counts:", paste0(opt$output, ".merged.gene_counts.tsv"), "\n")
cat("- Gene TPM:", paste0(opt$output, ".merged.gene_tpm.tsv"), "\n")
cat("- Transcript counts:", paste0(opt$output, ".merged.transcript_counts.tsv"), "\n")
cat("- Transcript TPM:", paste0(opt$output, ".merged.transcript_tpm.tsv"), "\n")