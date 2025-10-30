#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")

# R function to merge tximport counts from multiple samples with annotation integration
option_list <- list(
    make_option(c("--gene_counts"), type="character", default=NULL, help="Gene counts file from tximport", metavar="character"),
    make_option(c("--gene_tpm"), type="character", default=NULL, help="Gene TPM file from tximport", metavar="character"),
    make_option(c("--gene_counts_length_scaled"), type="character", default=NULL, help="Gene counts length scaled file from tximport", metavar="character"),
    make_option(c("--gene_counts_scaled"), type="character", default=NULL, help="Gene counts scaled file from tximport", metavar="character"),
    make_option(c("--gene_lengths"), type="character", default=NULL, help="Gene lengths file from tximport", metavar="character"),
    make_option(c("--transcript_counts"), type="character", default=NULL, help="Transcript counts file from tximport", metavar="character"),
    make_option(c("--transcript_tpm"), type="character", default=NULL, help="Transcript TPM file from tximport", metavar="character"),
    make_option(c("--transcript_lengths"), type="character", default=NULL, help="Transcript lengths file from tximport", metavar="character"),
    make_option(c("--annotation_matrix"), type="character", default=NULL, help="Annotation matrix file path", metavar="character"),
    make_option(c("--output_prefix"), type="character", default="tximport", help="Output file prefix", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=== TXIMPORT MERGE WITH ANNOTATION ===\n")

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

# Function to process and merge a tximport file with annotation
process_tximport_file <- function(file_path, file_type, annotation_data) {
    if (is.null(file_path) || !file.exists(file_path)) {
        cat("Warning: File", file_type, "not found or NULL, skipping...\n")
        return(NULL)
    }
    
    cat("Processing", file_type, "file:", file_path, "\n")
    
    # Read the tximport file
    data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Check if gene_id column exists, if not use first column as gene_id
    if (!"gene_id" %in% colnames(data)) {
        if (ncol(data) > 0) {
            colnames(data)[1] <- "gene_id"
            cat("Renamed first column to 'gene_id' for", file_type, "\n")
        } else {
            cat("Warning: Empty data file for", file_type, ", skipping...\n")
            return(NULL)
        }
    }
    
    # Merge with annotation data if available
    if (!is.null(annotation_data)) {
        merged_data <- merge(annotation_data, data, by = "gene_id", all.y = TRUE, all.x = FALSE)
        cat("Merged", file_type, "with annotation data\n")
        return(merged_data)
    } else {
        return(data)
    }
}

# Process each tximport file type
file_types <- list(
    list(path = opt$gene_counts, type = "gene_counts", suffix = "gene_counts.tsv"),
    list(path = opt$gene_tpm, type = "gene_tpm", suffix = "gene_tpm.tsv"),
    list(path = opt$gene_counts_length_scaled, type = "gene_counts_length_scaled", suffix = "gene_counts_length_scaled.tsv"),
    list(path = opt$gene_counts_scaled, type = "gene_counts_scaled", suffix = "gene_counts_scaled.tsv"),
    list(path = opt$gene_lengths, type = "gene_lengths", suffix = "gene_lengths.tsv"),
    list(path = opt$transcript_counts, type = "transcript_counts", suffix = "transcript_counts.tsv"),
    list(path = opt$transcript_tpm, type = "transcript_tpm", suffix = "transcript_tpm.tsv"),
    list(path = opt$transcript_lengths, type = "transcript_lengths", suffix = "transcript_lengths.tsv")
)

# Initialize summary data
summary_data <- data.frame(
    file_type = character(0),
    input_file = character(0),
    output_file = character(0),
    num_genes = numeric(0),
    num_samples = numeric(0),
    annotation_included = logical(0),
    annotation_columns = numeric(0),
    stringsAsFactors = FALSE
)

# Process each file type
for (file_info in file_types) {
    file_path <- file_info$path
    file_type <- file_info$type
    suffix <- file_info$suffix
    
    # Process gene-level files with annotation, transcript-level files without
    use_annotation <- grepl("^gene_", file_type)
    processed_data <- process_tximport_file(file_path, file_type, if(use_annotation) annotation_data else NULL)
    
    if (!is.null(processed_data)) {
        # Create output file
        output_file <- paste0(opt$output_prefix, ".merged.", suffix)
        write.table(processed_data, 
                   file = output_file, 
                   sep = "\t", 
                   quote = FALSE, 
                   row.names = FALSE)
        
        cat("Saved merged", file_type, "to:", output_file, "\n")
        
        # Add to summary
        summary_data <- rbind(summary_data, data.frame(
            file_type = file_type,
            input_file = basename(file_path),
            output_file = basename(output_file),
            num_genes = nrow(processed_data),
            num_samples = ncol(processed_data) - if(use_annotation && !is.null(annotation_data)) ncol(annotation_data) else 1,
            annotation_included = use_annotation && !is.null(annotation_data),
            annotation_columns = if(use_annotation && !is.null(annotation_data)) ncol(annotation_data) - 1 else 0,
            stringsAsFactors = FALSE
        ))
    }
}

# Save summary
summary_file <- paste0(opt$output_prefix, "_merge_summary.txt")
write.table(summary_data, 
           file = summary_file, 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

cat("\nTximport merge completed successfully!\n")
if (!is.null(annotation_data)) {
    cat("Gene-level files include annotation with", ncol(annotation_data) - 1, "annotation columns\n")
    cat("Annotation columns included:", paste(colnames(annotation_data)[-1], collapse = ", "), "\n")
    cat("Note: Excluded has_transcript, has_exon, has_intron, has_5utr, has_3utr columns. All other annotation and metadata columns are included\n")
}

cat("Summary saved to:", summary_file, "\n")
cat("Summary:\n")
print(summary_data)