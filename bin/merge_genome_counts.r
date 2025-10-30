#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")

# R function to merge genome-based counts from multiple samples
option_list <- list(
    make_option(c("-i", "--input_dir"), type="character", default=NULL, help="Input directory containing count files", metavar="character"),
    make_option(c("-p", "--pattern"), type="character", default="_combined_counts.txt", help="File pattern to match", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file prefix", metavar="character"),
    make_option(c("-t", "--feature_types"), type="character", default="transcript,intron,exon,5utr,3utr", help="Comma-separated list of feature types to merge", metavar="character"),
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
    
    # Define standard annotation matrix columns (for reference)
    standard_annotation_cols <- c("gene_id", "has_transcript", "has_exon", "has_intron", "has_5utr", "has_3utr",
                                 "transcript_count", "exon_count", "intron_count", 
                                 "total_transcript_length", "total_exon_length", "total_intron_length",
                                 "chromosome", "strand")
    
    # Get all available columns and exclude the unwanted ones
    available_cols <- colnames(full_annotation)
    selected_cols <- setdiff(available_cols, excluded_cols)
    
    annotation_data <- full_annotation[, selected_cols, drop = FALSE]
    
    cat("Selected annotation columns:", paste(colnames(annotation_data), collapse = ", "), "\n")
    cat("Filtered annotation to", ncol(annotation_data), "columns (", ncol(annotation_data) - 1, "annotation columns)\n")
}

# Get list of count files
count_files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)

if (length(count_files) == 0) {
    stop("No count files found matching pattern: ", opt$pattern)
}

cat("Found", length(count_files), "count files\n")

# Parse feature types
feature_types <- strsplit(opt$feature_types, ",")[[1]]

# Initialize merged data frames for each feature type
merged_counts <- list()

for (feature_type in feature_types) {
    merged_counts[[feature_type]] <- data.frame(gene_id = character(0))
}

# Process each count file
for (count_file in count_files) {
    cat("Processing:", basename(count_file), "\n")
    
    # Read count data
    count_data <- read.table(count_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Extract sample ID from filename
    sample_id <- gsub(paste0(opt$pattern, "$"), "", basename(count_file))
    sample_id <- gsub("^.*_", "", sample_id)  # Remove prefix if present
    
    # Process each feature type
    for (feature_type in feature_types) {
        # Find columns for this feature type
        feature_cols <- grep(paste0("_", feature_type, "$"), colnames(count_data), value = TRUE)
        
        if (length(feature_cols) > 0) {
            # Create a data frame for this sample and feature type
            sample_data <- count_data[, c("gene_id", feature_cols), drop = FALSE]
            
            # Rename the count column to sample ID
            colnames(sample_data)[2] <- sample_id
            
            # Merge with existing data
            if (nrow(merged_counts[[feature_type]]) == 0) {
                merged_counts[[feature_type]] <- sample_data
            } else {
                merged_counts[[feature_type]] <- merge(merged_counts[[feature_type]], sample_data, 
                                                     by = "gene_id", all = TRUE)
            }
        }
    }
}

# Replace NA values with 0
for (feature_type in feature_types) {
    merged_counts[[feature_type]][is.na(merged_counts[[feature_type]])] <- 0
}

# Create output directory
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Save merged count matrices with annotation information
for (feature_type in feature_types) {
    if (nrow(merged_counts[[feature_type]]) > 0) {
        # Merge with annotation data if available
        final_output <- merged_counts[[feature_type]]
        
        if (!is.null(annotation_data)) {
            # Merge count data with annotation information
            final_output <- merge(annotation_data, merged_counts[[feature_type]], 
                                by = "gene_id", all.y = TRUE, all.x = FALSE)
            cat("Merged", feature_type, "counts with annotation data\n")
        }
        
        output_file <- paste0(opt$output, "_", feature_type, "_counts_merged.txt")
        write.table(final_output, 
                   file = output_file, 
                   sep = "\t", 
                   quote = FALSE, 
                   row.names = FALSE)
        cat("Saved merged", feature_type, "counts to:", output_file, "\n")
    }
}

# Create a summary of merged data
summary_data <- data.frame(
    feature_type = feature_types,
    num_genes = sapply(feature_types, function(ft) nrow(merged_counts[[ft]])),
    num_samples = sapply(feature_types, function(ft) ncol(merged_counts[[ft]]) - 1),
    annotation_included = !is.null(annotation_data),
    annotation_columns = if (!is.null(annotation_data)) ncol(annotation_data) - 1 else 0,
    stringsAsFactors = FALSE
)

write.table(summary_data, 
           file = paste0(opt$output, "_merge_summary.txt"), 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE)

cat("Merge completed successfully!\n")
if (!is.null(annotation_data)) {
    cat("Annotation matrix integrated with", ncol(annotation_data) - 1, "annotation columns\n")
    cat("Annotation columns included:", paste(colnames(annotation_data)[-1], collapse = ", "), "\n")
    cat("Note: Excluded has_transcript, has_exon, has_intron, has_5utr, has_3utr columns. All other annotation and metadata columns are included\n")
}
cat("Summary:\n")
print(summary_data)