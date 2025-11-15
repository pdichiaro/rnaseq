#!/usr/bin/env Rscript

library("tidyverse")
library("optparse")

# R function to merge tximport counts from multiple samples with annotation integration
option_list <- list(
    make_option(c("--gene_counts"), type="character", default=NULL, help="Gene counts file from tximport", metavar="character"),
    make_option(c("--gene_tpm"), type="character", default=NULL, help="Gene TPM file from tximport", metavar="character"),
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
    
    # Check for gene name columns in annotation
    gene_name_cols_found <- grep("gene.*name", colnames(annotation_data), value = TRUE, ignore.case = TRUE)
    if (length(gene_name_cols_found) > 0) {
        cat("Gene name columns found in annotation:", paste(gene_name_cols_found, collapse = ", "), "\n")
    } else {
        cat("WARNING: No gene name columns found in annotation data\n")
    }
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
    cat("  Count data columns:", paste(head(colnames(data), 10), collapse = ", "))
    if (ncol(data) > 10) cat(" ... +", ncol(data) - 10, "more")
    cat("\n")
    
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
    
    # Check for gene name columns in count data
    count_gene_name_cols <- grep("gene.*name", colnames(data), value = TRUE, ignore.case = TRUE)
    if (length(count_gene_name_cols) > 0) {
        cat("  Gene name columns found in count data:", paste(count_gene_name_cols, collapse = ", "), "\n")
    } else {
        cat("  No gene name columns found in count data\n")
    }
    
    # Merge with annotation data if available
    if (!is.null(annotation_data)) {
        cat("Merging", file_type, "with annotation data...\n")
        
        # Check for duplicate gene_ids in both datasets
        data_dupes <- sum(duplicated(data$gene_id))
        if (data_dupes > 0) {
            cat("Warning: Found", data_dupes, "duplicate gene_ids in", file_type, "data\n")
            cat("Removing duplicates, keeping first occurrence\n")
            data <- data[!duplicated(data$gene_id), ]
        }
        
        annot_dupes <- sum(duplicated(annotation_data$gene_id))
        if (annot_dupes > 0) {
            cat("Warning: Found", annot_dupes, "duplicate gene_ids in annotation data\n")
            cat("Removing annotation duplicates, keeping first occurrence\n")
            annotation_data <- annotation_data[!duplicated(annotation_data$gene_id), ]
        }
        
        # Note: Keeping ENSG* patterns in gene_names (no cleaning applied)
        
        # Identify sample columns (numeric columns that are not in annotation)
        sample_cols <- setdiff(colnames(data), "gene_id")
        annotation_cols <- setdiff(colnames(annotation_data), "gene_id")
        
        cat("Sample columns in data (", length(sample_cols), "):", paste(head(sample_cols, 3), collapse = ", "))
        if (length(sample_cols) > 3) cat(" ...")
        cat("\n")
        cat("Annotation columns (", length(annotation_cols), "):", paste(head(annotation_cols, 3), collapse = ", "))
        if (length(annotation_cols) > 3) cat(" ...")
        cat("\n")
        
        # Debug: Show gene_name content in both datasets
        if ("gene_name" %in% colnames(annotation_data)) {
            annot_gene_examples <- annotation_data$gene_name[1:min(3, nrow(annotation_data))]
            cat("Annotation gene_name examples:", paste(annot_gene_examples, collapse = ", "), "\n")
        } else {
            cat("No gene_name column in annotation data\n")
        }
        
        if ("gene_name" %in% colnames(data)) {
            count_gene_examples <- data$gene_name[1:min(3, nrow(data))]
            cat("Count data gene_name examples:", paste(count_gene_examples, collapse = ", "), "\n")
        } else {
            cat("No gene_name column in count data\n")
        }
        
        # Perform the merge
        merged_data <- merge(annotation_data, data, by = "gene_id", all.y = TRUE, all.x = FALSE)
        cat("Columns after merge:", paste(colnames(merged_data), collapse = ", "), "\n")
        
        # Check for .x/.y columns immediately after merge
        xy_cols_immediate <- grep("\\.[xy]$", colnames(merged_data), value = TRUE)
        if (length(xy_cols_immediate) > 0) {
            cat("WARNING: .x/.y columns created by merge:", paste(xy_cols_immediate, collapse = ", "), "\n")
        } else {
            cat("INFO: No .x/.y columns created - merge had no column conflicts\n")
        }
        
        # Create gene_ids column (remove version suffix from gene_id)  
        # This handles cases like ENSG00000123456.1 -> ENSG00000123456
        merged_data$gene_ids <- gsub("\\.[0-9]+$", "", merged_data$gene_id)
        cat("Created gene_ids column (unversioned):", sum(merged_data$gene_ids != merged_data$gene_id), "of", nrow(merged_data), "genes had version suffixes removed\n")
        
        # Check for any remaining issues
        final_dupes <- sum(duplicated(merged_data$gene_id))
        if (final_dupes > 0) {
            cat("Warning: Still", final_dupes, "duplicate gene_ids after merge. Removing duplicates.\n")
            merged_data <- merged_data[!duplicated(merged_data$gene_id), ]
        }
        
        # Handle ALL .x/.y columns by consolidating them BEFORE column reordering
        # First find all base names that have .x/.y pairs
        all_xy_cols <- grep("\\.[xy]$", colnames(merged_data), value = TRUE)
        xy_base_names <- unique(gsub("\\.[xy]$", "", all_xy_cols))
        
        if (length(xy_base_names) > 0) {
            cat("CONSOLIDATION: Found .x/.y column pairs for:", paste(xy_base_names, collapse = ", "), "\n")
        } else {
            cat("INFO: No .x/.y column conflicts found - merge was clean\n")
        }
        
        # Track source choices globally for consistent consolidation
        source_choice <- character(nrow(merged_data))
        
        # Process gene_name-like columns first to establish source preference
        gene_name_patterns <- c("gene_name", "gene_names", "gene_symbol")
        gene_base_found <- NULL
        
        for (base_name in gene_name_patterns) {
            if (base_name %in% xy_base_names) {
                gene_base_found <- base_name
                break
            }
        }
        
        # If we found a gene name column, use it to establish source preference
        if (!is.null(gene_base_found)) {
            gene_name_x <- paste0(gene_base_found, ".x")
            gene_name_y <- paste0(gene_base_found, ".y")
            
            cat("CONSOLIDATION: Consolidating", gene_base_found, ".x and", gene_base_found, ".y columns to establish source preference\n")
            
            x_vals <- merged_data[[gene_name_x]]
            y_vals <- merged_data[[gene_name_y]]
            
            # Create consolidated gene_name: prefer non-NA values, prioritize those with ENSG or gene_name patterns
            consolidated <- character(nrow(merged_data))
            
            for (i in 1:nrow(merged_data)) {
                    x_val <- if(is.na(x_vals[i]) || x_vals[i] == "") NA else x_vals[i]
                    y_val <- if(is.na(y_vals[i]) || y_vals[i] == "") NA else y_vals[i]
                    
                    if (!is.na(x_val) && !is.na(y_val)) {
                        # Both have values - prefer one with ENSG or gene name pattern
                        x_has_ensg <- grepl("ENSG[0-9]{11}", x_val)
                        y_has_ensg <- grepl("ENSG[0-9]{11}", y_val)
                        x_has_gene <- grepl("[A-Z][A-Z0-9]+", x_val) && !grepl("^ENSG", x_val)
                        y_has_gene <- grepl("[A-Z][A-Z0-9]+", y_val) && !grepl("^ENSG", y_val)
                        
                        if ((x_has_ensg || x_has_gene) && !(y_has_ensg || y_has_gene)) {
                            consolidated[i] <- x_val
                            source_choice[i] <- "x"
                        } else if ((y_has_ensg || y_has_gene) && !(x_has_ensg || x_has_gene)) {
                            consolidated[i] <- y_val
                            source_choice[i] <- "y"
                        } else {
                            # If both or neither have patterns, prefer x (annotation usually)
                            consolidated[i] <- x_val
                            source_choice[i] <- "x"
                        }
                    } else if (!is.na(x_val)) {
                        consolidated[i] <- x_val
                        source_choice[i] <- "x"
                    } else if (!is.na(y_val)) {
                        consolidated[i] <- y_val
                        source_choice[i] <- "y"
                    } else {
                        consolidated[i] <- NA
                        source_choice[i] <- "none"
                    }
                }
                
                # Replace with consolidated column (use base_name, not hardcoded gene_name)
                merged_data[[base_name]] <- consolidated
                merged_data[[gene_name_x]] <- NULL
                merged_data[[gene_name_y]] <- NULL
                
                # Count non-NA values
                non_na_count <- sum(!is.na(consolidated) & consolidated != "")
                cat("SUCCESS: Consolidated", base_name, "column created:", non_na_count, "genes with names\n")
                
                # Show some examples
                example_indices <- which(!is.na(consolidated) & consolidated != "")[1:min(3, non_na_count)]
                for (i in example_indices) {
                    x_example <- if(i <= length(x_vals)) x_vals[i] else NA
                    y_example <- if(i <= length(y_vals)) y_vals[i] else NA
                    cat("  Example: .x='", x_example, "' .y='", y_example, "' -> '", consolidated[i], "'\n", sep = "")
                }
                
                # Now consolidate strand columns based on the same source choice as gene names
                strand_x <- "strand.x"
                strand_y <- "strand.y"
                
                if (strand_x %in% colnames(merged_data) && strand_y %in% colnames(merged_data)) {
                    cat("CONSOLIDATION: Consolidating strand.x and strand.y based on", base_name, "source choice\n")
                    
                    strand_x_vals <- merged_data[[strand_x]]
                    strand_y_vals <- merged_data[[strand_y]]
                    strand_consolidated <- character(nrow(merged_data))
                    
                    for (i in 1:nrow(merged_data)) {
                        if (source_choice[i] == "x" && !is.na(strand_x_vals[i]) && strand_x_vals[i] != "") {
                            strand_consolidated[i] <- strand_x_vals[i]
                        } else if (source_choice[i] == "y" && !is.na(strand_y_vals[i]) && strand_y_vals[i] != "") {
                            strand_consolidated[i] <- strand_y_vals[i]
                        } else if (source_choice[i] == "x") {
                            strand_consolidated[i] <- strand_x_vals[i]  # Keep even if NA/empty for consistency
                        } else if (source_choice[i] == "y") {
                            strand_consolidated[i] <- strand_y_vals[i]  # Keep even if NA/empty for consistency
                        } else {
                            # If no gene name was chosen, prefer x
                            strand_consolidated[i] <- strand_x_vals[i]
                        }
                    }
                    
                    # Replace with consolidated strand column
                    merged_data$strand <- strand_consolidated
                    merged_data[[strand_x]] <- NULL
                    merged_data[[strand_y]] <- NULL
                    
                    # Count non-NA strand values
                    strand_non_na_count <- sum(!is.na(strand_consolidated) & strand_consolidated != "")
                    cat("SUCCESS: Consolidated strand column created:", strand_non_na_count, "genes with strand info\n")
                    
                    # Show some examples
                    strand_example_indices <- which(!is.na(strand_consolidated) & strand_consolidated != "")[1:min(3, strand_non_na_count)]
                    for (i in strand_example_indices) {
                        source_used <- source_choice[i]
                        strand_x_example <- if(i <= length(strand_x_vals)) strand_x_vals[i] else NA
                        strand_y_example <- if(i <= length(strand_y_vals)) strand_y_vals[i] else NA
                        cat("  Example: .x='", strand_x_example, "' .y='", strand_y_example, "' -> '", strand_consolidated[i], "' (used .", source_used, " to match gene name)\n", sep = "")
                    }
                }
        }
        
        # Check for any remaining .x/.y column names and remove them
        remaining_xy_cols <- grep("\\.[xy]$", colnames(merged_data), value = TRUE)
        if (length(remaining_xy_cols) > 0) {
            cat("INFO: Other .x/.y columns found (not gene_name or strand):", paste(remaining_xy_cols, collapse = ", "), "\n")
            cat("REMOVING: Dropping remaining .x/.y duplicate columns...\n")
            
            # For each .x/.y pair, keep only the .x version (or .y if .x doesn't exist)
            base_names <- unique(gsub("\\.[xy]$", "", remaining_xy_cols))
            for (base_name in base_names) {
                x_col <- paste0(base_name, ".x")
                y_col <- paste0(base_name, ".y")
                
                if (x_col %in% colnames(merged_data) && y_col %in% colnames(merged_data)) {
                    # Both exist - consolidate by preferring non-NA values
                    x_vals <- merged_data[[x_col]]
                    y_vals <- merged_data[[y_col]]
                    consolidated <- ifelse(!is.na(x_vals) & x_vals != "", x_vals, y_vals)
                    merged_data[[base_name]] <- consolidated
                    merged_data[[x_col]] <- NULL
                    merged_data[[y_col]] <- NULL
                    cat("  Consolidated:", x_col, "+", y_col, "->", base_name, "\n")
                } else if (x_col %in% colnames(merged_data)) {
                    # Only .x exists
                    merged_data[[base_name]] <- merged_data[[x_col]]
                    merged_data[[x_col]] <- NULL
                    cat("  Renamed:", x_col, "->", base_name, "\n")
                } else if (y_col %in% colnames(merged_data)) {
                    # Only .y exists
                    merged_data[[base_name]] <- merged_data[[y_col]]
                    merged_data[[y_col]] <- NULL
                    cat("  Renamed:", y_col, "->", base_name, "\n")
                }
            }
            cat("SUCCESS: All .x/.y columns consolidated\n")
        }
        
        # Check and ensure gene_name column exists
        cat("CHECKING: gene_name column existence...\n")
        cat("  Current columns after consolidation:", paste(colnames(merged_data), collapse = ", "), "\n")
        
        if (!"gene_name" %in% colnames(merged_data)) {
            cat("WARNING: No gene_name column found after merge - creating fallback using gene_id\n")
            merged_data$gene_name <- merged_data$gene_id
            cat("SUCCESS: Created fallback gene_name column as copy of gene_id\n")
        } else {
            # Check if gene_name column has actual gene names or just gene_ids
            sample_gene_names <- merged_data$gene_name[1:min(5, nrow(merged_data))]
            gene_ids_only <- all(grepl("^ENSG[0-9]{11}", sample_gene_names, ignore.case = TRUE), na.rm = TRUE)
            
            if (gene_ids_only) {
                cat("WARNING: gene_name column contains only Ensembl IDs - may need better source data\n")
                cat("  Examples:", paste(head(sample_gene_names, 3), collapse = ", "), "\n")
            } else {
                cat("SUCCESS: gene_name column contains actual gene names\n")
                cat("  Examples:", paste(head(sample_gene_names, 3), collapse = ", "), "\n")
            }
        }
        
        # Final column reordering after all consolidation is complete
        # Specific order: gene_id, gene_ids, gene_name, chromosome, strand, other annotations, samples
        core_cols <- c("gene_id", "gene_ids", "gene_name")
        
        # Order annotation columns with chromosome before strand
        annotation_order <- c("chromosome", "strand", "gene_biotype", "gene_symbol", "description")
        ordered_annotation_cols <- c()
        
        # Add annotation columns in preferred order
        for (col in annotation_order) {
            if (col %in% colnames(merged_data)) {
                ordered_annotation_cols <- c(ordered_annotation_cols, col)
            }
        }
        
        # Add any remaining annotation columns not in the preferred order
        final_annotation_cols <- intersect(colnames(merged_data), annotation_cols)
        remaining_annotation_cols <- setdiff(final_annotation_cols, c(core_cols, ordered_annotation_cols))
        ordered_annotation_cols <- c(ordered_annotation_cols, remaining_annotation_cols)
        
        # Sample columns are everything else
        all_sample_cols <- setdiff(colnames(merged_data), c(core_cols, ordered_annotation_cols))
        
        final_ordered_cols <- c(core_cols, ordered_annotation_cols, all_sample_cols)
        merged_data <- merged_data[, final_ordered_cols, drop = FALSE]
        
        cat("Successfully merged", file_type, "with annotation:", nrow(merged_data), "genes\n")
        cat("Final column structure: gene_id, gene_ids, gene_name +", length(ordered_annotation_cols), "annotations +", length(all_sample_cols), "samples\n")
        cat("Annotation column order:", paste(ordered_annotation_cols, collapse = ", "), "\n")
        return(merged_data)
    } else {
        # For transcript-level data, still check for duplicates
        data_dupes <- sum(duplicated(data$gene_id))
        if (data_dupes > 0) {
            cat("Warning: Found", data_dupes, "duplicate gene_ids in", file_type, "data\n")
            cat("Removing duplicates, keeping first occurrence\n")
            data <- data[!duplicated(data$gene_id), ]
        }
        
        # Create gene_ids column (remove version suffix from gene_id) for non-annotation files too
        data$gene_ids <- gsub("\\.[0-9]+$", "", data$gene_id)
        cat("Created gene_ids column (unversioned):", sum(data$gene_ids != data$gene_id), "of", nrow(data), "genes had version suffixes removed\n")
        
        # Ensure gene_name column exists for non-annotation files too
        if (!"gene_name" %in% colnames(data)) {
            cat("WARNING: No gene_name column found in non-annotation file - creating fallback using gene_id\n")
            data$gene_name <- data$gene_id
            cat("SUCCESS: Created gene_name column as copy of gene_id\n")
        }
        
        # Reorder columns: gene_id, gene_ids, gene_name, then sample columns
        sample_cols <- setdiff(colnames(data), c("gene_id", "gene_ids", "gene_name"))
        ordered_cols <- c("gene_id", "gene_ids", "gene_name", sample_cols)
        data <- data[, ordered_cols, drop = FALSE]
        
        return(data)
    }
}

# Process each tximport file type
file_types <- list(
    list(path = opt$gene_counts, type = "gene_counts", suffix = "gene_counts.tsv"),
    list(path = opt$gene_tpm, type = "gene_tpm", suffix = "gene_tpm.tsv"),
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
cat("SUCCESS: Added 'gene_ids' column with unversioned gene IDs (e.g., ENSG00000123456 from ENSG00000123456.1)\n")
cat("SUCCESS: Consolidated .x/.y columns: gene_name and strand columns unified based on intelligent selection\n")
if (!is.null(annotation_data)) {
    cat("Gene-level files include annotation with", ncol(annotation_data) - 1, "annotation columns\n")
    cat("Annotation columns included:", paste(colnames(annotation_data)[-1], collapse = ", "), "\n")
    cat("Note: Excluded has_transcript, has_exon, has_intron, has_5utr, has_3utr columns. All other annotation and metadata columns are included\n")
}

# Check all output files for .x/.y column issues and column structure
cat("\n=== DIAGNOSTIC CHECK FOR COLUMN STRUCTURE ===\n")
output_files <- list.files(pattern = paste0(opt$output_prefix, "\\.merged\\..*\\.tsv"))
for (file in output_files) {
    if (file.exists(file)) {
        col_names <- colnames(read.table(file, header = TRUE, nrows = 1, sep = "\t"))
        xy_cols <- grep("\\.[xy]$", col_names, value = TRUE)
        
        cat("\nFILE:", file, "\n")
        cat("   Total columns:", length(col_names), "\n")
        cat("   First 5 columns:", paste(head(col_names, 5), collapse = ", "), "\n")
        cat("   Last 5 columns:", paste(tail(col_names, 5), collapse = ", "), "\n")
        
        if (length(xy_cols) > 0) {
            cat("   WARNING: .x/.y suffixes found:", paste(xy_cols, collapse = ", "), "\n")
        } else {
            cat("   SUCCESS: No column name conflicts detected\n")
        }
        
        # Check for gene_ids column
        if ("gene_ids" %in% col_names) {
            cat("   INFO: gene_ids column detected (unversioned gene IDs)\n")
        }
        
        # Try to identify annotation vs sample columns
        likely_annotation <- grep("^(gene_|transcript_|exon_|intron_|chromosome|strand|symbol|name|biotype)", col_names, value = TRUE)
        if (length(likely_annotation) > 0) {
            cat("   ANNOTATION columns (", length(likely_annotation), "):", paste(head(likely_annotation, 3), collapse = ", "))
            if (length(likely_annotation) > 3) cat(" ...")
            cat("\n")
        }
        
        # Remaining columns are likely samples
        likely_samples <- setdiff(col_names, c("gene_id", "gene_ids", likely_annotation))
        if (length(likely_samples) > 0) {
            cat("   SAMPLE columns (", length(likely_samples), "):", paste(head(likely_samples, 3), collapse = ", "))
            if (length(likely_samples) > 3) cat(" ...")
            cat("\n")
        }
    }
}

cat("Summary saved to:", summary_file, "\n")
cat("Summary:\n")
print(summary_data)