#!/usr/bin/env Rscript

library("GenomicFeatures")
library("GenomicAlignments")
library(GenomicRanges)
library("Rsamtools")
library("BiocParallel")
library("openxlsx")
library("tidyverse")
library("yaml")
library(biomaRt)
library(rtracklayer)
library(optparse)

# Custom functions:
`%!in%` <- Negate(`%in%`)

# R function to parse parameters using optparse.
# This will take in: 
# 1) TxDB sqlite file path 
# 2) metadata file path
# 3) output file path
# 4) number of cores

option_list <- list(
    make_option(c("-q", "--txdb_sqlite"), type="character", default=NULL, help="sqlite txdb file path", metavar="character"),
    make_option(c("-m", "--meta"), type="character", default=NULL, help="meta file txdb", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character"),
    make_option(c("-c", "--cpus"), type="integer", default=1, help="Import n cores cpus", metavar="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

cpus = as.numeric(opt$cpus)

# Load the txdb_sqlite file
txdb <- loadDb(opt$txdb_sqlite)
Meta <- read.delim(opt$meta, header=TRUE, sep="\t")

# Extract genomic features from TxDB
TxbyGene <- transcriptsBy(txdb, by="gene")
ExonByGene <- exonsBy(txdb, by="gene")
IntronByGene <- intronsByTranscript(txdb, use.names=TRUE)
UTR5ByTranscript <- fiveUTRsByTranscript(txdb, use.names=TRUE)
UTR3ByTranscript <- threeUTRsByTranscript(txdb, use.names=TRUE)

# Create annotation matrices for different genomic features
# Get gene IDs from each feature type
gene_ids_tx <- names(TxbyGene)
gene_ids_exon <- names(ExonByGene)
gene_ids_intron <- names(IntronByGene)
gene_ids_utr5 <- names(UTR5ByTranscript)
gene_ids_utr3 <- names(UTR3ByTranscript)

# Create comprehensive gene annotation matrix
all_gene_ids <- unique(c(gene_ids_tx, gene_ids_exon, gene_ids_intron, gene_ids_utr5, gene_ids_utr3))

# Initialize annotation matrix
annotation_matrix <- data.frame(
    gene_id = all_gene_ids,
    has_transcript = all_gene_ids %in% gene_ids_tx,
    has_exon = all_gene_ids %in% gene_ids_exon,
    has_intron = all_gene_ids %in% gene_ids_intron,
    has_5utr = all_gene_ids %in% gene_ids_utr5,
    has_3utr = all_gene_ids %in% gene_ids_utr3,
    stringsAsFactors = FALSE
)

# Add transcript counts per gene
annotation_matrix$transcript_count <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(TxbyGene)) {
        length(TxbyGene[[x]])
    } else {
        0
    }
})

# Add exon counts per gene
annotation_matrix$exon_count <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(ExonByGene)) {
        length(ExonByGene[[x]])
    } else {
        0
    }
})

# Add intron counts per gene (sum across all transcripts)
annotation_matrix$intron_count <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(IntronByGene)) {
        length(IntronByGene[[x]])
    } else {
        0
    }
})

# Calculate total transcript length per gene
annotation_matrix$total_transcript_length <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(TxbyGene)) {
        tx_ranges <- TxbyGene[[x]]
        if (length(tx_ranges) > 0) {
            sum(width(GenomicRanges::reduce(tx_ranges)))
        } else {
            0
        }
    } else {
        0
    }
})

# Calculate total exon length per gene
annotation_matrix$total_exon_length <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(ExonByGene)) {
        exon_ranges <- ExonByGene[[x]]
        if (length(exon_ranges) > 0) {
            sum(width(GenomicRanges::reduce(exon_ranges)))
        } else {
            0
        }
    } else {
        0
    }
})

# Calculate total intron length per gene
annotation_matrix$total_intron_length <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(IntronByGene)) {
        intron_ranges <- IntronByGene[[x]]
        if (length(intron_ranges) > 0) {
            sum(width(intron_ranges))
        } else {
            0
        }
    } else {
        0
    }
})

# Add chromosome information
annotation_matrix$chromosome <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(TxbyGene)) {
        tx_ranges <- TxbyGene[[x]]
        if (length(tx_ranges) > 0) {
            unique(as.character(seqnames(tx_ranges)))[1]
        } else {
            NA
        }
    } else if (x %in% names(ExonByGene)) {
        exon_ranges <- ExonByGene[[x]]
        if (length(exon_ranges) > 0) {
            unique(as.character(seqnames(exon_ranges)))[1]
        } else {
            NA
        }
    } else {
        NA
    }
})

# Add strand information
annotation_matrix$strand <- sapply(annotation_matrix$gene_id, function(x) {
    if (x %in% names(TxbyGene)) {
        tx_ranges <- TxbyGene[[x]]
        if (length(tx_ranges) > 0) {
            unique(as.character(strand(tx_ranges)))[1]
        } else {
            NA
        }
    } else if (x %in% names(ExonByGene)) {
        exon_ranges <- ExonByGene[[x]]
        if (length(exon_ranges) > 0) {
            unique(as.character(strand(exon_ranges)))[1]
        } else {
            NA
        }
    } else {
        NA
    }
})

# Merge with metadata if available
if (nrow(Meta) > 0) {
    # Try to match gene IDs with metadata
    # Assuming metadata has gene information in one of its columns
    meta_gene_cols <- grep("gene", colnames(Meta), ignore.case = TRUE)
    if (length(meta_gene_cols) > 0) {
        # Use the first gene column found
        gene_col <- colnames(Meta)[meta_gene_cols[1]]
        annotation_matrix <- merge(annotation_matrix, Meta, 
                                 by.x = "gene_id", by.y = gene_col, 
                                 all.x = TRUE, all.y = FALSE)
    }
}

# Create output directory if it doesn't exist
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Save annotation matrix
write.table(annotation_matrix, 
            file = paste0(opt$output, "_annotation_matrix.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# Save genomic ranges objects for downstream analysis
genomic_features <- list(
    TxbyGene = TxbyGene,
    ExonByGene = ExonByGene,
    IntronByGene = IntronByGene,
    UTR5ByTranscript = UTR5ByTranscript,
    UTR3ByTranscript = UTR3ByTranscript
)

saveRDS(genomic_features, file = paste0(opt$output, "_genomic_features.rds"))

# Create summary statistics
summary_stats <- data.frame(
    feature_type = c("Genes", "Transcripts", "Exons", "Introns", "5'UTRs", "3'UTRs"),
    count = c(
        length(all_gene_ids),
        sum(annotation_matrix$transcript_count),
        sum(annotation_matrix$exon_count),
        sum(annotation_matrix$intron_count),
        sum(annotation_matrix$has_5utr),
        sum(annotation_matrix$has_3utr)
    ),
    stringsAsFactors = FALSE
)

write.table(summary_stats, 
            file = paste0(opt$output, "_summary_stats.txt"), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

cat("Annotation matrix created successfully!\n")
cat("Output files:\n")
cat("- Annotation matrix:", paste0(opt$output, "_annotation_matrix.txt"), "\n")
cat("- Genomic features:", paste0(opt$output, "_genomic_features.rds"), "\n")
cat("- Summary statistics:", paste0(opt$output, "_summary_stats.txt"), "\n")