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
# 1) Genomic features RDS file (from annotation_matrix.r)
# 2) BAM file path
# 3) Sample ID
# 4) Strand information
# 5) Paired-end information
# 6) Output path
# 7) Number of cores

option_list <- list(
    make_option(c("-f", "--features"), type="character", default=NULL, help="Genomic features RDS file path", metavar="character"),
    make_option(c("-b", "--bam"), type="character", default=NULL, help="BAM file path", metavar="character"),
    make_option(c("-i", "--id"), type="character", default=NULL, help="Sample ID", metavar="character"),
    make_option(c("-s", "--strand"), type="character", default=NULL, help="Sample strand (forward/reverse/non_specific)", metavar="character"),
    make_option(c("-p", "--paired_end"), type="character", default=NULL, help="single/paired", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character"),
    make_option(c("-c", "--cpus"), type="integer", default=1, help="Number of cores", metavar="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

cpus = as.numeric(opt$cpus)

# Load genomic features from RDS file
genomic_features <- readRDS(opt$features)
TxbyGene <- genomic_features$TxbyGene
ExonByGene <- genomic_features$ExonByGene
IntronByGene <- genomic_features$IntronByGene
UTR5ByTranscript <- genomic_features$UTR5ByTranscript
UTR3ByTranscript <- genomic_features$UTR3ByTranscript

# Set up BAM file reading
bamfile <- BamFile(opt$bam, yieldSize=1000000)
param <- NULL

# Define the strandedness and single/paired end parameters
if(opt$paired_end == "single"){
    singleEnd <- TRUE
    ignorestrand <- FALSE
}else{
    singleEnd <- FALSE
    ignorestrand <- FALSE
}

# Set strand-specific parameters
if(opt$strand == "reverse"){
    strandfun <- invertStrand
    ignorestrand <- FALSE
}else{
    if(opt$strand == "non_specific"){
        strandfun <- NULL
        ignorestrand <- TRUE
    }else{
        if(opt$strand == "forward"){
            strandfun <- NULL
            ignorestrand <- FALSE
        }
    }
}

cat("Starting read counting for sample:", opt$id, "\n")
cat("Strand:", opt$strand, "\n")
cat("Library type:", opt$paired_end, "\n")

# Count reads in transcripts (gene-level)
cat("Counting reads in transcripts...\n")
se_all <- summarizeOverlaps(features=TxbyGene,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Count reads in introns
cat("Counting reads in introns...\n")
se_int <- summarizeOverlaps(features=IntronByGene,
    reads=bamfile,
    mode="IntersectionNotEmpty",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Count reads in exons
cat("Counting reads in exons...\n")
se_exonic <- summarizeOverlaps(features=ExonByGene,
    reads=bamfile,
    mode="IntersectionStrict",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Count reads in 5' UTRs
cat("Counting reads in 5' UTRs...\n")
se_5_utr <- summarizeOverlaps(features=UTR5ByTranscript,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Count reads in 3' UTRs
cat("Counting reads in 3' UTRs...\n")
se_3_utr <- summarizeOverlaps(features=UTR3ByTranscript,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Extract counts from summarized experiments
dat_all <- as.data.frame(assay(se_all)) %>%
    rownames_to_column(var = "gene_id")

dat_int <- as.data.frame(assay(se_int)) %>%
    rownames_to_column(var = "gene_id")

dat_ex <- as.data.frame(assay(se_exonic)) %>%
    rownames_to_column(var = "gene_id")

dat_5_utr <- as.data.frame(assay(se_5_utr)) %>%
    rownames_to_column(var = "gene_id")

dat_3_utr <- as.data.frame(assay(se_3_utr)) %>%
    rownames_to_column(var = "gene_id")

# Rename count columns with sample ID
colnames(dat_all) <- c("gene_id", paste0(opt$id, "_transcript"))
colnames(dat_int) <- c("gene_id", paste0(opt$id, "_intron"))
colnames(dat_ex) <- c("gene_id", paste0(opt$id, "_exon"))
colnames(dat_5_utr) <- c("gene_id", paste0(opt$id, "_5utr"))
colnames(dat_3_utr) <- c("gene_id", paste0(opt$id, "_3utr"))

# Set rownames
rownames(dat_all) <- dat_all$gene_id
rownames(dat_int) <- dat_int$gene_id
rownames(dat_ex) <- dat_ex$gene_id
rownames(dat_5_utr) <- dat_5_utr$gene_id
rownames(dat_3_utr) <- dat_3_utr$gene_id

# Create output directory
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Save individual count matrices
write.table(dat_all, 
            file = paste0(opt$output, "_", opt$id, "_transcript_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(dat_int, 
            file = paste0(opt$output, "_", opt$id, "_intron_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(dat_ex, 
            file = paste0(opt$output, "_", opt$id, "_exon_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(dat_5_utr, 
            file = paste0(opt$output, "_", opt$id, "_5utr_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(dat_3_utr, 
            file = paste0(opt$output, "_", opt$id, "_3utr_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Combine all counts into a single matrix for this sample
combined_counts <- dat_all %>%
    full_join(dat_int, by = "gene_id") %>%
    full_join(dat_ex, by = "gene_id") %>%
    full_join(dat_5_utr, by = "gene_id") %>%
    full_join(dat_3_utr, by = "gene_id")

# Replace NA values with 0
combined_counts[is.na(combined_counts)] <- 0

# Save combined counts
write.table(combined_counts, 
            file = paste0(opt$output, "_", opt$id, "_combined_counts.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save as RDS for easy loading in R
res_list <- list(
    transcript = dat_all, 
    intron = dat_int, 
    exon = dat_ex, 
    utr5 = dat_5_utr, 
    utr3 = dat_3_utr,
    combined = combined_counts
)

saveRDS(res_list, file = paste0(opt$output, "_", opt$id, "_counts.rds"))

# Generate summary statistics
total_reads_transcript <- sum(dat_all[, 2])
total_reads_intron <- sum(dat_int[, 2])
total_reads_exon <- sum(dat_ex[, 2])
total_reads_5utr <- sum(dat_5_utr[, 2])
total_reads_3utr <- sum(dat_3_utr[, 2])

summary_stats <- data.frame(
    sample_id = opt$id,
    feature_type = c("transcript", "intron", "exon", "5utr", "3utr"),
    total_reads = c(total_reads_transcript, total_reads_intron, total_reads_exon, total_reads_5utr, total_reads_3utr),
    genes_with_reads = c(
        sum(dat_all[, 2] > 0),
        sum(dat_int[, 2] > 0),
        sum(dat_ex[, 2] > 0),
        sum(dat_5_utr[, 2] > 0),
        sum(dat_3_utr[, 2] > 0)
    ),
    stringsAsFactors = FALSE
)

write.table(summary_stats, 
            file = paste0(opt$output, "_", opt$id, "_summary.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Read counting completed successfully!\n")
cat("Output files:\n")
cat("- Combined counts:", paste0(opt$output, "_", opt$id, "_combined_counts.txt"), "\n")
cat("- Individual count files saved with prefixes\n")
cat("- RDS file:", paste0(opt$output, "_", opt$id, "_counts.rds"), "\n")
cat("- Summary:", paste0(opt$output, "_", opt$id, "_summary.txt"), "\n")

cat("\nSummary statistics:\n")
print(summary_stats)