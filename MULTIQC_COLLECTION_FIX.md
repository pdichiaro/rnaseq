# MultiQC Collection Fix - DESeq2 QC Files

## Problem Identified

DESeq2 QC TSV files were **not being collected by MultiQC** despite:
- ✅ Files being generated correctly by R scripts (without quotes)
- ✅ Files being published to results directories
- ✅ MultiQC config correctly defining custom content patterns
- ✅ Workflow mixing files into `ch_multiqc_files` channel

### Root Cause

The issue was **premature channel collection** in the workflow. The TSV outputs from DESeq2 normalization processes were being collected using `.collect()` immediately after being emitted, which caused issues with how Nextflow stages files for downstream processes.

**Original Pattern:**
```groovy
ch_multiqc_files = ch_multiqc_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.sample_distances_txt.collect().ifEmpty([]))
```

**Problem:** Using `.collect()` creates a list containing all files, but when combined with `.ifEmpty([])`, empty channels become `[[]]` (a list containing an empty list), which prevents proper file staging in MultiQC's work directory.

## Solution

**Remove all premature `.collect()` operators** from individual file outputs. Let files flow through channels as individual items, and rely on the single `.collect()` at the MultiQC process call to gather all files properly.

**Fixed Pattern:**
```groovy
ch_multiqc_files = ch_multiqc_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.sample_distances_txt)
```

## Changes Made

### Modified File
- `workflows/rnaseq/main.nf`

### Scope
Removed `.collect()` from **all** DESeq2 QC file outputs across **all** quantifiers:

1. **Pseudo-alignment** (Kallisto/Salmon)
   - `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO`
   - `NORMALIZE_DESEQ2_QC_ALL_GENES_PSEUDO`

2. **STAR + RSEM**
   - `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_ALIGNMENT` (RSEM)
   - `NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT` (RSEM)

3. **STAR + Salmon**
   - `NORMALIZE_DESEQ2_QC_INVARIANT_GENES` (Salmon)
   - `NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT` (Salmon)

4. **STAR + Genome counts**
   - `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_ALIGNMENT` (Genome)
   - `NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT` (Genome)

5. **HISAT2 + Genome counts**
   - `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_ALIGNMENT` (HISAT2)
   - `NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT` (HISAT2)

### File Types Fixed
For each quantifier and normalization method:
- `sample_distances_txt` (sample distance matrices)
- `pca_all_genes_txt` (PCA using all genes)
- `pca_top_genes_txt` (PCA using top 500 genes)
- `read_dist_norm_txt` (normalized read distribution)

**Total:** 40 `.collect()` operators removed

## How It Works Now

### Channel Flow
```
DESeq2 Process
    ↓ (emits individual files)
mix into ch_multiqc_files
    ↓ (files remain as individual items)
ch_multiqc_files.collect()  ← Single collection point
    ↓ (all files gathered into list)
MULTIQC process
    ↓ (files properly staged)
MultiQC scans and matches patterns
```

### Key Point
The **single `.collect()` at line 1176** in the MULTIQC process call:
```groovy
MULTIQC (
    ch_multiqc_files.collect(),  // ← This collects everything once
    ...
)
```

This ensures:
1. All files remain as individual channel items until MultiQC
2. Files are staged correctly with `stageAs: "?/*"` in MultiQC process
3. MultiQC can find and match files using regex patterns from config

## Expected Outcome

After this fix, MultiQC should now properly collect and display:

### For Kallisto (or any quantifier)
1. **Read Distribution (All Genes)** - Box plot
2. **Read Distribution (Invariant Genes)** - Box plot
3. **Sample Distance (All Genes)** - Heatmap
4. **Sample Distance (Invariant Genes)** - Heatmap
5. **PCA All Genes (All Genes Normalization)** - Scatter plot
6. **PCA All Genes (Invariant Genes Normalization)** - Scatter plot
7. **PCA Top 500 (All Genes Normalization)** - Scatter plot
8. **PCA Top 500 (Invariant Genes Normalization)** - Scatter plot

**Total: 8 interactive plots per quantifier**

## Testing

To verify the fix:

1. Run the pipeline with Kallisto:
```bash
nextflow run . -profile test,docker --pseudo_aligner kallisto
```

2. Check MultiQC report sections for "DESeq2 Kallisto QC"

3. Verify that all 8 plots appear in the report

4. Check MultiQC logs for successful file matching:
```bash
grep "deseq2-kallisto" results/multiqc/multiqc_data/multiqc.log
```

## Related Issues

This fix resolves the issue where DESeq2 QC sections were missing from MultiQC reports despite:
- Files being generated (commit 23d9114)
- Quotes removed from TSV files (commit a98ca65)
- MultiQC configuration correct

The only remaining issue was the channel collection pattern preventing proper file staging.

## Commit
- **Hash:** 315cc35
- **Date:** 2025-11-09
- **Files Changed:** 1 (workflows/rnaseq/main.nf)
- **Lines Changed:** 40 insertions(+), 40 deletions(-)
