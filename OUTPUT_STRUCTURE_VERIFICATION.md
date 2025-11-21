# Output Structure Verification: STAR/HISAT2 + Genome Quantification

## Overview
This document verifies the complete output structure after merge for genome quantification with STAR or HISAT2 aligners.

## Complete Output Structure

### STAR + Genome Quantification
```
results/
└── star/
    └── genome/
        ├── ──────────────────────────────────────
        │  PER-SAMPLE COUNT FILES (from GENOME_COUNT)
        ├── ──────────────────────────────────────
        ├── sample1_transcript_counts.txt
        ├── sample1_intron_counts.txt
        ├── sample1_exon_counts.txt
        ├── sample1_5utr_counts.txt
        ├── sample1_3utr_counts.txt
        ├── sample1_combined_counts.txt
        ├── sample1_summary.txt
        ├── sample2_transcript_counts.txt
        ├── sample2_intron_counts.txt
        ├── ... (all samples)
        │
        ├── ──────────────────────────────────────
        │  MERGED COUNT FILES (from MERGE_GENOME_COUNTS)
        ├── ──────────────────────────────────────
        ├── genome_transcript_counts_merged.txt    ✅ At root level
        ├── genome_intron_counts_merged.txt        ✅ At root level
        ├── genome_exon_counts_merged.txt          ✅ At root level
        ├── genome_5utr_counts_merged.txt          ✅ At root level
        ├── genome_3utr_counts_merged.txt          ✅ At root level
        ├── genome_transcript_merge_summary.txt
        ├── genome_intron_merge_summary.txt
        ├── genome_exon_merge_summary.txt
        ├── genome_5utr_merge_summary.txt
        ├── genome_3utr_merge_summary.txt
        │
        ├── ──────────────────────────────────────
        │  DESEQ2 NORMALIZATION & QC
        ├── ──────────────────────────────────────
        └── deseq2/
            ├── all_genes/
            │   ├── Quality_Control/
            │   │   ├── Sample_Distance_Heatmap.pdf
            │   │   ├── Sample_Correlation_Heatmap.pdf
            │   │   ├── PCA_Plot_All_Genes.pdf
            │   │   ├── PCA_Plot_Top500_Genes.pdf
            │   │   └── PCA_Decomposed_*.pdf
            │   ├── Read_Distribution/
            │   │   └── Read_Distribution.pdf
            │   ├── scaling_dat.txt
            │   ├── normalized_counts.txt
            │   ├── rlog_counts.txt
            │   ├── size_factors
            │   ├── *.pca.vals.txt
            │   ├── *.sample.dists.txt
            │   └── *_mqc.tsv files
            │
            └── invariant_genes/
                ├── Quality_Control/
                │   ├── Sample_Distance_Heatmap.pdf
                │   ├── Sample_Correlation_Heatmap.pdf
                │   ├── PCA_Plot_All_Genes.pdf
                │   ├── PCA_Plot_Top500_Genes.pdf
                │   └── PCA_Decomposed_*.pdf
                ├── Read_Distribution/
                │   └── Read_Distribution.pdf
                ├── normalization/
                ├── scaling_dat.txt
                ├── normalized_counts.txt
                ├── normalized_filt.txt
                └── *_mqc.tsv files
        
        ├── ──────────────────────────────────────
        │  DEEPTOOLS NORMALIZED BIGWIG (when skip_deeptools_norm=FALSE)
        ├── ──────────────────────────────────────
        └── deeptools/
            ├── all_genes/
            │   ├── sample1.norm.bw
            │   ├── sample2.norm.bw
            │   └── sample3.norm.bw
            └── invariant_genes/
                ├── sample1.norm.bw
                ├── sample2.norm.bw
                └── sample3.norm.bw
```

### HISAT2 + Genome Quantification
```
results/
└── hisat2/
    └── genome/
        ├── [Same structure as STAR above]
        ├── Per-sample count files
        ├── genome_*_counts_merged.txt             ✅ At root level
        ├── genome_*_merge_summary.txt
        ├── deseq2/
        │   ├── all_genes/
        │   └── invariant_genes/
        └── deeptools/
            ├── all_genes/
            └── invariant_genes/
```

## Comparison with RSEM Structure

### RSEM Quantification Structure
```
results/
└── star/
    └── rsem/
        ├── sample1.genes.results              ← Per-sample results at root
        ├── sample2.genes.results
        ├── sample1.isoforms.results
        ├── sample2.isoforms.results
        ├── merged_gene_counts.txt             ← Merged counts at root
        ├── merged_transcript_counts.txt
        ├── deseq2/
        │   ├── all_genes/
        │   └── invariant_genes/
        └── deeptools/
            ├── all_genes/
            └── invariant_genes/
```

### Genome Quantification Structure
```
results/
└── star/
    └── genome/
        ├── sample1_exon_counts.txt            ← Per-sample results at root
        ├── sample2_exon_counts.txt
        ├── sample1_transcript_counts.txt
        ├── sample2_transcript_counts.txt
        ├── genome_exon_counts_merged.txt      ← Merged counts at root ✅
        ├── genome_transcript_counts_merged.txt
        ├── deseq2/
        │   ├── all_genes/
        │   └── invariant_genes/
        └── deeptools/
            ├── all_genes/
            └── invariant_genes/
```

**✅ STRUCTURE IS CONSISTENT!** Both RSEM and Genome quantification follow the same pattern:
- Per-sample results at root level of `{aligner}/{quantification}/`
- Merged counts at root level of `{aligner}/{quantification}/`
- DESeq2 outputs in `{aligner}/{quantification}/deseq2/{method}/`
- DeepTools bigWigs in `{aligner}/{quantification}/deeptools/{method}/`

## File Details

### Merged Count Files (Primary Analysis Files)
These are the key output files that contain merged counts across all samples:

1. **`genome_exon_counts_merged.txt`**
   - Used as input for DESeq2 normalization
   - Contains exon-level counts for all genes across all samples
   - **Primary file for downstream analysis**

2. **`genome_transcript_counts_merged.txt`**
   - Transcript-level counts across all samples
   - Alternative to exon counts for isoform-level analysis

3. **`genome_intron_counts_merged.txt`**
   - Intron retention analysis
   - Useful for alternative splicing studies

4. **`genome_5utr_counts_merged.txt`** & **`genome_3utr_counts_merged.txt`**
   - UTR-specific counts
   - Useful for 3' RNA-seq or UTR-specific analyses

### Merge Summary Files
These files document the merge process:
- `genome_exon_merge_summary.txt`
- `genome_transcript_merge_summary.txt`
- etc.

They contain information about:
- Number of samples merged
- Feature counts per sample
- Any merge warnings or issues

## Configuration Verification

### GENOME_COUNT Process
**PublishDir Configuration:**
```groovy
withName: '.*:GENOME_COUNT' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> 
            if (filename.equals('versions.yml') || filename.endsWith('_counts.rds')) {
                return null  // Exclude versions and RDS files
            }
            return filename  // Publish everything else
        }
    ]
}
```

**Published files:**
- ✅ `{sample}_transcript_counts.txt`
- ✅ `{sample}_intron_counts.txt`
- ✅ `{sample}_exon_counts.txt`
- ✅ `{sample}_5utr_counts.txt`
- ✅ `{sample}_3utr_counts.txt`
- ✅ `{sample}_combined_counts.txt`
- ✅ `{sample}_summary.txt`
- ❌ `{sample}_counts.rds` (excluded)
- ❌ `versions.yml` (excluded)

### MERGE_GENOME_COUNTS Process
**PublishDir Configuration:**
```groovy
withName: '.*:MERGE_GENOME_COUNTS' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

**Published files:**
- ✅ `genome_transcript_counts_merged.txt`
- ✅ `genome_intron_counts_merged.txt`
- ✅ `genome_exon_counts_merged.txt`
- ✅ `genome_5utr_counts_merged.txt`
- ✅ `genome_3utr_counts_merged.txt`
- ✅ `genome_*_merge_summary.txt`
- ❌ `versions.yml` (excluded)

### MERGE_GENOME_COUNTS_HISAT2 Process
**PublishDir Configuration:**
```groovy
withName: '.*:MERGE_GENOME_COUNTS_HISAT2' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Same as `MERGE_GENOME_COUNTS` - publishes all merged files to `hisat2/genome/` root.

## Downstream Usage

### DESeq2 Normalization Input
The DESeq2 normalization processes use the exon-level merged counts:

```groovy
// From main.nf
NORMALIZE_DESEQ2_QC_ALL_GENES_GENOME (
    MERGE_GENOME_COUNTS.out.merged_counts
        .flatten()
        .filter { it.name.contains('genome_exon_counts_merged.txt') }
        .first(),
    "STAR_Genome"
)
```

**Input file:** `genome_exon_counts_merged.txt` (from the root of `star/genome/`)

### DeepTools Normalization Input
DeepTools uses the scaling factors generated by DESeq2 normalization, which are derived from the exon counts.

## Summary

### ✅ Verified: Structure is Correct
1. **Merged counts at root level**: `{aligner}/genome/genome_*_counts_merged.txt`
2. **DESeq2 outputs organized**: `{aligner}/genome/deseq2/{method}/`
3. **DeepTools outputs organized**: `{aligner}/genome/deeptools/{method}/`
4. **Consistent with RSEM**: Same organizational pattern

### ✅ Verified: Both Aligners Supported
1. **STAR + genome**: Fully configured and working
2. **HISAT2 + genome**: Now fully configured (as of latest commit)

### ✅ Verified: All Features Published
1. **Per-sample counts**: Published at root level
2. **Merged counts**: Published at root level
3. **Merge summaries**: Published at root level
4. **DESeq2 QC**: Published in subdirectories with PDFs
5. **DeepTools bigWigs**: Published in subdirectories for both normalization methods

## Notes

### Why Multiple Feature Types?
Genome quantification provides counts at multiple feature levels:
- **Exon**: Most commonly used for gene expression analysis
- **Transcript**: For isoform-level analysis
- **Intron**: For intron retention and splicing analysis
- **5' UTR / 3' UTR**: For UTR-specific analyses

This is different from RSEM which provides:
- **Genes**: Gene-level counts
- **Isoforms**: Transcript-level counts

### Primary Analysis File
For most RNA-seq analyses, **`genome_exon_counts_merged.txt`** is the primary file, equivalent to RSEM's merged gene counts. This is why DESeq2 normalization specifically uses the exon-level merged counts.
