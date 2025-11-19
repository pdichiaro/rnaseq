# RNA-seq Pipeline Output Folder Structure

This document describes the organized output folder structure for the RNA-seq pipeline, structured by aligner and quantification method.

## Complete Folder Structure

```
results/
└── star/                                    # Aligner folder
    ├── rsem/                               # RSEM quantification method
    │   ├── <sample_id>/                    # Individual sample results
    │   │   ├── *.genes.results
    │   │   ├── *.isoforms.results
    │   │   └── *.stat/
    │   ├── merged_gene_counts.txt          # Merged count matrices
    │   ├── merged_gene_tpm.txt
    │   ├── merged_transcript_counts.txt
    │   ├── merged_transcript_tpm.txt
    │   └── deseq2/                         # DESeq2 analysis
    │       ├── all_genes/
    │       │   ├── size_factors/
    │       │   ├── R_sessionInfo.log
    │       │   └── *.RData
    │       └── invariant_genes/
    │           ├── size_factors/
    │           ├── R_sessionInfo.log
    │           └── *.RData
    │
    ├── genome/                             # Genome quantification method
    │   ├── <sample_id>.markdup.sorted.bam  # Individual BAM files
    │   ├── <sample_id>.markdup.sorted.bam.bai
    │   ├── merged_gene_counts.txt          # Merged count matrices
    │   └── deseq2/                         # DESeq2 analysis
    │       ├── all_genes/
    │       │   ├── size_factors/
    │       │   ├── R_sessionInfo.log
    │       │   └── *.RData
    │       └── invariant_genes/
    │           ├── size_factors/
    │           ├── R_sessionInfo.log
    │           └── *.RData
    │
    └── bigwig/                             # BigWig files (shared across quantification methods)
        ├── <sample_id>.forward.bigWig
        └── <sample_id>.reverse.bigWig
```

## Configuration Files Modified

### 1. `subworkflows/local/quantify_rsem/nextflow.config`
- **RSEM_CALCULATEEXPRESSION**: `${params.outdir}/${params.aligner}/rsem/${meta.id}`
- **RSEM_MERGE_COUNTS**: `${params.outdir}/${params.aligner}/rsem`

### 2. `workflows/rnaseq/nextflow.config`
- **DESEQ2_QC_RSEM** (all_genes): `${params.outdir}/${params.aligner}/rsem/deseq2/all_genes`
- **DESEQ2_QC_RSEM** (invariant_genes): `${params.outdir}/${params.aligner}/rsem/deseq2/invariant_genes`
- **GENOME_COUNT**: `${params.outdir}/${params.aligner}/genome`
- **MERGE_GENOME_COUNTS**: `${params.outdir}/${params.aligner}/genome`
- **DESEQ2_QC_STAR_GENOME** (all_genes): `${params.outdir}/${params.aligner}/genome/deseq2/all_genes`
- **DESEQ2_QC_STAR_GENOME** (invariant_genes): `${params.outdir}/${params.aligner}/genome/deseq2/invariant_genes`
- **NORMALIZE_BIGWIG** (all_genes): `${params.outdir}/${params.aligner}/${params.quantification}/bigwig_norm/all_genes`
- **NORMALIZE_BIGWIG** (invariant_genes): `${params.outdir}/${params.aligner}/${params.quantification}/bigwig_norm/invariant_genes`

## Key Design Principles

1. **Quantification Method Organization**: Each quantification method (rsem, genome) has its own subfolder under the aligner folder
2. **Consistent DESeq2 Structure**: All DESeq2 outputs follow the pattern: `<aligner>/<quantification>/deseq2/{all_genes,invariant_genes}`
3. **Normalized BigWig Files**: Use dynamic paths based on quantification method: `<aligner>/<quantification>/bigwig_norm/`
4. **Sample-Specific Outputs**: Individual sample results are organized under their respective quantification method folders

## Benefits

- **Clear Separation**: Each quantification method has dedicated folder space
- **Scalability**: Easy to add new quantification methods by creating new subfolders
- **Consistency**: All analysis outputs (DESeq2, normalized BigWig) follow the same hierarchical structure
- **Traceability**: Output paths clearly indicate aligner, quantification method, and analysis type
