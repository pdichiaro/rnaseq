# Improvement: Unified Output Structure for Genome Quantification (STAR and HISAT2)

## Overview
Extended genome quantification support to work consistently for both STAR and HISAT2 aligners, ensuring output structure matches RSEM quantification for consistency across the pipeline.

## Problem Description

### Issue: HISAT2 Genome Quantification Not Configured
The pipeline configuration only supported STAR with genome quantification. When using HISAT2 with genome quantification, the merged counts and DESeq2 normalization outputs were not properly configured, even though the workflow code supported it via `MERGE_GENOME_COUNTS_HISAT2`.

### Issue: Inconsistent Conditional Checks
- **Genome quantification**: Only checked `params.aligner == 'star'`
- **DeepTools normalization**: Only checked `params.aligner == 'star'`
- Both should support `hisat2` as well

## Solution

### Change 1: Unified Genome Quantification Configuration
**Extended the genome quantification block to support both STAR and HISAT2:**

**Before:**
```groovy
if (params.aligner == 'star' && params.quantification == 'genome') {
    process {
        withName: '.*:GENOME_COUNT' { ... }
        withName: '.*:MERGE_GENOME_COUNTS' { ... }
        // Missing: MERGE_GENOME_COUNTS_HISAT2
    }
}
```

**After:**
```groovy
if ((['star', 'hisat2'].contains(params.aligner)) && params.quantification == 'genome') {
    process {
        withName: '.*:GENOME_COUNT' { ... }
        withName: '.*:MERGE_GENOME_COUNTS' { ... }
        withName: '.*:MERGE_GENOME_COUNTS_HISAT2' {
            publishDir = [
                path: { "${params.outdir}/${params.aligner}/genome" },
                mode: params.publish_dir_mode,
                saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
            ]
        }
    }
}
```

### Change 2: Extended DeepTools Support
**Updated deepTools bigWig normalization to support both aligners:**

**Before:**
```groovy
if (params.aligner == 'star' && !params.skip_deeptools_norm) {
```

**After:**
```groovy
if ((['star', 'hisat2'].contains(params.aligner)) && !params.skip_deeptools_norm) {
```

## Expected Output Structure

### For STAR + Genome Quantification
```
results/
└── star/
    └── genome/
        ├── genome_exon_counts_merged.txt           ✅ Merged counts at root
        ├── genome_transcript_counts_merged.txt
        ├── genome_intron_counts_merged.txt
        ├── genome_5utr_counts_merged.txt
        ├── genome_3utr_counts_merged.txt
        ├── deseq2/
        │   ├── all_genes/
        │   │   ├── Quality_Control/
        │   │   │   ├── Sample_Distance_Heatmap.pdf
        │   │   │   ├── Sample_Correlation_Heatmap.pdf
        │   │   │   └── PCA_Plot_*.pdf
        │   │   └── Read_Distribution/
        │   └── invariant_genes/
        │       ├── Quality_Control/
        │       └── Read_Distribution/
        └── deeptools/                              ✅ When skip_deeptools_norm=FALSE
            ├── all_genes/
            │   └── *.norm.bw
            └── invariant_genes/
                └── *.norm.bw
```

### For HISAT2 + Genome Quantification
```
results/
└── hisat2/
    └── genome/
        ├── genome_exon_counts_merged.txt           ✅ Now properly configured!
        ├── genome_transcript_counts_merged.txt
        ├── genome_intron_counts_merged.txt
        ├── genome_5utr_counts_merged.txt
        ├── genome_3utr_counts_merged.txt
        ├── deseq2/
        │   ├── all_genes/
        │   │   ├── Quality_Control/
        │   │   └── Read_Distribution/
        │   └── invariant_genes/
        │       ├── Quality_Control/
        │       └── Read_Distribution/
        └── deeptools/                              ✅ Now supported!
            ├── all_genes/
            │   └── *.norm.bw
            └── invariant_genes/
                └── *.norm.bw
```

### Comparison with RSEM Structure
The structure now matches RSEM quantification for consistency:

**RSEM:**
```
star/rsem/
├── *.genes.results (merged counts)
├── deseq2/
│   ├── all_genes/
│   └── invariant_genes/
└── deeptools/
```

**Genome (STAR or HISAT2):**
```
star/genome/  or  hisat2/genome/
├── genome_*_counts_merged.txt (merged counts)
├── deseq2/
│   ├── all_genes/
│   └── invariant_genes/
└── deeptools/
```

## Files Modified

### Modified
- **Workflow Configuration:** `workflows/rnaseq/nextflow.config`
  - Extended genome quantification condition to include HISAT2
  - Added `MERGE_GENOME_COUNTS_HISAT2` configuration block
  - Extended deepTools condition to include HISAT2

## Key Benefits

### 1. Aligner Parity
- ✅ Both STAR and HISAT2 now have identical genome quantification support
- ✅ Configuration is aligner-agnostic where possible

### 2. Consistent Structure
- ✅ Genome quantification structure matches RSEM quantification structure
- ✅ Merged counts at `{aligner}/{quantification}/` root level
- ✅ DESeq2 outputs in `{aligner}/{quantification}/deseq2/{method}/`
- ✅ DeepTools bigWigs in `{aligner}/{quantification}/deeptools/{method}/`

### 3. Feature Completeness
- ✅ HISAT2 users can now use genome quantification with full DESeq2 QC and deepTools normalization
- ✅ All quantification methods (RSEM, genome) work consistently across aligners (STAR, HISAT2)

## Testing Recommendations

### Test 1: STAR + Genome Quantification
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner star \
  --quantification genome

# Verify structure
ls -lh results/star/genome/*.txt
ls -lh results/star/genome/deseq2/*/Quality_Control/*.pdf
ls -lh results/star/genome/deeptools/*/*.norm.bw
```

### Test 2: HISAT2 + Genome Quantification
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner hisat2 \
  --quantification genome

# Verify structure (should match STAR structure)
ls -lh results/hisat2/genome/*.txt
ls -lh results/hisat2/genome/deseq2/*/Quality_Control/*.pdf
ls -lh results/hisat2/genome/deeptools/*/*.norm.bw
```

### Test 3: Verify Merged Counts
```bash
# Verify merged counts are at root of genome/ directory
head results/star/genome/genome_exon_counts_merged.txt
head results/hisat2/genome/genome_exon_counts_merged.txt
```

## Technical Details

### Process Names in Workflow
- **STAR genome quantification** uses: `MERGE_GENOME_COUNTS`
- **HISAT2 genome quantification** uses: `MERGE_GENOME_COUNTS_HISAT2`
- Both are aliases of the same module, just with different names

### Configuration Pattern
The configuration uses Groovy's `.contains()` method with a list for cleaner multi-aligner support:

```groovy
if ((['star', 'hisat2'].contains(params.aligner)) && params.quantification == 'genome') {
    // Works for both aligners
}
```

This is more maintainable than:
```groovy
if ((params.aligner == 'star' || params.aligner == 'hisat2') && params.quantification == 'genome') {
```

## Related Changes

This improvement builds on previous fixes:
1. **Quality Control PDFs** - Already configured with `${params.aligner}` for dynamic paths
2. **DeepTools BigWig** - Already using `${params.quantification}` for dynamic paths
3. Now extended to support both STAR and HISAT2 aligners

## Impact

### No Breaking Changes
- ✅ STAR genome quantification behavior unchanged
- ✅ Existing outputs remain in same locations
- ✅ Purely additive - enables HISAT2 genome quantification

### Consistency Improvements
- ✅ Aligner-agnostic configuration where possible
- ✅ Unified structure across quantification methods
- ✅ Better code organization and maintainability

## Notes

### Why This Matters
The pipeline already had the workflow logic to handle HISAT2 genome quantification (via `MERGE_GENOME_COUNTS_HISAT2`), but the configuration wasn't set up to publish the outputs correctly. This fix ensures the configuration matches the workflow capabilities.

### Future-Proofing
Using the list-based `.contains()` pattern makes it easier to add support for additional aligners in the future without duplicating configuration blocks.
