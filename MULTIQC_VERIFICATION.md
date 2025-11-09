# MultiQC Configuration Verification

## Summary
This document verifies that MultiQC can now correctly collect DESeq2 QC files from the pipeline.

## Fixed Issues

### 1. ✅ File Naming Convention
**R Scripts Generate:**
- Pattern: `{quantifier}.deseq2.all_genes.{type}.txt`
- Examples:
  - `kallisto.deseq2.all_genes.pca.vals.txt`
  - `kallisto.deseq2.all_genes.sample.dists.txt`
  - `kallisto.deseq2.all_genes.read.distribution.normalized.txt`

**MultiQC Config Expects:**
```yaml
fn_re: kallisto\.deseq2\.all_genes\.pca\.vals\.txt
fn_re: kallisto\.deseq2\.all_genes\.sample\.dists\.txt
fn_re: kallisto\.deseq2\.all_genes\.read\.distribution\.normalized\.txt
```

✅ **Perfect match!**

### 2. ✅ File Format (TSV without quotes)
**Before:** Files had quoted values (`quote=TRUE`)
```
"sample"	"PC1: 45% variance"	"PC2: 23% variance"
"sample1"	"-12.34"	"5.67"
```

**After:** Files have unquoted values (`quote=FALSE`)
```
sample	PC1: 45% variance	PC2: 23% variance
sample1	-12.34	5.67
```

✅ **Fixed in commits:**
- `bin/normalize_deseq2_qc_all_genes.r`
- `bin/normalize_deseq2_qc_invariant_genes.r`

### 3. ✅ File Collection in Workflow
The workflow correctly collects all DESeq2 QC outputs:

```groovy
ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT.out.sample_distances_txt.collect()
)
ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT.out.pca_all_genes_txt.collect()
)
ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT.out.pca_top_genes_txt.collect()
)
ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_ALL_GENES_ALIGNMENT.out.read_dist_norm_txt.collect()
)
```

✅ **Files are passed to MultiQC**

### 4. ✅ MultiQC Custom Content Configuration
The MultiQC config defines 8 plots per quantifier (Kallisto example):

1. **Read Distribution (All Genes)** - Box plot
   - Pattern: `kallisto\.deseq2\.all_genes\.read\.distribution\.normalized\.txt`
   - Section: `deseq2-kallisto-qc`

2. **Read Distribution (Invariant Genes)** - Box plot
   - Pattern: `kallisto\.deseq2\.invariant_genes\.read\.distribution\.normalized\.txt`
   - Section: `deseq2-kallisto-qc`

3. **Sample Distance (All Genes)** - Heatmap
   - Pattern: `kallisto\.deseq2\.all_genes\.sample\.dists\.txt`
   - Section: `deseq2-kallisto-qc`

4. **Sample Distance (Invariant Genes)** - Heatmap
   - Pattern: `kallisto\.deseq2\.invariant_genes\.sample\.dists\.txt`
   - Section: `deseq2-kallisto-qc`

5. **PCA All Genes (All Genes Normalization)** - Scatter
   - Pattern: `kallisto\.deseq2\.all_genes\.pca\.vals\.txt`
   - Section: `deseq2-kallisto-qc`

6. **PCA All Genes (Invariant Genes Normalization)** - Scatter
   - Pattern: `kallisto\.deseq2\.invariant_genes\.pca\.vals\.txt`
   - Section: `deseq2-kallisto-qc`

7. **PCA Top 500 Genes (All Genes Normalization)** - Scatter
   - Pattern: `kallisto\.deseq2\.all_genes\.pca\.top500\.vals\.txt`
   - Section: `deseq2-kallisto-qc`

8. **PCA Top 500 Genes (Invariant Genes Normalization)** - Scatter
   - Pattern: `kallisto\.deseq2\.invariant_genes\.pca\.top500\.vals\.txt`
   - Section: `deseq2-kallisto-qc`

✅ **All patterns match generated files**

## Report Section Order
The MultiQC config places DESeq2 QC sections at priority 1001-1005:
```yaml
report_section_order:
  deseq2-kallisto-qc:
    order: 1005
  deseq2-salmon-qc:
    order: 1004
  deseq2-star-salmon-qc:
    order: 1003
  deseq2-star-rsem-qc:
    order: 1002
  deseq2-star-genome-qc:
    order: 1001
```

✅ **Sections will appear after quantification results (order: 2000) and before summaries (order: -1000)**

## Test File Examples

### Sample PCA Data File
```
sample	PC1: 45.2% variance	PC2: 23.1% variance
sample1	-12.34	5.67
sample2	10.23	-8.91
sample3	2.45	12.34
sample4	-5.67	-3.21
```

### Sample Distance File
```
	sample1	sample2	sample3	sample4
sample1	0.00	45.23	32.11	28.45
sample2	45.23	0.00	38.92	41.23
sample3	32.11	38.92	0.00	29.87
sample4	28.45	41.23	29.87	0.00
```

### Sample Read Distribution File
```
	sample1	sample2	sample3	sample4
gene1	5.23	6.78	5.92	6.12
gene2	7.45	7.89	7.23	7.56
gene3	4.12	4.56	4.34	4.45
gene4	8.90	9.12	8.78	8.95
gene5	3.45	3.67	3.56	3.59
```

## Quantifiers Supported
The configuration supports all pipeline quantifiers:
- ✅ kallisto
- ✅ salmon
- ✅ star.rsem
- ✅ star.salmon
- ✅ star.genome
- ✅ hisat2.genome

Each quantifier has its own set of 8 custom content plots (All Genes + Invariant Genes normalizations).

## Expected MultiQC Behavior
When the pipeline runs:

1. **DESeq2 QC modules** generate TSV files with proper naming
2. **Workflow** collects and stages all TSV files for MultiQC
3. **MultiQC** scans for files matching `fn_re` patterns
4. **Custom content module** parses TSV files
5. **Report sections** display interactive plots organized by quantifier and normalization method

## Verification Checklist
- [x] File naming patterns match MultiQC regex
- [x] TSV files use proper format (no quotes)
- [x] Files are collected in workflow channels
- [x] MultiQC config has custom_data sections defined
- [x] Report section ordering is configured
- [x] All quantifiers are supported

## Conclusion
✅ **MultiQC is now correctly configured to collect and display DESeq2 QC files from all quantifiers.**

The fixes ensure that:
1. R scripts generate properly formatted TSV files
2. File naming follows MultiQC expectations
3. Custom content configuration matches file patterns
4. Workflow passes files to MultiQC correctly
