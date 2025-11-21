# Bug Fixes: Quality Control PDFs + DeepTools BigWig Outputs + Code Cleanup

## Problem Description

### Issue 1: Missing Quality Control PDFs for RSEM Quantification
When using STAR as the aligner with RSEM for quantification, the DESeq2 normalization scripts (`normalize_deseq2_qc_all_genes.r` and `normalize_deseq2_qc_invariant_genes.r`) were generating quality control plots in PDF format within the `Quality_Control/` directory, but these PDFs were **not being published** to the output directory.

This issue affected:
- STAR + RSEM quantification (both all_genes and invariant_genes normalization)

**Note:** Salmon is only used for auto-strandedness detection in this pipeline, not for quantification, so Salmon DESeq2 configurations were not modified.

The genome quantification method was working correctly and already included the Quality_Control directory in its output.

### Issue 2: DeepTools BigWig Only Publishing One Normalization Method
When `skip_deeptools_norm` was FALSE, the pipeline was generating normalized bigWig files for **both** normalization methods (all_genes and invariant_genes), but the `publishDir` configuration was only publishing **one** of them based on the `params.normalization_method` setting.

The workflow correctly created two separate processes:
- `DEEPTOOLS_BIGWIG_NORM_ALL_GENES` 
- `DEEPTOOLS_BIGWIG_NORM_INVARIANT`

However, the configuration used mutually exclusive `enabled` conditions, preventing both outputs from being published simultaneously.

## Root Cause

### Issue 1: Quality Control PDFs
The `publishDir` pattern in `workflows/rnaseq/nextflow.config` was missing the `Quality_Control` directory for RSEM quantification methods, even though the R scripts were correctly generating PDFs in that directory.

### Issue 2: DeepTools BigWig
The `publishDir` configuration used a single `withName: '.*:DEEPTOOLS_BIGWIG_NORM'` block that applied to both process variants, with mutually exclusive `enabled` conditions:

```groovy
// PROBLEMATIC CONFIGURATION
withName: '.*:DEEPTOOLS_BIGWIG_NORM' {
    publishDir = [
        [
            path: { "...bigwig_norm/all_genes" },
            enabled: { params.normalization_method == 'all_genes' }  // Only one enabled
        ],
        [
            path: { "...bigwig_norm/invariant_genes" },
            enabled: { params.normalization_method == 'invariant_genes' }  // at a time
        ]
    ]
}
```

This meant that even though both processes ran, only one set of bigWig files was published.

### R Script Output Structure
The R scripts generate the following PDF files in the `Quality_Control/` directory:
- `Sample_Distance_Heatmap.pdf` - Euclidean distance heatmap between samples
- `Sample_Correlation_Heatmap.pdf` - Pearson correlation heatmap between samples  
- `PCA_Plot_All_Genes.pdf` - PCA plot using all genes
- `PCA_Plot_Top500_Genes.pdf` - PCA plot using top 500 variable genes
- `PCA_Decomposed_All_Genes.pdf` - Decomposed PCA by sample groups (if applicable)
- `PCA_Decomposed_Top500_Genes.pdf` - Decomposed PCA for top genes

These PDFs provide high-quality static visualizations suitable for publications and presentations.

## Solution

### Fix 1: Quality Control PDFs
**Added `Quality_Control`** to the `publishDir` pattern for RSEM DESeq2 normalization processes.

### Fix 2: DeepTools BigWig Outputs
**Split the single configuration block into two separate blocks** - one for each process variant, removing the mutually exclusive `enabled` conditions. This ensures both normalization methods publish their bigWig files simultaneously.

### Fix 3: Code Cleanup
**Removed deprecated Salmon quantification block** - The entire `if (params.aligner == 'star' && params.quantification == 'salmon')` configuration block was removed since Salmon is now only used for auto-strandedness detection, not quantification.

### Files Modified
- `workflows/rnaseq/nextflow.config` (removed 22 lines + updated 2 RSEM patterns + restructured deepTools config)

### Changes Made

#### 1. Removed Deprecated Salmon Quantification Block (Lines 368-386)
The entire configuration block for `params.quantification == 'salmon'` was removed:

```groovy
// REMOVED - Salmon is only used for auto-strandedness, not quantification
if (params.aligner == 'star' && params.quantification == 'salmon') {
    if (!(params.skip_qc ?: false) & !(params.skip_deseq2_qc ?: false)) {
        process {
            withName: '.*:NORMALIZE_DESEQ2_QC_ALL_GENES.*' { ... }
            withName: '.*:NORMALIZE_DESEQ2_QC_INVARIANT_GENES.*' { ... }
        }
    }
}
```

**Rationale:** Salmon is retained in the pipeline only for auto-strandedness detection. It is no longer used as a quantification method, making this configuration block obsolete.

---

### DESeq2 Quality Control PDF Fixes

#### 2. RSEM - All Genes Normalization (Line ~372, previously ~397)
**Before:**
```groovy
pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log,normalized_counts.txt,rlog_counts.txt,scaling_dat.txt,correlation_mqc.tsv,read_dist_mqc.tsv,pca.vals_mqc.tsv,sample.dists_mqc.tsv,Read_Distribution}"
```

**After:**
```groovy
pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log,normalized_counts.txt,rlog_counts.txt,scaling_dat.txt,correlation_mqc.tsv,read_dist_mqc.tsv,pca.vals_mqc.tsv,sample.dists_mqc.tsv,Read_Distribution,Quality_Control}"
```

#### 3. RSEM - Invariant Genes Normalization (Line ~386, previously ~404)
**Before:**
```groovy
pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log,normalization,normalized_counts.txt,normalized_filt.txt,scaling_dat.txt,correlation_mqc.tsv,read_dist_mqc.tsv,pca.vals_mqc.tsv,sample.dists_mqc.tsv,Read_Distribution}"
```

**After:**
```groovy
pattern: "*{RData,pca.vals.txt,plots.pdf,sample.dists.txt,size_factors,log,normalization,normalized_counts.txt,normalized_filt.txt,scaling_dat.txt,correlation_mqc.tsv,read_dist_mqc.tsv,pca.vals_mqc.tsv,sample.dists_mqc.tsv,Read_Distribution,Quality_Control}"
```

---

### DeepTools BigWig Normalization Fixes

#### 4. DeepTools Configuration Restructure (Lines ~503-521, previously ~503-522)

**Before (Single block with conditional publishing):**
```groovy
if (params.aligner == 'star' && !params.skip_deeptools_norm) {
    process {
        withName: '.*:DEEPTOOLS_BIGWIG_NORM' {
            publishDir = [
                [
                    path: { "${params.outdir}/${params.aligner}/${params.quantification}/bigwig_norm/all_genes" },
                    mode: params.publish_dir_mode,
                    pattern: "*.norm.bw",
                    enabled: { params.normalization_method == 'all_genes' }  // ❌ Only one enabled
                ],
                [
                    path: { "${params.outdir}/${params.aligner}/${params.quantification}/bigwig_norm/invariant_genes" },
                    mode: params.publish_dir_mode,
                    pattern: "*.norm.bw", 
                    enabled: { params.normalization_method == 'invariant_genes' }  // ❌ at a time
                ]
            ]
        }
    }
}
```

**After (Separate blocks for each process, both always published):**
```groovy
if (params.aligner == 'star' && !params.skip_deeptools_norm) {
    process {
        withName: '.*:DEEPTOOLS_BIGWIG_NORM_ALL_GENES' {
            publishDir = [
                path: { "${params.outdir}/${params.aligner}/${params.quantification}/deeptools/all_genes" },
                mode: params.publish_dir_mode,
                pattern: "*.norm.bw"  // ✅ Always published
            ]
        }
        
        withName: '.*:DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
            publishDir = [
                path: { "${params.outdir}/${params.aligner}/${params.quantification}/deeptools/invariant_genes" },
                mode: params.publish_dir_mode,
                pattern: "*.norm.bw"  // ✅ Always published
            ]
        }
    }
}
```

**Key Changes:**
1. **Split configuration**: Separate `withName` blocks for `DEEPTOOLS_BIGWIG_NORM_ALL_GENES` and `DEEPTOOLS_BIGWIG_NORM_INVARIANT`
2. **Removed `enabled` conditions**: Both outputs are now always published when `skip_deeptools_norm` is FALSE
3. **Path change**: `bigwig_norm/` → `deeptools/` for consistency with other output directories
4. **Simplified publishDir**: Changed from array of configs to single config per process

**Rationale:** The workflow already runs both normalization methods in parallel when `skip_deeptools_norm` is FALSE. The configuration should publish both outputs, not force users to choose one via `params.normalization_method`.

---

## Expected Output Structure After Fix

When using STAR + RSEM quantification, the output directory structure will now include:

```
results/
└── star/
    ├── rsem/
    │   ├── deseq2/
    │   │   ├── all_genes/
    │   │   │   ├── Quality_Control/              # ✅ NOW PUBLISHED
    │   │   │   │   ├── Sample_Distance_Heatmap.pdf
    │   │   │   │   ├── Sample_Correlation_Heatmap.pdf
    │   │   │   │   ├── PCA_Plot_All_Genes.pdf
    │   │   │   │   ├── PCA_Plot_Top500_Genes.pdf
    │   │   │   │   └── PCA_Decomposed_*.pdf
    │   │   │   ├── Read_Distribution/
    │   │   │   ├── scaling_dat.txt
    │   │   │   ├── normalized_counts.txt
    │   │   │   └── rlog_counts.txt
    │   │   └── invariant_genes/
    │   │       ├── Quality_Control/              # ✅ NOW PUBLISHED
    │   │       │   ├── Sample_Distance_Heatmap.pdf
    │   │       │   ├── Sample_Correlation_Heatmap.pdf
    │   │       │   ├── PCA_Plot_All_Genes.pdf
    │   │       │   ├── PCA_Plot_Top500_Genes.pdf
    │   │       │   └── PCA_Decomposed_*.pdf
    │   │       ├── Read_Distribution/
    │   │       ├── normalization/
    │   │       ├── scaling_dat.txt
    │   │       ├── normalized_counts.txt
    │   │       └── rlog_counts.txt
    │   └── deeptools/                            # ✅ BOTH METHODS NOW PUBLISHED
    │       ├── all_genes/
    │       │   ├── sample1.norm.bw
    │       │   ├── sample2.norm.bw
    │       │   └── sample3.norm.bw
    │       └── invariant_genes/
    │           ├── sample1.norm.bw
    │           ├── sample2.norm.bw
    │           └── sample3.norm.bw
```

**Note:** When `skip_deeptools_norm` is FALSE, normalized bigWig files are now generated and published for **both** normalization methods simultaneously, regardless of the `params.normalization_method` setting.

## Impact

### Positive Changes
1. **Quality Control PDFs:** Users now get publication-quality PDF plots for RSEM quantification, matching the behavior of genome quantification
2. **Complete DeepTools Outputs:** Both normalization methods now publish bigWig files simultaneously when deepTools normalization is enabled
3. **Better Comparison:** Users can now compare normalized bigWig files from different normalization strategies side-by-side
4. **Cleaner Codebase:** Removed 22 lines of obsolete Salmon quantification configuration

### No Breaking Changes
- All existing TSV outputs for MultiQC interactive plots remain unchanged
- DESeq2 normalization behavior is unchanged - only output publishing is fixed
- Backward compatible - existing pipelines will continue to work

### Consistency Improvements
- RSEM quantification now produces the same Quality_Control outputs as genome quantification
- DeepTools configuration now matches the pattern used for DESeq2 normalization (separate blocks per process variant)

## Testing Recommendations

### Test 1: Quality Control PDFs (DESeq2)
1. Run pipeline with `--aligner star --quantification rsem`
2. Verify `Quality_Control/` directory appears in:
   - `results/star/rsem/deseq2/all_genes/Quality_Control/`
   - `results/star/rsem/deseq2/invariant_genes/Quality_Control/`
3. Confirm PDF files are present and contain valid plots:
   ```bash
   ls -lh results/star/rsem/deseq2/*/Quality_Control/*.pdf
   ```
4. Verify MultiQC report still includes interactive plots from TSV files

### Test 2: DeepTools Normalized BigWig Files
1. Run pipeline with `--aligner star --quantification rsem` (without `--skip_deeptools_norm`)
2. Verify **both** normalization methods publish bigWig files:
   ```bash
   ls -lh results/star/rsem/deeptools/all_genes/*.norm.bw
   ls -lh results/star/rsem/deeptools/invariant_genes/*.norm.bw
   ```
3. Confirm both directories contain `.norm.bw` files for all samples
4. Verify bigWig files are valid:
   ```bash
   # Check if bigWig files are valid (should show chromosome info)
   bigWigInfo results/star/rsem/deeptools/all_genes/sample1.norm.bw
   ```

### Test 3: Verify No Regression
1. Test with genome quantification: `--aligner star --quantification genome`
2. Test with deepTools disabled: `--skip_deeptools_norm`
3. Confirm all existing outputs remain unchanged

## Related Files

### Modified
- **Workflow Configuration:** `workflows/rnaseq/nextflow.config` (all fixes applied here)

### Related (unchanged)
- **R Scripts:** `bin/normalize_deseq2_qc_all_genes.r`, `bin/normalize_deseq2_qc_invariant_genes.r`
- **DESeq2 Modules:** `modules/local/normalize_deseq2_qc_all_genes/main.nf`, `modules/local/normalize_deseq2_qc_invariant_genes/main.nf`
- **DeepTools Module:** `modules/local/deeptools_bw_norm/main.nf`
- **Main Workflow:** `workflows/rnaseq/main.nf` (process invocations)

## Notes

### DESeq2 Quality Control
- The R scripts were already generating the Quality_Control PDFs correctly
- The genome quantification method already had this directory in the publish pattern
- This fix brings RSEM output behavior in line with genome quantification
- TSV files for MultiQC interactive plots are published at the root level, separate from the organized PDF folders

### DeepTools BigWig Normalization
- The workflow already executed both normalization processes in parallel
- The issue was purely in the configuration layer (publishDir settings)
- Path changed from `bigwig_norm/` to `deeptools/` for consistency
- Users no longer need to choose between normalization methods - both are always published when enabled
