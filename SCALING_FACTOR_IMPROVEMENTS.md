# Scaling Factor Handling Improvements

## Summary
Improved the robustness of scaling factor handling for DeepTools BigWig normalization by implementing per-sample scaling factor files. This eliminates fragile grep-based parsing and makes the workflow more resilient to sample naming edge cases.

## Problem Addressed
Previously, the workflow used a single `scaling_dat.txt` file with all scaling factors, requiring grep parsing to extract the correct value for each sample. This approach had several issues:
- **Fragile**: Prone to failures with special characters or similar sample names
- **Error-prone**: grep could match partial sample names incorrectly
- **Not Nextflow-native**: Required bash text manipulation instead of channel operations

## Solution Implemented: Per-Sample Scaling Factor Files

### Overview
Each DESeq2 normalization now generates individual scaling factor files for each sample in a `scaling_factors/` directory:
- **File naming**: `{sample_name}_scaling_factor.txt`
- **File content**: Just the numeric scaling factor value (e.g., `0.8956234`)
- **Benefits**: 
  - Direct file-to-sample mapping via Nextflow channels
  - No text parsing required
  - Immune to sample naming issues
  - More maintainable and testable

### Changes Made

#### 1. R Scripts Updated
**Files modified:**
- `rnaseq/bin/normalize_deseq2_qc_all_genes.r`
- `rnaseq/bin/normalize_deseq2_qc_invariant_genes.r`

**Changes:**
Added code after writing `scaling_dat.txt` to create per-sample files:
```r
# Create per-sample scaling factor files for safe DeepTools normalization
# This avoids grep parsing issues and is more robust
cat("Creating per-sample scaling factor files...\n")
scaling_factors_dir <- "scaling_factors"
dir.create(scaling_factors_dir, showWarnings = FALSE, recursive = TRUE)

for (i in 1:nrow(scaling_dat_early)) {
    sample_name <- scaling_dat_early$sample[i]
    scaling_factor <- scaling_dat_early$scaling_factor[i]
    
    # Write individual scaling factor file (just the numeric value)
    sample_file <- file.path(scaling_factors_dir, paste0(sample_name, "_scaling_factor.txt"))
    writeLines(as.character(scaling_factor), sample_file)
    
    cat("  - Created:", sample_file, "with scaling factor:", scaling_factor, "\n")
}
```

**Note:** The existing `scaling_dat.txt` file is retained for backward compatibility and MultiQC reporting.

#### 2. Module Outputs Updated
**Files modified:**
- `rnaseq/modules/local/normalize_deseq2_qc_all_genes/main.nf`
- `rnaseq/modules/local/normalize_deseq2_qc_invariant_genes/main.nf`

**Changes:**
Added new output channel for individual scaling factor files:
```groovy
output:
    path "scaling_dat.txt"       , emit: scaling_factors
    path "scaling_factors/*_scaling_factor.txt", emit: scaling_factors_individual
    // ... other outputs
```

#### 3. Workflow Channel Logic Updated
**File modified:**
- `rnaseq/workflows/rnaseq/main.nf`

**Changes:**

a) **Added new channel initialization:**
```groovy
ch_scaling_factors_individual = Channel.empty()
```

b) **Updated normalization sections to emit individual files:**
```groovy
ch_normalization_scaling_factors_individual_pseudo = Channel.empty()
// ... then mix in the individual scaling factor outputs
ch_scaling_factors_individual = ch_scaling_factors_individual.mix(ch_normalization_scaling_factors_individual_pseudo)
```

c) **Replaced CSV parsing with file-based channel operations:**

**BEFORE (invariant_genes example):**
```groovy
ch_scaling_factors_invariant = ch_scaling_factors
    .collect()
    .map { files -> /* find invariant file */ }
    .first()

ch_scaling_per_sample_invariant = ch_scaling_factors_invariant
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row -> [row.sample, row.size_factor] }

ch_combined_input_invariant = ch_bam_for_deeptools_invariant
    .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
    .join(ch_scaling_per_sample_invariant, by: 0)
    .map { sample_id, meta, bam, bai, scaling_factor -> 
        [meta + [scaling_factor: scaling_factor], bam, bai] 
    }
```

**AFTER:**
```groovy
ch_scaling_per_sample_invariant = ch_scaling_factors_individual
    .flatten()
    .filter { file ->
        // Filter for invariant-specific files
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file ->
        // Extract sample name from filename: {sample_name}_scaling_factor.txt
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        [sample_name, file]
    }

ch_combined_input_invariant = ch_bam_for_deeptools_invariant
    .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
    .join(ch_scaling_per_sample_invariant, by: 0)
    .map { sample_id, meta, bam, bai, scaling_factor_file -> 
        [meta + [scaling_factor_file: scaling_factor_file], bam, bai] 
    }
```

**Key improvements:**
- Direct file path passed to meta map instead of numeric value
- No CSV parsing required
- File filtering based on directory structure (invariant vs all_genes)
- Sample name extracted from filename

d) **Same pattern applied to all_genes method:**
The all_genes normalization uses the same approach but filters for files NOT in invariant directories.

#### 4. DEEPTOOLS_BIGWIG_NORM Module Simplified
**File modified:**
- `rnaseq/modules/local/deeptools_bw_norm/main.nf`

**Changes:**

a) **Removed scaling_factors input parameter:**
```groovy
// BEFORE
input:
    tuple val(meta), path(bam), path(bai)
    path scaling_factors

// AFTER
input:
    tuple val(meta), path(bam), path(bai)
```

b) **Updated all script sections (6 variants) to read from file:**
```bash
# BEFORE
scaling=$(grep "^${meta.id}" $scaling_factors | cut -f2)
if [ -z "$scaling" ]; then
    scaling=1.0
fi

# AFTER
scaling=$(cat ${meta.scaling_factor_file})
if [ -z "$scaling" ]; then
    scaling=1.0
fi
```

c) **Added scaling factor file logging:**
```bash
echo "Scaling factor file: ${meta.scaling_factor_file}"
echo "Scaling factor: $scaling"
```

**Script variants updated:**
1. Unstranded libraries
2. Single-end forward-stranded libraries
3. Single-end reverse-stranded libraries
4. Paired-end forward-stranded libraries
5. Paired-end reverse-stranded libraries

#### 5. Publishing Configuration Updated
**File modified:**
- `rnaseq/nextflow.config`

**Changes:**
Updated both `NORMALIZE_DESEQ2_QC_ALL_GENES` and `NORMALIZE_DESEQ2_QC_INVARIANT_GENES` process configurations to exclude per-sample scaling factor files from publishing:

```groovy
withName: '.*NORMALIZE_DESEQ2_QC_ALL_GENES.*' {
    publishDir = [
        // ... path configuration ...
        saveAs: { filename -> 
            // ... other exclusions ...
            // NEW: Exclude individual per-sample scaling factor files
            if (filename.endsWith('_scaling_factor.txt')) {
                return null  // Don't publish these files
            }
            return filename
        }
    ]
}
```

**Effect:**
- Per-sample files are generated and used in workflow channels
- Files remain in work directory for DeepTools normalization
- Output directory stays clean without dozens of small scaling factor files
- Only the consolidated `scaling_dat.txt` is published

## Testing Recommendations

1. **Basic functionality**: Test with standard sample names
2. **Edge cases**: Test with samples containing:
   - Underscores in names
   - Hyphens in names
   - Similar prefixes (e.g., "sample1", "sample10", "sample1_rep2")
   - Special characters (if permitted by your workflow)
3. **Both normalization methods**: Test with `invariant_genes` and `all_genes`
4. **Multiple aligners**: Test with STAR, HISAT2, and Kallisto/Salmon pseudo-alignments

## Publishing Behavior

The per-sample scaling factor files are **NOT published** to the output directory:
- **Work directory only**: Files remain in the Nextflow work directory for channel operations
- **Not cluttering output**: The `saveAs` clause in `nextflow.config` excludes `*_scaling_factor.txt` files
- **Published files**: Only `scaling_dat.txt` is published (for backward compatibility and reporting)

This keeps the output directory clean while maintaining full functionality for the workflow.

## Backward Compatibility

- **Retained**: The `scaling_dat.txt` file is still generated and published for backward compatibility
- **MultiQC**: Existing MultiQC integration unchanged
- **Added**: New per-sample files are additional outputs used only internally, not replacements

## Benefits

1. **Robustness**: Eliminates grep parsing issues
2. **Clarity**: Direct file-to-sample mapping is more intuitive
3. **Maintainability**: Easier to debug and test
4. **Nextflow-native**: Uses channels and file operations instead of bash text parsing
5. **Scalability**: Better performance with large sample sets (no CSV parsing overhead)

## Files Changed Summary

### R Scripts (2 files)
- `rnaseq/bin/normalize_deseq2_qc_all_genes.r` - Added per-sample file generation
- `rnaseq/bin/normalize_deseq2_qc_invariant_genes.r` - Added per-sample file generation

### Modules (3 files)
- `rnaseq/modules/local/normalize_deseq2_qc_all_genes/main.nf` - Added scaling_factors_individual output
- `rnaseq/modules/local/normalize_deseq2_qc_invariant_genes/main.nf` - Added scaling_factors_individual output
- `rnaseq/modules/local/deeptools_bw_norm/main.nf` - Simplified to read from per-sample files

### Workflows (1 file)
- `rnaseq/workflows/rnaseq/main.nf` - Updated channel logic to use per-sample files

### Configuration (1 file)
- `rnaseq/nextflow.config` - Added saveAs filter to exclude per-sample scaling factor files from publishing

**Total: 7 files modified**

## Future Improvements

Potential enhancements could include:
1. Add validation to ensure scaling factor is numeric and within expected range
2. Add process to verify all samples have scaling factor files before running DeepTools
3. Create MultiQC module to display per-sample scaling factors in report
4. Add optional parameter to use original CSV-based approach if needed
