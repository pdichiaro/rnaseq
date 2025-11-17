# Implementation Summary: DESeq2 Normalization Strategy for BigWig Generation

## Overview
This document summarizes the changes made to integrate the `mmrnaseq` normalization strategy into the `nf-core/rnaseq` pipeline. The goal was to make the BigWig normalization process consistent with the `mmrnaseq` approach, where scaling factors are computed once by `DESEQ2_QC_SALMON` and then reused for BigWig generation.

## Key Changes

### 1. Workflow Integration (`workflows/rnaseq/main.nf`)

#### Parameter Updates
- Modified `normalization_method` parameter handling to accept multiple normalization methods
- Changed from single string to list of methods: `'invariant_genes,all_genes'` or individual values
- Updated logic to support both legacy single-method and new multi-method syntax

#### Channel Handling for Scaling Factors

**Previous approach (File-based):**
```groovy
ch_scaling_per_sample = ch_scaling_factors
    .map { meta, file -> [meta.id, file] }
    
ch_combined_input = ch_bam
    .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
    .join(ch_scaling_per_sample, by: 0)
    .map { sample_id, meta, bam, bai, scaling_factor_file -> 
        [meta + [scaling_factor_file: scaling_factor_file], bam, bai] 
    }
```

**New approach (Value-based, matching mmrnaseq):**
```groovy
ch_scaling_per_sample = ch_scaling_factors_individual
    .flatten()
    .filter { file ->
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file ->
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        def scaling_value = file.text.trim()  // Read the value directly
        [sample_name, scaling_value]
    }

ch_combined_input = ch_bam
    .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
    .join(ch_scaling_per_sample, by: 0)
    .map { sample_id, meta, bam, bai, scaling -> 
        [meta, bam, bai, scaling]  // Pass scaling value instead of file
    }
```

#### Key Differences:
1. **File vs Value**: The new approach reads the scaling factor value from the file in the channel operation and passes it as a value to the process
2. **Meta map**: No longer adding `scaling_factor_file` to the meta map
3. **Process input**: The process now receives a `val(scaling)` instead of accessing a file path from meta

### 2. Process Module Updates (`modules/local/deeptools_bw_norm/main.nf`)

#### Input Signature Change
**Before:**
```groovy
input:
tuple val(meta), path(bam), path(bai)
```

**After:**
```groovy
input:
tuple val(meta), path(bam), path(bai), val(scaling)
```

#### Script Changes
**Before (reading from file in meta):**
```bash
SCALING_FACTOR=$(cat "${meta.scaling_factor_file}")

bamCoverage \\
    --scaleFactor $SCALING_FACTOR \\
    ...
```

**After (using passed value):**
```bash
echo "Scaling factor: $scaling"

bamCoverage \\
    --scaleFactor $scaling \\
    ...
```

All script blocks (unstranded, single-end forward/reverse, paired-end forward/reverse) were updated with this pattern.

### 3. Separation of Normalization Methods

The workflow now properly handles two separate normalization paths:

1. **Invariant Genes Normalization:**
   - Uses scaling factors from `all_genes/invariant_genes/*_scaling_factor.txt`
   - Filtered by checking if parent or grandparent directory contains 'invariant'
   - Output: `DEEPTOOLS_BIGWIG_NORM_INVARIANT`

2. **All Genes Normalization (default):**
   - Uses scaling factors from `all_genes/*_scaling_factor.txt` (excluding invariant subdirectory)
   - Filtered by checking if parent directory does NOT contain 'invariant'
   - Output: `DEEPTOOLS_BIGWIG_NORM_ALL_GENES`

## Technical Implementation Details

### Directory Structure for Scaling Factors
```
results/deseq2/all_genes/
├── sample1_scaling_factor.txt
├── sample2_scaling_factor.txt
├── ...
└── invariant_genes/
    ├── sample1_scaling_factor.txt
    ├── sample2_scaling_factor.txt
    └── ...
```

### Channel Filtering Logic
```groovy
// For invariant genes
.filter { file ->
    def parent_dir = file.getParent()?.getName() ?: ""
    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
    parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
}

// For all genes (excluding invariant)
.filter { file ->
    def parent_dir = file.getParent()?.getName() ?: ""
    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
    !(parent_dir.contains('invariant') || grandparent_dir.contains('invariant'))
}
```

### Sample Name Extraction
```groovy
def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
```

This regex removes the `_scaling_factor.txt` suffix to get the sample ID for channel joining.

## Benefits of This Approach

1. **Consistency with mmrnaseq**: The BigWig generation now uses the same normalization methodology as the mmrnaseq pipeline
2. **Computational Efficiency**: Scaling factors are computed once by DESeq2 and reused, rather than being recomputed
3. **Flexibility**: Supports multiple normalization methods simultaneously
4. **Maintainability**: Clear separation between invariant genes and all genes approaches
5. **Simpler Process Logic**: The process receives a clean value rather than having to read from a file

## Testing Recommendations

1. **Syntax Validation**: ✅ Completed - `nextflow run . --help` succeeds
2. **Parameter Parsing**: ✅ Completed - `--normalization_method` and `--skip_deeptools_norm` appear in help
3. **Full Pipeline Test**: Requires compute environment with:
   - RNA-seq samples
   - Reference genome and annotation
   - `--normalization_method 'invariant_genes,all_genes'` (or specific method)
   - `--skip_deeptools_norm false`

## Migration Notes

For users upgrading from previous versions:

1. **Parameter Syntax**: The `--normalization_method` parameter now accepts comma-separated lists:
   - Old: `--normalization_method invariant_genes`
   - New: `--normalization_method 'invariant_genes,all_genes'` (or single method still works)

2. **Output Structure**: BigWigs are now organized by normalization method:
   - `results/bigwig/deeptools_norm/invariant_genes/*.bw`
   - `results/bigwig/deeptools_norm/all_genes/*.bw`

3. **Scaling Factor Source**: Scaling factors are now sourced from DESeq2 QC output rather than being computed during BigWig generation

## References

- **mmrnaseq pipeline**: https://github.com/Deena-G/mmrnaseq
- **Original PR discussion**: Focused on aligning normalization strategies between pipelines
- **nf-core/rnaseq**: https://github.com/nf-core/rnaseq

## Files Modified

1. `workflows/rnaseq/main.nf` - Main workflow logic and channel operations
2. `modules/local/deeptools_bw_norm/main.nf` - Process definition and script

## Commit Message Suggestion

```
feat: integrate mmrnaseq normalization strategy for BigWig generation

- Modified DEEPTOOLS_BIGWIG_NORM process to receive scaling factor as value
- Updated channel operations to read scaling factors from files and pass as values
- Added support for multiple normalization methods (invariant_genes, all_genes)
- Separated normalization paths for invariant genes and all genes
- Aligned with mmrnaseq approach: compute scaling factors once in DESeq2, reuse for BigWigs
- Improved computational efficiency by avoiding redundant calculations

This change makes the pipeline more consistent with the mmrnaseq workflow
and provides better flexibility for different normalization strategies.
```
