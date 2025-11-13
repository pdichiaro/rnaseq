# Publishing Configuration for Per-Sample Scaling Factors

## Overview
Per-sample scaling factor files (`*_scaling_factor.txt`) are generated but **NOT published** to the output directory.

## Rationale
1. **Internal use only**: These files are used exclusively for Nextflow channel operations
2. **Avoid clutter**: Publishing dozens of small text files would clutter the output directory
3. **Consolidated output**: The `scaling_dat.txt` file already contains all scaling factors in one place
4. **Work directory sufficient**: Files in the work directory are accessible during workflow execution

## Implementation

### Location: `nextflow.config`

Both normalization process configurations include the filter:

```groovy
withName: '.*NORMALIZE_DESEQ2_QC_ALL_GENES.*' {
    publishDir = [
        // ... path configuration ...
        saveAs: { filename -> 
            if (filename.endsWith('_scaling_factor.txt')) {
                return null  // Exclude from publishing
            }
            return filename
        }
    ]
}

withName: '.*NORMALIZE_DESEQ2_QC_INVARIANT_GENES.*' {
    publishDir = [
        // ... path configuration ...
        saveAs: { filename -> 
            if (filename.endsWith('_scaling_factor.txt')) {
                return null  // Exclude from publishing
            }
            return filename
        }
    ]
}
```

## Published vs Unpublished Files

### ✅ PUBLISHED (in deseq2/all_genes or deseq2/invariant_genes directories):
- `scaling_dat.txt` - Consolidated scaling factors table
- `*_normalized_counts.txt` - Normalized count matrices
- `*_rlog_counts.txt` - rlog-transformed counts
- `*.RData` - DESeq2 R objects
- `Quality_Control/` - QC plots (PCA, heatmaps)
- `Read_Distribution/` - Read distribution plots

### ❌ NOT PUBLISHED (remain in work directory only):
- `scaling_factors/*_scaling_factor.txt` - Individual per-sample scaling factor files
- `versions.yml` - Version information (published centrally)
- `*.pca.vals.txt` - Raw PCA data (transformed versions published to multiqc_data)
- `*.sample.dists.txt` - Raw distance data (transformed versions published to multiqc_data)
- Other intermediate files excluded via saveAs filters

## Benefits

1. **Clean output**: Users see only meaningful analysis results
2. **Efficient storage**: No redundant small files
3. **Workflow functionality**: Channel operations work seamlessly with work directory files
4. **Debugging capability**: Files still exist in work directory if needed for troubleshooting

## Accessing Unpublished Files (if needed)

If you need to access the per-sample scaling factor files for debugging:

1. **Find the work directory**:
   ```bash
   # Look at the Nextflow log or .nextflow.log
   grep "NORMALIZE_DESEQ2" .nextflow.log | grep "Submitted"
   ```

2. **Navigate to work directory**:
   ```bash
   cd work/ab/cd1234567890abcdef1234567890abcd/
   ls scaling_factors/
   ```

3. **View scaling factor**:
   ```bash
   cat scaling_factors/sample1_scaling_factor.txt
   ```

## Configuration Override (Advanced)

If you want to publish these files for a specific run, you can override the configuration:

```groovy
// In a custom config file
process {
    withName: '.*NORMALIZE_DESEQ2_QC_ALL_GENES.*' {
        publishDir = [
            path: "${params.outdir}/deseq2/all_genes",
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename }  // Publish everything
        ]
    }
}
```

Then run with:
```bash
nextflow run main.nf -c custom.config ...
```

However, this is **not recommended** as it will clutter your output directory.
