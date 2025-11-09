# MultiQC Files Investigation Summary

## Problem
The custom gene count QC data from `normalize_deseq2_qc` and `normalize_deseq2_qc_invariant_genes` modules was not appearing in MultiQC reports.

## Investigation

### Files Generated
Both normalization modules generate several files in the working directory:

**normalize_deseq2_qc:**
- `{prefix}norm_genes_detected_mqc.txt` - Number of genes detected per sample
- `{prefix}norm_genes_detected_plot_mqc.txt` - Stacked bar plot showing genes detected by count threshold
- `{prefix}norm_sample_correlation_mqc.txt` - Sample-to-sample correlation heatmap
- `{prefix}norm_sample_distances_dendrogram_mqc.txt` - Hierarchical clustering dendrogram
- `{prefix}norm_sample_distances_heatmap_mqc.txt` - Sample distance heatmap
- `{prefix}norm_size_factors_mqc.txt` - DESeq2 size factors used for normalization
-  Normalized counts files (.tsv)
- DESeq2 object (.rds)
- Log files

**normalize_deseq2_qc_invariant_genes:**
- `{prefix}invariant_genes_mqc.txt` - Plot showing coefficient of variation of top 1000 genes
- Normalized counts files (.tsv)
- List of invariant genes
- Log files

### Key Findings

1. **File Location**: The `*_mqc.txt` files are generated in the process working directory
2. **Output Declaration**: Both modules properly declare `*_mqc.txt` files in their `output:` section
3. **MultiQC Discovery**: MultiQC uses the custom_content module which looks for files matching specific patterns

4. **Pattern Matching**: 
   - The original modules used comments like `# Just let MultiQC find the files where they are`
   - This suggests the files should be automatically discovered by MultiQC
   - MultiQC's custom_content module searches for files with `_mqc` in their name

5. **Potential Issues**:
   - Files may not be in the correct location when MultiQC runs
   - MultiQC might not be searching the normalization process output directories
   - The file format might not match what MultiQC's custom_content module expects

### What We Did

1. **Added Debug Output**: Modified both normalization modules to list all generated `*_mqc.txt` files immediately after the R script runs:
   ```bash
   echo "=== FILES GENERATED FOR MULTIQC ==="
   ls -la *.txt 2>/dev/null || echo "No .txt files found"
   echo "===================================="
   ```

2. **Workflow Output Connection**: Verified that:
   - Both modules output `multiqc_files` as a named tuple element
   - The workflow (`workflows/rnaseq.nf`) collects these files into `ch_multiqc_custom_config`
   - Files are mixed with other QC outputs before being passed to MULTIQC

### Expected Behavior

When the normalization processes run, you should see output like:
```
=== FILES GENERATED FOR MULTIQC ===
-rw-r--r-- 1 user user  1234 Nov  9 13:00 genome_deseq2norm_genes_detected_mqc.txt
-rw-r--r-- 1 user user  2345 Nov  9 13:00 genome_deseq2norm_genes_detected_plot_mqc.txt
...
====================================
```

### Next Steps for Debugging

1. **Run the pipeline** with alignment enabled (not `--skip_alignment`)
2. **Check process logs** for the debug output showing which files were generated
3. **Examine MultiQC inputs**:
   ```bash
   # Look at the MultiQC work directory
   ls -la work/*/*/MULTIQC/
   ```
4. **Verify file content**: Check that the `*_mqc.txt` files have the correct format:
   ```bash
   # Example: check one of the generated files
   cat work/*/*/NORMALIZE_DESEQ2_QC/*_mqc.txt | head
   ```

### File Format Requirements

MultiQC custom_content expects specific formats (from GeneralNormalizer R package):
- Tab-separated text files
- Special headers starting with `#` for plot configuration  
- Data section with appropriate column headers

Example header format:
```
# id: 'custom_data_section'  
# section_name: 'My Custom Data'
# plot_type: 'bargraph'
```

### Resolution Recommendations

If files still don't appear in MultiQC after running with the debug output:

1. **Check file generation**: Confirm files are created in the process work directory
2. **Check file format**: Verify the R scripts generate properly formatted custom_content files
3. **Check MultiQC collection**: Ensure the files make it to MultiQC's input channel
4. **Check MultiQC config**: Verify MultiQC is configured to search for custom_content files
5. **Check file naming**: Ensure `*_mqc.txt` pattern is consistent

### Files Modified

1. `modules/local/normalize_deseq2_qc/main.nf` - Added debug echo statements
2. `modules/local/normalize_deseq2_qc_invariant_genes/main.nf` - Added debug echo statements

## Conclusion

The modules are correctly set up to output MultiQC custom content files. The debug statements added will help identify:
- Whether files are actually being generated
- What the exact filenames are
- Whether there are any file I/O issues during generation

Run the pipeline with alignment enabled and check the process logs to see the debug output and confirm file generation.
