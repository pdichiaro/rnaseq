# MultiQC Collection Debugging Guide

## Please run these commands and share the output:

### 1. Check if TSV files exist in results directory
```bash
find results -name "*.deseq2.*.txt" -type f 2>/dev/null
```

### 2. Check one sample file content (first 5 lines)
```bash
find results -name "*.pca.vals.txt" -type f 2>/dev/null | head -1 | xargs head -5
```

### 3. Check for quotes in files
```bash
find results -name "*.deseq2.*.txt" -type f 2>/dev/null | head -1 | xargs cat | head -3
```

### 4. Find the MultiQC work directory
```bash
find work -type d -name "*multiqc*" | grep -v "multiqc_data" | grep -v "multiqc_plots" | head -1
```

### 5. List files that were staged for MultiQC
```bash
# Replace PATH with the path from command 4
ls -la work/PATH_TO_MULTIQC_DIR/
```

### 6. Check what files match the patterns
```bash
cd work/PATH_TO_MULTIQC_DIR/
find . -name "*.deseq2.*.txt" -type f
```

### 7. Check the actual filenames being generated
```bash
find results/kallisto/deseq2 -name "*.txt" -type f 2>/dev/null
```

### 8. Check the MultiQC command that was executed
```bash
grep -A 5 "multiqc" .nextflow.log | tail -30
```

### 9. Check if files are being collected into the channel
```bash
grep "NORMALIZE_DESEQ2_QC" .nextflow.log | grep -i "publishing\|output" | tail -20
```

### 10. Most Important - Check the exact filenames
```bash
ls -1 results/kallisto/deseq2/all_genes/*.txt 2>/dev/null
ls -1 results/kallisto/deseq2/invariant_genes/*.txt 2>/dev/null
```

## What I'm Looking For:

1. Are the TSV files actually being created?
2. What are the EXACT filenames?
3. Are they reaching the MultiQC work directory?
4. Do the filenames match the regex patterns in multiqc_config.yml?

## Expected Filename Patterns:

According to the config, files should match:
- `kallisto\\.deseq2\\.all_genes\\.pca\\.vals\\.txt`
- `kallisto\\.deseq2\\.all_genes\\.sample\\.dists\\.txt`
- `kallisto\\.deseq2\\.all_genes\\.read\\.distribution\\.normalized\\.txt`
- `kallisto\\.deseq2\\.invariant_genes\\.pca\\.vals\\.txt`
- etc.

Are the actual files named correctly?
