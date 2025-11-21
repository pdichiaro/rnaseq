# Genome Quantification Merge Output Verification Plan

## Test Configuration
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner star \
  --quantification genome,rsem
```

This will run:
1. **STAR alignment**
2. **RSEM quantification** (transcriptome-based)
3. **Genome quantification** (genomic feature-based with featureCounts/Rsubread)

## Verification Checklist for Genome Quantification Merge

### 1. Per-Sample Organization (NEW STRUCTURE)
Check that per-sample count files are in individual folders:

```bash
# Should show sample folders
ls -lh results/star/genome/

# Expected output:
# drwxr-xr-x sample1/
# drwxr-xr-x sample2/
# drwxr-xr-x sample3/
# -rw-r--r-- genome_exon_counts_merged.txt
# -rw-r--r-- genome_transcript_counts_merged.txt
# ... (other merged files)
# drwxr-xr-x deseq2/
# drwxr-xr-x deeptools/
```

```bash
# Check individual sample folder contents
ls -lh results/star/genome/sample1/

# Expected files in each sample folder:
# sample1_exon_counts.txt
# sample1_transcript_counts.txt
# sample1_intron_counts.txt
# sample1_5utr_counts.txt
# sample1_3utr_counts.txt
# sample1_combined_counts.txt
# sample1_summary.txt
```

### 2. Merged Count Files at Root Level
Verify all merged count matrices exist at the genome root:

```bash
# List all merged files
ls -1 results/star/genome/*.txt

# Expected merged output files:
# 1. genome_exon_counts_merged.txt
# 2. genome_transcript_counts_merged.txt
# 3. genome_intron_counts_merged.txt
# 4. genome_5utr_counts_merged.txt
# 5. genome_3utr_counts_merged.txt
# 6. genome_exon_merge_summary.txt
# 7. genome_transcript_merge_summary.txt
# 8. genome_intron_merge_summary.txt
# 9. genome_5utr_merge_summary.txt
# 10. genome_3utr_merge_summary.txt
```

### 3. Merged File Content Verification

#### A. Check Matrix Dimensions
```bash
# Exon counts should have all samples as columns
head -n 5 results/star/genome/genome_exon_counts_merged.txt

# Expected format:
# Geneid    Chr    Start    End    Strand    Length    sample1    sample2    sample3
# ENSG001   chr1   1000     2000   +         1000      150        200        175
# ENSG002   chr1   3000     4000   -         1000      80         90         85
```

```bash
# Count columns (should be 6 fixed columns + N samples)
head -n 1 results/star/genome/genome_exon_counts_merged.txt | awk '{print NF}'

# For 3 samples: should print 9 (6 + 3)
# For 10 samples: should print 16 (6 + 10)
```

#### B. Verify All Samples Included
```bash
# Check header contains all sample names
head -n 1 results/star/genome/genome_exon_counts_merged.txt

# Verify each sample from samplesheet appears in header
```

#### C. Check Gene Coverage
```bash
# Count genes in merged file
wc -l results/star/genome/genome_exon_counts_merged.txt

# Compare to individual sample
wc -l results/star/genome/sample1/sample1_exon_counts.txt

# Should be identical (same genes, just more columns in merged)
```

### 4. Compare Feature Types

Check that all 5 feature types are properly merged:

```bash
# Quick comparison of file sizes
ls -lh results/star/genome/genome_*_counts_merged.txt

# Expected features (in order of typical size):
# 1. transcript (largest - includes UTRs and introns)
# 2. exon (protein coding regions)
# 3. intron (intronic regions)
# 4. 5utr (5' untranslated regions)
# 5. 3utr (3' untranslated regions)
```

### 5. Summary File Verification

Check merge summary files for quality metrics:

```bash
# View exon merge summary
cat results/star/genome/genome_exon_merge_summary.txt

# Expected content:
# - Sample names
# - Total reads assigned
# - Unassigned reads
# - Assignment rate %
```

```bash
# Check all summary files exist
ls -1 results/star/genome/genome_*_merge_summary.txt

# Should show 5 summary files (one per feature type)
```

### 6. DESeq2 QC Integration

Verify that DESeq2 QC uses the merged exon counts:

```bash
# Check DESeq2 output directory
ls -lh results/star/genome/deseq2/all_genes/

# Expected subdirectories:
# Quality_Control/
# Read_Distribution/
```

```bash
# Check Quality Control PDFs exist
ls -1 results/star/genome/deseq2/all_genes/Quality_Control/*.pdf

# Expected QC files:
# - Sample_correlation_heatmap.pdf
# - Sample_PCA_plot.pdf
# - Sample_MDS_plot.pdf
# - Sample_hierarchical_clustering.pdf
```

```bash
# Verify Read Distribution plots
ls -1 results/star/genome/deseq2/all_genes/Read_Distribution/*.pdf

# Expected files:
# - Read_distribution_barplot.pdf
# - Read_distribution_percentage.pdf
```

### 7. DeepTools BigWig Files

Check that normalized BigWig files are generated:

```bash
# Check for all_genes normalization
ls -lh results/star/genome/deeptools/all_genes/

# Expected: one .norm.bw file per sample
# sample1.norm.bw
# sample2.norm.bw
# sample3.norm.bw
```

```bash
# If invariant gene normalization was run
ls -lh results/star/genome/deeptools/invariant_genes/

# Expected: one .norm.bw file per sample
```

### 8. Comparison with RSEM Structure

Verify consistency between genome and RSEM output organization:

```bash
# Compare directory structures
tree -L 2 results/star/genome/
tree -L 2 results/star/rsem/

# Both should follow same pattern:
# ├── sample1/           (per-sample folder)
# ├── sample2/
# ├── sample3/
# ├── merged_*.txt       (merged at root)
# ├── deseq2/
# └── deeptools/
```

### 9. File Count Validation

```bash
# Count files in genome root (should be manageable)
ls results/star/genome/ | wc -l

# With 3 samples, expected:
# - 3 sample folders
# - 10 merged files (5 counts + 5 summaries)
# - 2 subdirectories (deseq2/, deeptools/)
# Total: ~15 items

# WITHOUT sample folders (old structure):
# - 21 per-sample files (7 × 3 samples)
# - 10 merged files
# - 2 subdirectories
# Total: ~33 items (much more cluttered!)
```

### 10. Data Integrity Checks

#### A. No Empty Files
```bash
# Check that merged files are not empty
for file in results/star/genome/genome_*_counts_merged.txt; do
    lines=$(wc -l < "$file")
    echo "$file: $lines lines"
    if [ "$lines" -lt 2 ]; then
        echo "⚠️  WARNING: File is too small!"
    fi
done
```

#### B. Consistent Gene IDs Across Features
```bash
# Extract gene IDs from exon counts
tail -n +2 results/star/genome/genome_exon_counts_merged.txt | cut -f1 | sort > /tmp/exon_genes.txt

# Extract gene IDs from transcript counts
tail -n +2 results/star/genome/genome_transcript_counts_merged.txt | cut -f1 | sort > /tmp/transcript_genes.txt

# Compare - should be identical
diff /tmp/exon_genes.txt /tmp/transcript_genes.txt
# No output = files are identical ✅
```

#### C. Count Sum Verification
```bash
# For one sample, verify counts are consistent
# Sum from individual file should match column sum in merged file

# Get sample1 total from individual file
tail -n +2 results/star/genome/sample1/sample1_exon_counts.txt | awk '{sum+=$7} END {print sum}'

# Get sample1 total from merged file (column 7 if it's first sample)
tail -n +2 results/star/genome/genome_exon_counts_merged.txt | awk '{sum+=$7} END {print sum}'

# These should match ✅
```

### 11. MultiQC Integration

Check that genome quantification is included in MultiQC report:

```bash
# Check MultiQC report exists
ls -lh results/multiqc/star_genome/multiqc_report.html

# Open report and verify sections:
# - featureCounts / Rsubread statistics
# - Assignment rates per sample
# - Read distribution across features
```

## Expected Issues to Watch For

### ⚠️ Potential Problems

1. **Missing Samples in Merge**
   - Check if all per-sample folders have output files
   - Verify merge includes all samples in header

2. **Feature Type Mismatch**
   - All 5 feature types should have same number of genes
   - Check if any feature failed to merge

3. **Sample Order Inconsistency**
   - Verify sample order is consistent across all merged files
   - DESeq2 requires consistent sample ordering

4. **Empty or Truncated Files**
   - Check file sizes are reasonable
   - Verify no files are 0 bytes

5. **Permission Issues**
   - Verify all files are readable
   - Check publishDir worked correctly

## Success Criteria

✅ **All checks pass if:**

1. Per-sample files are in `genome/{sample_id}/` folders
2. 10 merged files exist at `genome/` root (5 counts + 5 summaries)
3. All merged files contain all samples as columns
4. Gene IDs are consistent across feature types
5. Count sums match between individual and merged files
6. DESeq2 QC PDFs are generated successfully
7. DeepTools BigWig files exist for all samples
8. Structure matches RSEM organization pattern
9. No empty or truncated files
10. MultiQC includes genome quantification metrics

## Quick Verification Script

Save this as `verify_genome_merge.sh`:

```bash
#!/bin/bash

OUTDIR="results/star/genome"

echo "=== Genome Quantification Merge Verification ==="
echo

echo "1. Checking directory structure..."
ls -lh $OUTDIR/ | head -20

echo
echo "2. Checking per-sample folders..."
ls -d $OUTDIR/*/ | grep -v "deseq2\|deeptools" || echo "⚠️  No sample folders found!"

echo
echo "3. Checking merged count files..."
ls -1 $OUTDIR/genome_*_counts_merged.txt 2>/dev/null || echo "⚠️  No merged count files found!"

echo
echo "4. Checking merged summary files..."
ls -1 $OUTDIR/genome_*_merge_summary.txt 2>/dev/null || echo "⚠️  No summary files found!"

echo
echo "5. Checking sample1 folder contents..."
ls -lh $OUTDIR/sample1/ 2>/dev/null || echo "⚠️  sample1 folder not found!"

echo
echo "6. Verifying merged file dimensions..."
for file in $OUTDIR/genome_*_counts_merged.txt; do
    if [ -f "$file" ]; then
        lines=$(wc -l < "$file")
        cols=$(head -n 1 "$file" | awk '{print NF}')
        echo "$(basename $file): $lines genes × $cols columns"
    fi
done

echo
echo "7. Checking DESeq2 QC outputs..."
ls -1 $OUTDIR/deseq2/all_genes/Quality_Control/*.pdf 2>/dev/null | wc -l | xargs echo "QC PDFs found:"

echo
echo "8. Checking DeepTools BigWig files..."
ls -1 $OUTDIR/deeptools/all_genes/*.bw 2>/dev/null | wc -l | xargs echo "BigWig files found:"

echo
echo "=== Verification Complete ==="
```

Run with:
```bash
bash verify_genome_merge.sh
```

## Next Steps After Verification

Once all checks pass:

1. ✅ Confirm structure matches expectations
2. ✅ Validate with biological QC (PCA, correlation)
3. ✅ Compare with RSEM results for consistency
4. ✅ Proceed with downstream differential expression analysis
5. ✅ Document any observations or issues

## Notes

- This verification assumes the pipeline completed successfully
- All checks should be performed on actual pipeline output
- Any warnings should be investigated before downstream analysis
- Keep this document for future pipeline runs as a QC checklist
