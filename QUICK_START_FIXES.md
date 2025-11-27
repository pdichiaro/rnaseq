# Quick Start: Common Parameter Fixes

## 🚨 Two Most Common Issues (FIXED IN 2 MINUTES)

### Issue 1: Missing MultiQC Report

**Symptom:** No `multiqc_report.html` generated

**Cause:** Using invalid `--publish_dir` parameter

**Fix:**
```bash
# ❌ WRONG - Remove this line
--publish_dir $outdir \

# ✅ CORRECT - Only use --outdir
--outdir $outdir \
```

---

### Issue 2: Missing Invariant Genes Normalization

**Symptom:** No files in `star/genome/deseq2/invariant_genes/deeptools_normalize/`

**Cause:** The `--normalization_method` parameter requires both methods to be specified correctly

**Note:** Both formats work (the pipeline trims whitespace automatically):
```bash
# ✅ BOTH ARE CORRECT
--normalization_method 'all_genes,invariant_genes'
--normalization_method 'all_genes, invariant_genes'
```

**Common Issues:**
- Missing one of the methods
- Typo in method names
- Check your actual command to ensure both methods are specified

---

## ✅ Minimal Working Example

```bash
nextflow run pdichiaro/rnaseq \
    -r main \
    --input samplesheet.csv \
    --outdir results/ \
    --fasta genome.fa \
    --gtf annotation.gtf \
    --star_index star_idx/ \
    --gene_bed genes.bed \
    --normalization_method 'all_genes,invariant_genes' \
    -profile singularity \
    -resume
```

**Key Points:**
- ✅ No `--publish_dir`
- ✅ No spaces in comma-separated values
- ✅ Always use `-resume` for efficiency

---

## 📊 Verify Success

After running, check these files exist:

```bash
# MultiQC report
ls results/multiqc/star/multiqc_report.html

# Both normalization methods
ls results/star/genome/deseq2/all_genes/deeptools_normalize/
ls results/star/genome/deseq2/invariant_genes/deeptools_normalize/
```

---

## 🔗 More Information

For detailed troubleshooting, see: [TROUBLESHOOTING_COMMON_ISSUES.md](TROUBLESHOOTING_COMMON_ISSUES.md)
