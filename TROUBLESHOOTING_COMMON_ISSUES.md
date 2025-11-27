# Troubleshooting Common Issues

This guide addresses common issues encountered when running the pdichiaro/rnaseq pipeline.

## Table of Contents
- [Missing MultiQC Report](#missing-multiqc-report)
- [Missing DeepTools Invariant Genes Normalization](#missing-deeptools-invariant-genes-normalization)
- [Parameter Issues](#parameter-issues)
- [Output Structure](#output-structure)

---

## Missing MultiQC Report

### Symptom
The MultiQC report is not generated at `$outdir/multiqc/star/multiqc_report.html`

### Common Causes

#### 1. Invalid `--publish_dir` Parameter ⚠️ CRITICAL

**Problem**: Using `--publish_dir` parameter which does not exist in this pipeline.

```bash
# ❌ INCORRECT - This parameter does not exist
nextflow run pdichiaro/rnaseq \
    --publish_dir $outdir \
    --outdir $outdir \
    ...
```

**Solution**: Remove the `--publish_dir` line completely. The pipeline uses only `--outdir`.

```bash
# ✅ CORRECT
nextflow run pdichiaro/rnaseq \
    --outdir $outdir \
    ...
```

**Why this matters**: The non-existent `--publish_dir` parameter interferes with the pipeline's internal publishDir logic, preventing proper file organization and MultiQC report generation.

#### 2. Pipeline Failure Before MultiQC Step

Check the pipeline logs to ensure all upstream processes completed successfully:

```bash
# Check execution trace
less $outdir/pipeline_info/execution_trace*.txt

# Check Nextflow log
less .nextflow.log
```

#### 3. Missing Input Data

Ensure all required QC processes have outputs to aggregate:

```bash
# Verify star alignment outputs exist
ls -lh $outdir/star/

# Check for individual QC outputs
find $outdir -name "*.html" -type f
```

### Verification

After fixing, verify the MultiQC report exists:

```bash
ls -lh $outdir/multiqc/star/multiqc_report.html
```

---

## Missing DeepTools Invariant Genes Normalization

### Symptom
DeepTools normalized BigWig files are missing from:
`$outdir/star/genome/deseq2/invariant_genes/deeptools_normalize/`

### Common Causes

#### 1. Space in Normalization Method Parameter ⚠️ CRITICAL

**Problem**: Including a space after the comma in the normalization method parameter.

```bash
# ❌ INCORRECT - Has space after comma
--normalization_method 'all_genes, invariant_genes'
```

**Solution**: Remove the space after the comma.

```bash
# ✅ CORRECT - No space
--normalization_method 'all_genes,invariant_genes'
```

**Why this matters**: The space prevents proper parsing of the comma-separated list. The pipeline treats it as a single malformed value instead of two separate normalization methods.

#### 2. DESeq2 Normalization Not Completing

The DeepTools normalization depends on DESeq2 QC outputs. Check that DESeq2 processes completed:

```bash
# Check for DESeq2 outputs
ls -lh $outdir/star/genome/deseq2/invariant_genes/

# Look for scaling factor files
find $outdir -name "*scaling_factor.txt" -path "*invariant_genes*"
```

#### 3. Insufficient Samples

Invariant genes normalization requires adequate sample size for robust estimation. Ensure you have:
- At least 3 samples per condition (recommended)
- Sufficient gene coverage for invariant gene detection

#### 4. Skip Flags Enabled

Verify that normalization is not being skipped:

```bash
# ✅ CORRECT - Do not skip
--skip_deseq2_qc False \
--skip_deeptools_norm False \
```

### Pipeline Logic Flow

Understanding the dependency chain helps troubleshoot:

```
STAR Alignment
    ↓
Feature Counts
    ↓
DESeq2 QC (all_genes) → Scaling Factors → DeepTools Norm (all_genes)
    ↓
DESeq2 QC (invariant_genes) → Scaling Factors → DeepTools Norm (invariant_genes)
    ↓
MultiQC
```

### Verification

After fixing, verify both normalization outputs exist:

```bash
# All genes normalization
ls -lh $outdir/star/genome/deseq2/all_genes/deeptools_normalize/

# Invariant genes normalization
ls -lh $outdir/star/genome/deseq2/invariant_genes/deeptools_normalize/
```

---

## Parameter Issues

### Valid Normalization Methods

The `--normalization_method` parameter accepts:

1. **Single method:**
   ```bash
   --normalization_method 'all_genes'
   # OR
   --normalization_method 'invariant_genes'
   ```

2. **Both methods (comma-separated, NO SPACE):**
   ```bash
   --normalization_method 'all_genes,invariant_genes'
   ```

### Parameters That Do NOT Exist

The following parameters do **not** exist in this pipeline and should not be used:

- ❌ `--publish_dir` - Use `--outdir` instead
- ❌ `--output` - Use `--outdir` instead
- ❌ `--results` - Use `--outdir` instead

### Required vs Optional Parameters

**Always Required:**
- `--input` - Path to samplesheet CSV
- `--outdir` - Output directory path
- `--fasta` or genome reference
- `--gtf` - Gene annotation file

**Commonly Used Optional:**
- `--aligner` - Default: 'star'
- `--quantification` - Default: 'genome'
- `--normalization_method` - Default: 'all_genes'
- `--skip_*` flags - Various QC steps

### Checking Available Parameters

To see all available parameters:

```bash
nextflow run pdichiaro/rnaseq --help
```

Or check the schema:

```bash
cat nextflow_schema.json | jq '.definitions'
```

---

## Output Structure

### Expected Directory Structure

After successful completion, your output should look like:

```
$outdir/
├── star/
│   ├── genome/
│   │   ├── deseq2/
│   │   │   ├── all_genes/
│   │   │   │   ├── deeptools_normalize/
│   │   │   │   │   ├── sample1.bw
│   │   │   │   │   ├── sample2.bw
│   │   │   │   │   └── ...
│   │   │   │   ├── *.counts.normalized.txt
│   │   │   │   ├── *.pca.plot.pdf
│   │   │   │   ├── *.sample.dists.plot.pdf
│   │   │   │   └── *.rds
│   │   │   └── invariant_genes/
│   │   │       ├── deeptools_normalize/
│   │   │       │   ├── sample1.bw
│   │   │       │   ├── sample2.bw
│   │   │       │   └── ...
│   │   │       ├── *.counts.normalized.txt
│   │   │       ├── *.pca.plot.pdf
│   │   │       ├── *.sample.dists.plot.pdf
│   │   │       └── *.rds
│   │   ├── qualimap/
│   │   └── ...
│   └── log/
├── multiqc/
│   └── star/
│       ├── multiqc_report.html
│       └── multiqc_data/
├── pipeline_info/
│   ├── execution_report*.html
│   ├── execution_timeline*.html
│   └── execution_trace*.txt
└── ...
```

### Missing Output Directories

If expected directories are missing:

1. **Check process completion:**
   ```bash
   grep "COMPLETED" $outdir/pipeline_info/execution_trace*.txt
   ```

2. **Check for errors:**
   ```bash
   grep "FAILED\|ERROR" $outdir/pipeline_info/execution_trace*.txt
   ```

3. **Review work directory:**
   ```bash
   # Find failed process directories
   find work/ -name ".exitcode" -exec grep -l "^[^0]" {} \;
   ```

---

## Example: Corrected Command

Here's a complete example with all common issues fixed:

```bash
#!/bin/bash

# Paths
sample_file=/path/to/samplesheet.csv
outdir=/path/to/output/
star_index=/path/to/star/index/
kallisto_index=/path/to/kallisto/index
GTF=/path/to/annotation.gtf
bed_file=/path/to/annotation.bed
genome_fasta=/path/to/genome.fa
transcriptome_fasta=/path/to/transcriptome.fa
conf_file=/path/to/local.config

# Run pipeline
nextflow run pdichiaro/rnaseq \
    -r main \
    --input $sample_file \
    --outdir $outdir \
    --fasta $genome_fasta \
    --transcript_fasta $transcriptome_fasta \
    --star_index $star_index \
    --kallisto_index $kallisto_index \
    --gene_bed $bed_file \
    --gtf $GTF \
    --gencode True \
    --gtf_extra_attributes 'gene_name' \
    --trimmer trimgalore \
    --extra_trimgalore_args '--quality 20 --stringency 3 --length 20' \
    --min_trimmed_reads 0 \
    --aligner star \
    --pseudo_aligner kallisto \
    --min_mapped_reads 0 \
    --extra_kallisto_quant_args '-b 50 --single-overhang' \
    --pseudo_aligner_kmer_size 31 \
    --with_umi False \
    --remove_ribo_rna False \
    --quantification genome \
    --normalization_method 'all_genes,invariant_genes' \
    --deseq2_vst True \
    --skip_linting True \
    --skip_gtf_filter True \
    --skip_fastqc True \
    --skip_rseqc False \
    --skip_trimming True \
    --skip_alignment False \
    --skip_pseudo_alignment True \
    --skip_markduplicates True \
    --skip_quantification_method False \
    --skip_bigwig False \
    --skip_deseq2_qc False \
    --skip_deeptools_norm False \
    --save_merged_fastq False \
    --save_trimmed False \
    --save_reference False \
    --save_align_intermeds False \
    --multiqc_title "My_RNA_Seq_Project" \
    -profile singularity \
    -c $conf_file \
    -process.echo \
    -resume
```

**Key Points:**
- ✅ No `--publish_dir` parameter
- ✅ Normalization method with no spaces: `'all_genes,invariant_genes'`
- ✅ All skip flags explicitly set
- ✅ Using `-resume` for efficient re-runs

---

## Getting Help

### Check Logs First

1. **Nextflow log:**
   ```bash
   tail -100 .nextflow.log
   ```

2. **Execution trace:**
   ```bash
   less $outdir/pipeline_info/execution_trace*.txt
   ```

3. **Process-specific logs:**
   ```bash
   # Find the work directory for a failed process
   grep "FAILED" $outdir/pipeline_info/execution_trace*.txt
   
   # Then check its log
   cat work/XX/XXXXXXX/.command.log
   ```

### Common Log Messages

**"Unknown parameter: publish_dir"**
- Remove the `--publish_dir` line from your command

**"Error parsing normalization_method"**
- Check for spaces in comma-separated values
- Use: `'all_genes,invariant_genes'` not `'all_genes, invariant_genes'`

**"Not enough samples for invariant genes normalization"**
- Ensure you have at least 3 samples per condition
- Consider using only `'all_genes'` normalization

**"MultiQC found no samples"**
- Check that upstream QC processes completed
- Verify files exist in expected locations

---

## Performance Tips

### Using -resume Effectively

Always use `-resume` when re-running after fixing issues:

```bash
nextflow run pdichiaro/rnaseq -resume ...
```

This will:
- ✅ Skip completed processes
- ✅ Only re-run failed or modified steps
- ✅ Save significant computation time

### Work Directory Management

The work directory can grow large. After successful completion:

```bash
# Remove work directory
rm -rf work/

# Or keep only specific process outputs
nextflow clean -n  # Preview what will be removed
nextflow clean -f  # Actually remove
```

### Resource Optimization

If processes are failing due to resources:

1. Check your config file for resource settings
2. Adjust memory/CPU in `local.config`:
   ```groovy
   process {
       withName: 'STAR_ALIGN' {
           cpus = 16
           memory = '64.GB'
       }
   }
   ```

---

## Version Information

This troubleshooting guide is current as of:
- Pipeline version: main branch
- Nextflow version: 25.04.7
- Last updated: 2025-11-27

For the latest updates, check the repository:
https://github.com/pdichiaro/rnaseq

---

## Quick Reference Checklist

Before running the pipeline, verify:

- [ ] No `--publish_dir` parameter in command
- [ ] Normalization method has no spaces: `'all_genes,invariant_genes'`
- [ ] All input files exist and are readable
- [ ] Output directory has write permissions
- [ ] Sufficient disk space available
- [ ] Required reference files (GTF, FASTA) are correct
- [ ] Samplesheet CSV is properly formatted
- [ ] Config file paths are correct
- [ ] Using correct Nextflow version (25.04.7+)

After running, verify outputs:

- [ ] MultiQC report exists: `$outdir/multiqc/star/multiqc_report.html`
- [ ] All genes normalization: `$outdir/star/genome/deseq2/all_genes/deeptools_normalize/`
- [ ] Invariant genes normalization: `$outdir/star/genome/deseq2/invariant_genes/deeptools_normalize/`
- [ ] Pipeline info reports: `$outdir/pipeline_info/`
- [ ] No failed processes in execution trace

---

*For additional support, please open an issue on the GitHub repository.*
