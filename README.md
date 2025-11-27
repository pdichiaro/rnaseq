# pdichiaro/rnaseq

RNA-seq analysis pipeline built with Nextflow DSL2.

## Quick Links

- 🚀 **[Quick Start Fixes](QUICK_START_FIXES.md)** - Fix the 2 most common issues in 2 minutes
- 🔧 **[Troubleshooting Guide](TROUBLESHOOTING_COMMON_ISSUES.md)** - Comprehensive troubleshooting documentation
- 📁 **[Folder Structure](FOLDER_STRUCTURE.md)** - Understanding the output structure

## Overview

This pipeline performs comprehensive RNA-seq analysis including:

- Quality control (FastQC, RSeQC, Qualimap)
- Read trimming (TrimGalore or Fastp)
- Alignment (STAR or HISAT2)
- Pseudo-alignment (Kallisto)
- Quantification (featureCounts, RSEM, Kallisto)
- Normalization (DESeq2 with all_genes and/or invariant_genes methods)
- Normalized BigWig generation (DeepTools)
- Comprehensive reporting (MultiQC)

## Quick Start

### Installation

Requires:
- Nextflow ≥ 25.04.7
- Singularity or Docker
- Reference genome files (FASTA, GTF, indices)

### Basic Usage

```bash
nextflow run pdichiaro/rnaseq \
    -r main \
    --input samplesheet.csv \
    --outdir results/ \
    --fasta genome.fa \
    --gtf annotation.gtf \
    --star_index star_index/ \
    --gene_bed genes.bed \
    -profile singularity
```

### Important: Common Parameter Issues ⚠️

**DO NOT use these parameters** (they don't exist):
- ❌ `--publish_dir`
- ❌ `--output`
- ❌ `--results`

**Use only:** ✅ `--outdir`

**Comma-separated parameters** must have NO SPACES:
```bash
# ✅ CORRECT
--normalization_method 'all_genes,invariant_genes'

# ❌ WRONG
--normalization_method 'all_genes, invariant_genes'
```

See [QUICK_START_FIXES.md](QUICK_START_FIXES.md) for immediate solutions to common issues.

## Input Format

### Samplesheet CSV

Required columns: `sample`, `fastq_1`, `fastq_2`, `strandedness`

```csv
sample,fastq_1,fastq_2,strandedness
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,auto
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,auto
```

## Key Parameters

### Required
- `--input` - Path to samplesheet CSV
- `--outdir` - Output directory
- `--fasta` - Genome FASTA file
- `--gtf` - Gene annotation GTF file

### Alignment & Quantification
- `--aligner` - Alignment tool: `star` (default), `hisat2`
- `--pseudo_aligner` - Pseudo-alignment: `kallisto`
- `--quantification` - Quantification method: `genome` (default), `rsem`, `genome,rsem`
- `--normalization_method` - Normalization: `all_genes` (default), `invariant_genes`, `all_genes,invariant_genes`

### Reference Files
- `--star_index` - STAR index directory
- `--kallisto_index` - Kallisto index file
- `--gene_bed` - Gene BED12 file
- `--transcript_fasta` - Transcript FASTA file

### Skip Options
- `--skip_fastqc` - Skip FastQC (default: false)
- `--skip_trimming` - Skip read trimming (default: false)
- `--skip_alignment` - Skip alignment (default: false)
- `--skip_pseudo_alignment` - Skip pseudo-alignment (default: false)
- `--skip_deseq2_qc` - Skip DESeq2 QC (default: false)
- `--skip_deeptools_norm` - Skip DeepTools normalization (default: false)
- `--skip_multiqc` - Skip MultiQC (default: false)

## Output Structure

```
results/
├── star/
│   ├── genome/
│   │   ├── deseq2/
│   │   │   ├── all_genes/
│   │   │   │   ├── deeptools_normalize/     # Normalized BigWig files
│   │   │   │   ├── *.counts.normalized.txt  # Normalized counts
│   │   │   │   ├── *.pca.plot.pdf           # PCA plots
│   │   │   │   └── *.sample.dists.plot.pdf  # Sample distance plots
│   │   │   └── invariant_genes/
│   │   │       ├── deeptools_normalize/     # Normalized BigWig files
│   │   │       ├── *.counts.normalized.txt
│   │   │       ├── *.pca.plot.pdf
│   │   │       └── *.sample.dists.plot.pdf
│   │   ├── qualimap/
│   │   └── ...
│   └── log/
├── kallisto/
│   └── deseq2/
├── multiqc/
│   └── star/
│       ├── multiqc_report.html             # Main QC report
│       └── multiqc_data/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## Normalization Methods

### All Genes Normalization
Standard DESeq2 normalization using all genes. Suitable for most RNA-seq datasets.

### Invariant Genes Normalization
Alternative normalization using only invariant (housekeeping) genes. Useful when:
- Global transcriptional changes are expected
- Treatment affects many genes
- Standard normalization assumptions may not hold

**Requirements:**
- At least 3 samples per condition (recommended)
- Sufficient gene coverage
- Adequate sequencing depth

### Using Both Methods

```bash
--normalization_method 'all_genes,invariant_genes'
```

This generates both normalizations for comparison.

## Configuration

### Custom Config File

Create a `local.config` file to customize resource allocation:

```groovy
process {
    executor = 'local'
    
    withName: 'STAR_ALIGN' {
        cpus = 16
        memory = '64.GB'
    }
    
    withName: 'DESEQ2_*' {
        cpus = 4
        memory = '32.GB'
    }
}
```

Use with: `-c local.config`

### Profiles

Available profiles:
- `singularity` - Use Singularity containers (recommended for HPC)
- `docker` - Use Docker containers
- `conda` - Use Conda environments

## Advanced Usage

### Complete Example with All Options

```bash
#!/bin/bash

# Paths
sample_file=samplesheet.csv
outdir=results/
star_index=references/star_index/
kallisto_index=references/kallisto_index/kallisto
GTF=references/annotation.gtf
bed_file=references/genes.bed
genome_fasta=references/genome.fa
transcriptome_fasta=references/transcriptome.fa
conf_file=local.config

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
    --skip_fastqc False \
    --skip_trimming False \
    --skip_alignment False \
    --skip_pseudo_alignment True \
    --skip_deseq2_qc False \
    --skip_deeptools_norm False \
    --multiqc_title "My_RNA_Seq_Project" \
    -profile singularity \
    -c $conf_file \
    -resume
```

### Resume Functionality

Always use `-resume` to continue from cached steps:

```bash
nextflow run pdichiaro/rnaseq -resume ...
```

This saves time by skipping successfully completed processes.

## Troubleshooting

### Most Common Issues

1. **No MultiQC report**
   - Remove `--publish_dir` parameter (doesn't exist)
   - See: [QUICK_START_FIXES.md](QUICK_START_FIXES.md)

2. **Missing invariant genes normalization**
   - Remove space in: `'all_genes,invariant_genes'`
   - See: [QUICK_START_FIXES.md](QUICK_START_FIXES.md)

3. **Pipeline fails**
   - Check: `.nextflow.log`
   - Check: `results/pipeline_info/execution_trace.txt`
   - See: [TROUBLESHOOTING_COMMON_ISSUES.md](TROUBLESHOOTING_COMMON_ISSUES.md)

### Getting Help

1. Check [TROUBLESHOOTING_COMMON_ISSUES.md](TROUBLESHOOTING_COMMON_ISSUES.md)
2. Review `.nextflow.log` for error messages
3. Check execution trace in `pipeline_info/`
4. Open an issue on GitHub with:
   - Command used
   - Error message
   - Relevant log excerpts

## Performance Tips

### Resource Optimization

- Use `-resume` to avoid re-computing
- Adjust process resources in config file
- Use appropriate profile (singularity/docker)
- Consider work directory cleanup after success

### Disk Space Management

```bash
# Check work directory size
du -sh work/

# Clean after successful completion
rm -rf work/

# Or use Nextflow clean
nextflow clean -f
```

## Citation

If you use this pipeline, please cite:

- Nextflow: https://www.nextflow.io/
- nf-core tools and modules used
- Individual tools (STAR, DESeq2, MultiQC, etc.)

## Version Information

- Pipeline version: main branch
- Nextflow version: ≥25.04.7
- DSL: 2

## License

See LICENSE file in repository.

## Credits

This pipeline uses components from:
- nf-core modules and subworkflows
- STAR aligner
- DESeq2
- MultiQC
- And many other open-source tools

## Contact

For issues and questions:
- GitHub Issues: https://github.com/pdichiaro/rnaseq/issues
- Email: [Maintainer contact]

---

**Quick Links:**
- 🚀 [Quick Start Fixes](QUICK_START_FIXES.md)
- 🔧 [Troubleshooting Guide](TROUBLESHOOTING_COMMON_ISSUES.md)
- 📁 [Folder Structure](FOLDER_STRUCTURE.md)
