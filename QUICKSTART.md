# 🚀 RNA-seq Pipeline Quick Start Guide

## ⚡ Minimum Command

```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile docker
```

## 🔴 MANDATORY Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Sample sheet CSV file | `samples.csv` |
| `--outdir` | Output directory | `results/` |
| `-profile` | Container system | `docker` |

## 🔶 Reference Genome (Choose ONE)

### Option A: iGenomes (Recommended)
```bash
--genome GRCh38    # Human (latest)
--genome GRCm39    # Mouse (latest)  
--genome <ID>      # Other organisms
```

### Option B: Custom Files
```bash
--fasta genome.fa  # Genome sequence
--gtf genes.gtf    # Gene annotations
```

## 📝 Sample Sheet Format

Create `samplesheet.csv`:
```csv
sample,fastq_1,fastq_2,strandedness
SAMPLE1,reads_1.fq.gz,reads_2.fq.gz,auto
SAMPLE2,reads_1.fq.gz,reads_2.fq.gz,auto
```

## 🎯 Common Scenarios

### Standard Paired-end Analysis
```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile docker
```

### UMI-based Analysis  
```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --with_umi \
    --umitools_bc_pattern 'NNNNNNNN' \
    -profile docker
```

### Fast Analysis (Skip Some QC)
```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --skip_dupradar \
    --skip_preseq \
    --skip_rseqc \
    -profile docker
```

### Custom Genome + Salmon Only
```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --fasta genome.fasta \
    --gtf annotation.gtf \
    --pseudo_aligner salmon \
    --skip_alignment \
    -profile docker
```

## ⚠️ Quick Troubleshooting

| Issue | Solution |
|-------|----------|
| "Profile not specified" | Add `-profile docker` |
| "Input file not found" | Check `--input` path |
| "Genome not found" | Use valid `--genome` ID or provide `--fasta`+`--gtf` |
| "UMI pattern required" | Add `--umitools_bc_pattern` when using `--with_umi` |
| "Permission denied" | Check write access to `--outdir` |

## 🧬 Popular Genome IDs

| Organism | ID | Version |
|----------|----|---------| 
| Human | `GRCh38` | hg38 |
| Human | `GRCh37` | hg19 |
| Mouse | `GRCm39` | mm39 |
| Mouse | `GRCm38` | mm10 |
| Fly | `dm6` | Release 6 |
| Worm | `ce11` | WBcel235 |

## 📊 Expected Outputs

```
results/
├── multiqc/                 # Quality reports
├── star_salmon/             # Quantification  
├── fastqc/                  # Read QC
├── trimgalore/              # Trimmed reads
└── pipeline_info/           # Run metadata
```

## 🚀 Next Steps

1. Check `results/multiqc/multiqc_report.html` for QC overview
2. Find gene counts in `results/star_salmon/`
3. Review logs in `results/pipeline_info/`
4. Customize analysis with additional parameters (see `docs/usage.md`)