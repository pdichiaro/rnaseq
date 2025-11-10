# pdichiaro/rnaseq: Usage

## Quick Start

### Standard RNA-seq Analysis
```bash
nextflow run pdichiaro/rnaseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile singularity
```

| Parameter | Status | Description | Default |
|-----------|--------|-------------|---------|
| `--input` | ** MANDATORY** | Path to samplesheet CSV file | `null` |
| `--outdir` | ** MANDATORY** | Output directory path | `null` |
| `-profile` | ** MANDATORY** | Execution environment (docker/singularity/conda) | None |


## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto
```


## Complete Parameter Reference

### Mandatory Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--input` | path | Path to samplesheet CSV | `samplesheet.csv` |
| `--outdir` | path | Output directory | `results/` |

###  Conditionally Mandatory Parameters

| Parameter | Condition | Type | Description |
|-----------|-----------|------|-------------|
| `--genome` | If no custom refs | string | iGenomes reference ID |
| `--fasta` | If no --genome | path | Genome FASTA file |
| `--gtf` | If no --genome | path | Gene annotation GTF |
| `--umitools_bc_pattern` | If --with_umi | string | UMI barcode pattern |
| `--bbsplit_fasta_list` | If not --skip_bbsplit | path | BBSplit reference list |
| `--kraken_db` | If contaminant screening | path | Kraken2 database |

###  Key Optional Parameters

| Category | Parameter | Default | Description |
|----------|-----------|---------|-------------|
| **Alignment** | `--aligner` | `star` | Alignment method (star/hisat2) |
| | `--quantification` | `salmon` | Quantification method |
| | `--pseudo_aligner` | `null` | Pseudo-aligner (salmon/kallisto) |
| **Quality** | `--trimmer` | `trimgalore` | Trimming tool |
| | `--min_trimmed_reads` | `10000` | Min reads after trimming |
| | `--skip_qc` | `false` | Skip all QC steps |
| **Analysis** | `--with_umi` | `false` | Enable UMI analysis |
| | `--stranded_threshold` | `0.8` | Strandedness detection |
| | `--normalization_method` | `all_genes` | DESeq2 normalization |

### Skip Options (All default to false)

- `--skip_trimming` - Skip read trimming
- `--skip_alignment` - Skip alignment steps  
- `--skip_pseudo_alignment` - Skip pseudo-alignment
- `--skip_markduplicates` - Skip duplicate marking
- `--skip_stringtie` - Skip StringTie assembly
- `--skip_fastqc` - Skip FastQC reports
- `--skip_multiqc` - Skip MultiQC report
- `--skip_deseq2_qc` - Skip DESeq2 QC plots
- `--skip_biotype_qc` - Skip biotype QC
- `--skip_dupradar` - Skip dupRadar analysis
- `--skip_preseq` - Skip Preseq analysis
- `--skip_qualimap` - Skip Qualimap QC
- `--skip_rseqc` - Skip RSeQC analysis

### Save Options (All default to false)

- `--save_reference` - Save generated indices
- `--save_trimmed` - Save trimmed reads
- `--save_unaligned` - Save unaligned reads
- `--save_align_intermeds` - Save intermediate BAMs
- `--save_merged_fastq` - Save merged FASTQ files

