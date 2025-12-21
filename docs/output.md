# pdichiaro/rnaseq: Output

## Directory Structure by Aligner and Quantification Method

This document shows the complete output directory structure for each aligner, quantification, and normalization method combination.

## Common Preprocessing Outputs (all methods)
```
results/
├── fastq/                        # Merged FastQ files (if --save_merged_fastq)
├── fq_lint/                      # FASTQ linting reports  
├── fastqc/                       # FastQC quality reports
│   ├── <SAMPLE>_fastqc.html
│   └── <SAMPLE>_fastqc.zip
├── umitools/                     # UMI extraction outputs (if --with_umi)
├── trimgalore/                   # TrimGalore outputs (if trimmer=trimgalore)
│   ├── <SAMPLE>_val_*.fq.gz     # Trimmed reads
│   ├── logs/
│   └── fastqc/                   # Post-trim FastQC
├── fastp/                        # fastp outputs (if trimmer=fastp)
│   ├── <SAMPLE>_*.fastp.fastq.gz # Trimmed reads  
│   ├── logs/
│   └── fastqc/                   # Post-trim FastQC
├── bbsplit/                      # BBSplit contamination filtering
├── sortmerna/                    # rRNA removal outputs
└── samplesheets/                 # Modified samplesheets
```

## STAR (--aligner star)
```
star/                          # Aligner directory
├── [common BAM files, logs, general QC] # Shared across all quantification methods
├── rsem/                      # Quantification subdirectory
│   ├── [rsem quantification files]
│   ├── deseq2/                # Analysis type subdirectory
│   │   ├── all_genes/         # Normalization method subdirectory
│   │   │   ├── *.plots.pdf    # DESeq2 results for all_genes method
│   │   │   ├── *.dds.RData
│   │   │   └── size_factors/
│   │   └── invariant_genes/   # Normalization method subdirectory  
│   │       └── [Same DESeq2 structure]
│   └── bigwig_norm/           # Analysis type subdirectory
│       ├── all_genes/         # Normalization method subdirectory
│       │   └── *.norm.bw      # Normalized BigWigs using all_genes factors
│       └── invariant_genes/   # Normalization method subdirectory
│           └── *.norm.bw      # Normalized BigWigs using invariant_genes factors
├── salmon/                    # Quantification subdirectory
│   ├── [salmon quantification files]
│   ├── deseq2/                # Analysis type subdirectory
│   │   ├── all_genes/         # Normalization method subdirectory
│   │   └── invariant_genes/   # Normalization method subdirectory
│   └── bigwig_norm/           # Analysis type subdirectory
│       ├── all_genes/         # Normalization method subdirectory
│       └── invariant_genes/   # Normalization method subdirectory
└── genome/                    # Quantification subdirectory
    ├── [genome quantification files] 
    ├── deseq2/                # Analysis type subdirectory
    │   ├── all_genes/         # Normalization method subdirectory
    │   └── invariant_genes/   # Normalization method subdirectory
    └── bigwig_norm/           # Analysis type subdirectory
        ├── all_genes/         # Normalization method subdirectory  
        └── invariant_genes/   # Normalization method subdirectory
```

## HISAT2 (--aligner hisat2)

**Note**: HISAT2 only supports `--quantification genome` using GENOME_COUNT module.


## Pseudoalignment (--pseudo_aligner kallisto)

```
kallisto /
├── [common BAM files, logs, general QC] # Shared across all quantification methods
├── tximport/                    # tximport processing outputs
    ├── [quantification files] 
    ├── deseq2/                # Analysis type subdirectory
    │   ├── all_genes/         # Normalization method subdirectory
    │   └── invariant_genes/   # Normalization method subdirectory
    └── bigwig_norm/           # Analysis type subdirectory
        ├── all_genes/         # Normalization method subdirectory  
        └── invariant_genes/   # Normalization method subdirectory
```

## Pseudoalignment (--pseudo_aligner salmon)

**Note**: SALMON does not generate BAM files.


## Reference Genome Files (if --save_reference)
```
genome/
├── *.fa                         # Reference genome FASTA
├── *.gtf                        # Annotation GTF file
├── *.gff                        # Annotation GFF file (if provided)  
├── *.bed                        # Gene BED file
├── *.tsv                        # Various annotation files
└── index/                       # Pre-built indices
    ├── star/                    # STAR index
    ├── hisat2/                  # HISAT2 index
    ├── rsem/                    # RSEM index  
    ├── salmon/                  # Salmon index
    └── kallisto/                # Kallisto index
```

## Pipeline Information
```
pipeline_info/
├── nf_core_rnaseq_software_mqc_versions.yml # Software versions
├── execution_report.html        # Nextflow execution report
├── execution_timeline.html      # Execution timeline  
├── execution_trace.txt          # Execution trace
└── pipeline_dag.svg             # Pipeline DAG (if -with-dag)
```

## Key Notes
- **NEW STRUCTURE**: STAR uses `star/` base directory with `rsem/`, `salmon/`, `genome/` subdirectories
- `<SAMPLE>` represents individual sample names from your samplesheet
- Many outputs are optional and depend on pipeline parameters (e.g., --save_align_intermeds, --skip_qc)
- All alignment methods include the same QC outputs (RSeQC, Qualimap, etc.) unless skipped
- DESeq2 normalization is available for all quantification methods
- MultiQC aggregates all QC results into a single HTML report
