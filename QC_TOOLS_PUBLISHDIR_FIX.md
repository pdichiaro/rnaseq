# QC Tools PublishDir Organization Fix

## Issue
QC tools (Qualimap, DupRadar, Preseq) were publishing to top-level folders instead of being organized under the aligner folder, creating an inconsistent directory structure.

## Problem
The default publishDir configuration uses:
```groovy
path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" }
```

For QC tools, this resulted in:
```
results/
├── qualimap/       ← Top level, not organized
├── dupradar/       ← Top level, not organized
├── preseq/         ← Top level, not organized
└── star/
    ├── genome/
    └── rsem/
```

## Expected Structure
QC tools should be organized under the aligner folder since they analyze alignment-specific BAM files:

```
results/
└── star/           (or hisat2/kallisto)
    ├── genome/
    │   ├── deseq2/
    │   └── deeptools/
    ├── rsem/
    │   ├── deseq2/
    │   └── deeptools/
    ├── qualimap/    ← Organized under aligner
    │   ├── sample1/
    │   │   ├── qualimapReport.html
    │   │   └── rnaseq_qc_results.txt
    │   └── sample2/
    ├── dupradar/    ← Organized under aligner
    │   ├── box_plot/
    │   ├── gene_data/
    │   ├── histogram/
    │   └── scatter_plot/
    └── preseq/      ← Organized under aligner
        ├── sample1.lc_extrap.txt
        └── sample2.lc_extrap.txt
```

## Solution

Added explicit publishDir configurations for QC tools in `nextflow.config`:

### 1. Qualimap Configuration
```groovy
withName: 'QUALIMAP_RNASEQ' {
    publishDir = [
        path: { 
            def aligner = 'unknown'
            if (!params.skip_alignment && params.aligner == 'star') {
                aligner = 'star'
            } else if (!params.skip_alignment && params.aligner == 'hisat2') {
                aligner = 'hisat2'
            } else if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'kallisto') {
                aligner = 'kallisto'
            }
            return "${params.outdir}/${aligner}/qualimap"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

### 2. DupRadar Configuration
```groovy
withName: 'DUPRADAR' {
    publishDir = [
        path: { 
            def aligner = 'unknown'
            if (!params.skip_alignment && params.aligner == 'star') {
                aligner = 'star'
            } else if (!params.skip_alignment && params.aligner == 'hisat2') {
                aligner = 'hisat2'
            } else if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'kallisto') {
                aligner = 'kallisto'
            }
            return "${params.outdir}/${aligner}/dupradar"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

### 3. Preseq Configuration
```groovy
withName: 'PRESEQ_LCEXTRAP' {
    publishDir = [
        path: { 
            def aligner = 'unknown'
            if (!params.skip_alignment && params.aligner == 'star') {
                aligner = 'star'
            } else if (!params.skip_alignment && params.aligner == 'hisat2') {
                aligner = 'hisat2'
            } else if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'kallisto') {
                aligner = 'kallisto'
            }
            return "${params.outdir}/${aligner}/preseq"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

## Rationale

1. **Logical Grouping**: QC tools analyze BAM files from specific aligners, so they should be grouped with the aligner output

2. **Consistency**: Matches the structure already used for:
   - DESeq2 normalization (`star/genome/deseq2/`, `star/rsem/deseq2/`)
   - DeepTools (`star/genome/deeptools/`, `star/rsem/deeptools/`)

3. **Clarity**: Users can find all STAR-related outputs under `star/`, all HISAT2 outputs under `hisat2/`, etc.

4. **Scalability**: When running multiple aligners, outputs remain organized

## Benefits

### Before:
```
results/
├── qualimap/       ← Which aligner?
├── dupradar/       ← Which aligner?
├── preseq/         ← Which aligner?
├── star/
│   └── ...
└── hisat2/
    └── ...
```

### After:
```
results/
├── star/
│   ├── qualimap/   ← Clearly STAR QC
│   ├── dupradar/   ← Clearly STAR QC
│   └── preseq/     ← Clearly STAR QC
└── hisat2/
    ├── qualimap/   ← Clearly HISAT2 QC
    ├── dupradar/   ← Clearly HISAT2 QC
    └── preseq/     ← Clearly HISAT2 QC
```

## About the QC Tools

### Qualimap
- **Purpose**: RNA-seq quality control for alignment data
- **Analyzes**: BAM files
- **Outputs**: 
  - `qualimapReport.html` - Interactive HTML report
  - `rnaseq_qc_results.txt` - Text summary
  - Coverage plots, genomic origin plots
- **Metrics**: Gene coverage, read distribution, junction analysis

### DupRadar
- **Purpose**: Detect technical duplication in RNA-seq
- **Analyzes**: BAM files + GTF annotation
- **Outputs**:
  - Scatter plots of duplication vs expression
  - Box plots of duplicate rates
  - Per-gene duplication metrics
- **Use**: Identify experiments with high technical duplication

### Preseq
- **Purpose**: Estimate library complexity and predict future sequencing depth
- **Analyzes**: BAM files
- **Outputs**: 
  - `.lc_extrap.txt` - Complexity curve data
- **Use**: Determine if deeper sequencing would yield more unique reads

## Migration Note

If you have existing results with QC outputs at the top level:
1. Old outputs remain in place (no automatic migration)
2. New runs will use the organized structure
3. MultiQC will still collect data from both locations

## Files Modified
- `nextflow.config` - Added publishDir configurations for QUALIMAP_RNASEQ, DUPRADAR, and PRESEQ_LCEXTRAP

## Testing
Run with STAR aligner:
```bash
nextflow run workflows/rnaseq/main.nf \
  --aligner star \
  --quantification genome \
  ...
```

Verify outputs appear in:
- `results/star/qualimap/`
- `results/star/dupradar/`
- `results/star/preseq/`

## Related Issues
- Complements deeptools publishDir organization (see DEEPTOOLS_QUANTIFICATION_FIX.md)
- Provides consistent structure across all QC tools
- Resolves "empty qualimap folder" confusion by placing it in expected location
