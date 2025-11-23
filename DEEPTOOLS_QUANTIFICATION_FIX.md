# DeepTools Quantification-Aware Publishing Fix

## Issue
DeepTools bigwig files were always published to `${aligner}/rsem/deeptools/` regardless of the quantification method used. When running with `--quantification genome`, files should be published to `${aligner}/genome/deeptools/` instead.

## Root Cause
The publishDir configuration in `nextflow.config` was hardcoded to publish to the `rsem` subdirectory, not taking into account the actual quantification method being used.

## Solution

### 1. Updated Module Output (modules/local/deeptools_bw_norm/main.nf)

**Changed from:**
```groovy
output:
    path "*.unstranded.norm.bw" , emit: unstranded_bw
    path "*.fwd.norm.bw"        , optional:true, emit: fw_bw
    path "*.rev.norm.bw"        , optional:true, emit: rev_bw
    path "versions.yml"         , emit: versions
```

**Changed to:**
```groovy
output:
    tuple val(meta), path("*.unstranded.norm.bw"), emit: unstranded_bw
    tuple val(meta), path("*.fwd.norm.bw")       , optional:true, emit: fw_bw
    tuple val(meta), path("*.rev.norm.bw")       , optional:true, emit: rev_bw
    path "versions.yml"                          , emit: versions
```

**Rationale:** Now the module emits metadata along with files, allowing the publishDir to access the quantification method stored in `meta.quantification`.

### 2. Updated PublishDir Configuration (nextflow.config)

**Changed from:**
```groovy
withName: 'DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
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
            return "${params.outdir}/${aligner}/rsem/deeptools/invariant_genes"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

**Changed to:**
```groovy
withName: 'DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
    publishDir = [
        path: { meta ->
            def aligner = 'unknown'
            if (!params.skip_alignment && params.aligner == 'star') {
                aligner = 'star'
            } else if (!params.skip_alignment && params.aligner == 'hisat2') {
                aligner = 'hisat2'
            } else if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'kallisto') {
                aligner = 'kallisto'
            }
            
            // Get quantification method from meta (set by workflow)
            def quant_method = meta.quantification ?: 'rsem'
            
            return "${params.outdir}/${aligner}/${quant_method}/deeptools/invariant_genes"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

**Same change applied to:** `DEEPTOOLS_BIGWIG_NORM_ALL_GENES`

**Rationale:** The publishDir closure now:
1. Receives the `meta` map from the tuple output
2. Extracts the `quantification` field (set by the workflow in lines 1126 and 1149 of main.nf)
3. Uses it to construct the correct output path dynamically

### 3. How It Works

The workflow already sets the quantification method in the metadata:

```groovy
// From workflows/rnaseq/main.nf lines ~1126 and ~1149
ch_combined_input_all_genes = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample_all_genes)
    .map { meta, bam, bai, sample_id, scaling, quant_method -> 
        if (meta.id == sample_id) {
            def new_meta = meta.clone()
            new_meta.quantification = quant_method  // ← Sets quantification method
            [new_meta, bam, bai, scaling]
        } else {
            null
        }
    }
```

The quantification method is detected from the scaling factor file path:
- Files in `/rsem/` → `quant_method = 'rsem'`
- Files in `/genome/` → `quant_method = 'genome'`
- Files in `/salmon/` → `quant_method = 'salmon'`

## Expected Directory Structure

### With `--quantification genome`:
```
results/
└── star/
    └── genome/
        ├── sample1/
        │   └── *_combined_counts.txt
        ├── sample2/
        │   └── *_combined_counts.txt
        ├── deseq2/
        │   ├── all_genes/
        │   │   └── normalization outputs
        │   └── invariant_genes/
        │       └── normalization outputs
        └── deeptools/                    ← NEW: deeptools in genome folder
            ├── all_genes/
            │   ├── sample1.unstranded.norm.bw
            │   ├── sample1.fwd.norm.bw
            │   └── sample1.rev.norm.bw
            └── invariant_genes/
                ├── sample1.unstranded.norm.bw
                ├── sample1.fwd.norm.bw
                └── sample1.rev.norm.bw
```

### With `--quantification rsem`:
```
results/
└── star/
    └── rsem/
        ├── sample1/
        │   └── *.results
        ├── sample2/
        │   └── *.results
        ├── deseq2/
        │   ├── all_genes/
        │   └── invariant_genes/
        └── deeptools/                    ← deeptools in rsem folder
            ├── all_genes/
            └── invariant_genes/
```

### With `--quantification genome,rsem`:
```
results/
└── star/
    ├── genome/
    │   ├── sample1/
    │   ├── deseq2/
    │   └── deeptools/                   ← genome deeptools
    │       ├── all_genes/
    │       └── invariant_genes/
    └── rsem/
        ├── sample1/
        ├── deseq2/
        └── deeptools/                   ← rsem deeptools
            ├── all_genes/
            └── invariant_genes/
```

## Testing

Test with:
```bash
nextflow run workflows/rnaseq/main.nf \
  --aligner star \
  --quantification genome \
  --skip_deeptools_norm false \
  ...other params...
```

Verify that deeptools outputs appear in `results/star/genome/deeptools/` instead of `results/star/rsem/deeptools/`.

## Files Modified
1. `modules/local/deeptools_bw_norm/main.nf` - Updated outputs to include meta
2. `nextflow.config` - Updated publishDir for both DEEPTOOLS_BIGWIG_NORM processes

## Related Issues
- This fix complements the earlier deeptools configuration fix documented in `SUMMARY_DEEPTOOLS_FIX.md`
- Ensures consistent directory structure across all quantification methods
