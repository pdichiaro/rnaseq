# DeepTools Workflow Config Fix

## Root Cause Analysis

The issue you were seeing was caused by **configuration override priority** in Nextflow:

1. **Main `nextflow.config`** had the correct logic (commit `d7c3bef`)
2. **BUT** `workflows/rnaseq/nextflow.config` is loaded AFTER and contains process-specific configurations
3. The workflow config had DeepTools publishDir settings **only for STAR and HISAT2**:
   ```groovy
   if ((['star', 'hisat2'].contains(params.aligner)) && !params.skip_deeptools_norm) {
       // DeepTools config here
   }
   ```
4. This conditional **excluded pseudo-aligners** completely!
5. So for Kallisto runs, it would fall back to the main config, which we had fixed
6. BUT the workflow config was still being loaded and might have caused conflicts

## The Fix

Added a **separate conditional block** in `workflows/rnaseq/nextflow.config` specifically for pseudo-aligners:

```groovy
//
// Deeptools BigWig normalization options for pseudo-aligners (Kallisto)
//

if (!params.skip_pseudo_alignment && params.pseudo_aligner && !params.skip_deeptools_norm) {
    process {
        withName: '.*:DEEPTOOLS_BIGWIG_NORM_ALL_GENES' {
            publishDir = [
                path: { 
                    // Pseudo-aligners don't have quantification methods like rsem/genome
                    "${params.outdir}/${params.pseudo_aligner}/deeptools/all_genes" 
                },
                mode: params.publish_dir_mode,
                pattern: "*.norm.bw"
            ]
        }
        
        withName: '.*:DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
            publishDir = [
                path: { 
                    // Pseudo-aligners don't have quantification methods like rsem/genome
                    "${params.outdir}/${params.pseudo_aligner}/deeptools/invariant_genes" 
                },
                mode: params.publish_dir_mode,
                pattern: "*.norm.bw"
            ]
        }
    }
}
```

## Configuration Loading Order

Nextflow loads configurations in this order (later overrides earlier):

1. **Main `nextflow.config`** ← Our first fix
2. **Profile configs** (test.config, docker.config, etc.)
3. **Workflow configs** via `includeConfig` ← This was overriding!
4. **Command-line parameters** ← Highest priority

## Files Modified

- ✅ `nextflow.config` (commit `d7c3bef`) - Added logic to handle pseudo-aligners vs aligners
- ✅ `workflows/rnaseq/nextflow.config` (commit `f0e741b`) - Added pseudo-aligner specific block

## Expected Output Structure

### For Kallisto (Pseudo-aligner)
```
results/
└── kallisto/
    ├── deseq2/
    │   ├── all_genes/
    │   └── invariant_genes/
    └── deeptools/
        ├── all_genes/
        │   └── *.norm.bw
        └── invariant_genes/
            └── *.norm.bw
```

### For STAR/HISAT2 (Aligners)
```
results/
└── star/
    ├── rsem/
    │   ├── deseq2/
    │   └── deeptools/
    │       ├── all_genes/
    │       │   └── *.norm.bw
    │       └── invariant_genes/
    │           └── *.norm.bw
    └── genome/
        ├── deseq2/
        └── deeptools/
```

## Lesson Learned

When dealing with Nextflow configuration issues:
1. **Check ALL config files** - not just the main one!
2. **Look for `includeConfig` statements** - these can override main settings
3. **Understand conditional blocks** - they might exclude certain parameter combinations
4. **Check workflow-specific configs** - they load last and have high priority

## Testing Recommendation

To verify the fix works:
```bash
# For Kallisto
nextflow run . --pseudo_aligner kallisto --skip_alignment --outdir results

# Expected structure:
# results/kallisto/deeptools/all_genes/*.norm.bw
# results/kallisto/deeptools/invariant_genes/*.norm.bw
```
