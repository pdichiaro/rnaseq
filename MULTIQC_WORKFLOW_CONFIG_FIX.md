# MultiQC Workflow Config Fix

## Issue Found

The workflow configuration file (`workflows/rnaseq/nextflow.config`) was using **dynamic parameter references** instead of hardcoded paths for MultiQC publishDir settings.

### Before (WRONG):
```groovy
withName: '.*:MULTIQC_STAR' {
    publishDir = [
        path: { "${params.outdir}/multiqc/${params.aligner}" }  // ❌ Dynamic
    ]
}

withName: '.*:MULTIQC_HISAT2' {
    publishDir = [
        path: { "${params.outdir}/multiqc/${params.aligner}" }  // ❌ Dynamic
    ]
}

withName: '.*:MULTIQC_KALLISTO' {
    publishDir = [
        path: { "${params.outdir}/multiqc/${params.pseudo_aligner}" }  // ❌ Dynamic
    ]
}
```

### After (CORRECT):
```groovy
withName: '.*:MULTIQC_STAR' {
    publishDir = [
        path: { "${params.outdir}/multiqc/star" }  // ✅ Explicit
    ]
}

withName: '.*:MULTIQC_HISAT2' {
    publishDir = [
        path: { "${params.outdir}/multiqc/hisat2" }  // ✅ Explicit
    ]
}

withName: '.*:MULTIQC_KALLISTO' {
    publishDir = [
        path: { "${params.outdir}/multiqc/kallisto" }  // ✅ Explicit
    ]
}
```

## Why This Matters

### Problem with Dynamic Parameters:
1. **Parameter dependency**: Path depends on `params.aligner` or `params.pseudo_aligner` being set correctly
2. **Inconsistent behavior**: If parameters aren't set as expected, reports could go to wrong locations
3. **Process name mismatch**: The process name is `MULTIQC_STAR` but the path uses `${params.aligner}` which might not always be "star"
4. **Debugging difficulty**: Harder to trace where outputs actually end up

### Benefits of Hardcoded Paths:
1. **Predictable**: Always know where each report will be published
2. **Process-specific**: Path matches the process name (MULTIQC_STAR → multiqc/star)
3. **Parameter-independent**: Works regardless of parameter values
4. **Consistent with main config**: Matches the fix we made in `nextflow.config`

## Configuration Override Chain

Understanding why this fix was needed:

```
1. Main nextflow.config
   └─> Sets publishDir for MULTIQC_STAR, MULTIQC_HISAT2, MULTIQC_KALLISTO
       (Using hardcoded paths: multiqc/star, multiqc/hisat2, multiqc/kallisto)

2. workflows/rnaseq/nextflow.config (loaded via includeConfig)
   └─> OVERRIDES the above with dynamic parameter-based paths
       (Was using: multiqc/${params.aligner}, multiqc/${params.pseudo_aligner})

3. Result: Workflow config wins, causing the issue you saw
```

## Files Modified

- ✅ `nextflow.config` (commit `59c44d1`) - Initial fix with hardcoded paths
- ✅ `workflows/rnaseq/nextflow.config` (commit `b2ab566`) - Fixed to match main config

## Expected Output Structure

```
results/
├── multiqc/
│   ├── star/
│   │   └── star_multiqc_report.html
│   ├── hisat2/
│   │   └── hisat2_multiqc_report.html
│   ├── kallisto/
│   │   └── kallisto_multiqc_report.html
│   └── multiqc_report.html (combined report)
├── star/
│   ├── rsem/
│   │   ├── deseq2/
│   │   └── deeptools/
│   └── genome/
│       ├── deseq2/
│       └── deeptools/
├── hisat2/
│   └── ... (similar structure)
└── kallisto/
    ├── deseq2/
    └── deeptools/
```

## Complete Fix Summary

All publishDir issues have now been addressed across both config files:

1. **MultiQC aligner-specific reports** - Fixed in both configs
2. **DeepTools for pseudo-aligners** - Fixed in both configs  
3. **DESeq2 outputs** - Already correct (used as reference pattern)

## Testing

To verify everything works:

```bash
# Test with Kallisto
nextflow run . --pseudo_aligner kallisto --skip_alignment --outdir results

# Expected:
# results/multiqc/kallisto/kallisto_multiqc_report.html
# results/kallisto/deseq2/
# results/kallisto/deeptools/

# Test with STAR
nextflow run . --aligner star --outdir results

# Expected:
# results/multiqc/star/star_multiqc_report.html
# results/star/rsem/deseq2/
# results/star/rsem/deeptools/
```

## Related Commits

- `59c44d1` - Fix MultiQC publishDir for aligner-specific reports (main config)
- `b2ab566` - Fix MultiQC publishDir in workflow config to use hardcoded paths
- `d7c3bef` - Fix DeepTools publishDir for pseudo-aligners (main config)
- `f0e741b` - Add DeepTools publishDir for pseudo-aligners in workflow config
