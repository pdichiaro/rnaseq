# Quick Reference: mmrnaseq Integration

## What Changed?

We integrated the mmrnaseq normalization strategy into nf-core/rnaseq for BigWig generation, making the pipelines consistent in their approach to scaling factor computation and usage.

## Key Points

### 🎯 Main Goal
- **Before**: BigWig normalization computed scaling factors independently
- **After**: Reuse scaling factors computed by DESeq2, exactly like mmrnaseq

### 🔧 Technical Changes

1. **Process Input** (`modules/local/deeptools_bw_norm/main.nf`)
   ```groovy
   // Before
   tuple val(meta), path(bam), path(bai)
   
   // After  
   tuple val(meta), path(bam), path(bai), val(scaling)
   ```

2. **Channel Pattern** (`workflows/rnaseq/main.nf`)
   ```groovy
   // mmrnaseq strategy: combine + filter
   ch_bam_for_deeptools
       .combine(ch_scaling_per_sample)
       .map { meta, bam, bai, sample_id, scaling -> 
           meta.id == sample_id ? [meta, bam, bai, scaling] : null
       }
       .filter { it != null }
   ```

3. **Scaling Factor Extraction**
   ```groovy
   // Read value from file directly in channel
   ch_scaling_factors_individual
       .flatten()
       .map { file ->
           def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
           def scaling_value = file.text.trim()
           [sample_name, scaling_value]
       }
   ```

### 📊 Comparison with mmrnaseq

| Feature | mmrnaseq | nf-core/rnaseq |
|---------|----------|----------------|
| Channel pattern | `.combine()` + filter | ✅ Same |
| Scaling value passing | `val(scaling)` | ✅ Same |
| Process signature | `(meta, bam, bai, scaling)` | ✅ Same |
| Multiple methods | ❌ Single output | ✅ invariant + all genes |

## Usage

### Default Behavior
```bash
# Uses 'all_genes' normalization (default)
nextflow run nf-core/rnaseq \\
    --input samplesheet.csv \\
    --outdir results \\
    --genome GRCh38
```

### Multiple Normalization Methods
```bash
# Generate BigWigs with both methods
nextflow run nf-core/rnaseq \\
    --input samplesheet.csv \\
    --normalization_method 'invariant_genes,all_genes' \\
    --outdir results
```

### Skip BigWig Normalization
```bash
# Skip entirely
nextflow run nf-core/rnaseq \\
    --skip_deeptools_norm
```

## File Locations

### Input (Scaling Factors)
```
results/deseq2/all_genes/
├── sample1_scaling_factor.txt
├── sample2_scaling_factor.txt
└── invariant_genes/
    ├── sample1_scaling_factor.txt
    └── sample2_scaling_factor.txt
```

### Output (BigWigs)
```
results/bigwig/deeptools_norm/
├── invariant_genes/
│   └── *.norm.bw
└── all_genes/
    └── *.norm.bw
```

## Benefits

✅ **Efficiency**: Compute scaling factors once, reuse multiple times  
✅ **Consistency**: Same strategy as mmrnaseq  
✅ **Flexibility**: Support for multiple normalization methods  
✅ **Maintainability**: Cleaner code with explicit filtering  
✅ **Robustness**: Better error handling with combine + filter pattern  

## Testing

### Syntax Check
```bash
cd rnaseq
nextflow run . --help | grep normalization
```

Expected output:
```
  --normalization_method        [string]  Normalization method(s)...
  --skip_deeptools_norm         [boolean] Skip deeptools bigwig normalization.
```

### Full Test (requires Docker/Singularity)
```bash
nextflow run . -profile test,docker --outdir test_results
```

## Documentation

- **Technical Details**: See `IMPLEMENTATION_SUMMARY.md`
- **Complete Changes**: See `CHANGE_SUMMARY.md`
- **mmrnaseq Reference**: https://github.com/fgualdr/mmrnaseq

## Commits

```bash
git log --oneline -2

a6e2309 refactor: adopt mmrnaseq channel combination strategy
8f1d134 feat: integrate mmrnaseq normalization strategy for BigWig generation
```

## Code Changes Summary

```bash
git diff cde0742..HEAD --stat

 IMPLEMENTATION_SUMMARY.md               | 234 +++++++++++++++++++
 modules/local/deeptools_bw_norm/main.nf |  87 ++++----
 workflows/rnaseq/main.nf                |  44 ++--
 3 files changed, 277 insertions(+), 88 deletions(-)
```

## Next Steps

1. ✅ Implementation complete
2. ✅ Documentation complete  
3. 🔄 Ready for full pipeline testing
4. 🔄 Ready for comparison with mmrnaseq outputs
5. 🔄 Ready for PR submission (if desired)

## Questions?

- Review `IMPLEMENTATION_SUMMARY.md` for technical implementation details
- Review `CHANGE_SUMMARY.md` for comprehensive change documentation
- Compare with mmrnaseq implementation in `mmrnaseq/workflows/rnaseq.nf` lines 699-729

---

**Last Updated**: 2025-11-17  
**Version**: nf-core/rnaseq 3.21.0+  
**Status**: ✅ Ready for testing
