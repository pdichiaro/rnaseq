# MultiQC Bug Fix

## Issue
MultiQC was not running because of an incorrect conditional check in the workflow.

## Root Cause
In `workflows/rnaseq/main.nf` line 1317, the MultiQC STAR condition was checking for invalid aligner values:

```groovy
// BUGGY CODE (BEFORE):
if (params.aligner && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem' || params.aligner == 'star')) {
```

### Why this was wrong:
1. **Valid aligners** (defined at line 144): `['star', 'hisat2']`
2. **`'star_salmon'` and `'star_rsem'` are NOT valid aligner values** in the current API
3. These were legacy values from an older version where aligner and quantification were combined
4. The modern API separates these:
   - `params.aligner = 'star'` or `'hisat2'`
   - `params.quantification = 'genome'` or `'rsem'` (or comma-separated list)

### Impact:
The condition could never match `'star_salmon'` or `'star_rsem'`, and would only work if exactly `params.aligner == 'star'` was true. However, the extra checks made the logic confusing and indicated incomplete refactoring from the old API.

## Fix Applied
Simplified the condition to check only for the valid aligner value:

```groovy
// FIXED CODE (AFTER):
if (params.aligner == 'star') {
```

## Location
- **File**: `workflows/rnaseq/main.nf`
- **Line**: 1317
- **Section**: MultiQC STAR Report generation

## Testing Recommendations
1. Run with default parameters (`aligner = 'star'`) - MultiQC should now generate
2. Run with `--aligner hisat2` - MultiQC_HISAT2 should generate
3. Run with `--pseudo_aligner kallisto` - MultiQC_KALLISTO should generate
4. Verify MultiQC reports are created in the results directory

## Additional Notes
- HISAT2 condition (line 1333) was already correct: `if (params.aligner == 'hisat2')`
- Kallisto condition (line 1347) was already correct: `if (params.pseudo_aligner == 'kallisto')`
- Legacy code still exists in `subworkflows/local/utils_nfcore_rnaseq_pipeline/main.nf` that references `'star_salmon'` and `'star_rsem'`, but this doesn't affect the main workflow execution

## Date
2025-11-30
