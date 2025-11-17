# Fix: Missing Salmon Index Parameter

## Issue
The `FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS` subworkflow call was missing the required `ch_salmon_index` parameter, causing a parameter count mismatch error:
- **Expected:** 24 parameters
- **Received:** 23 parameters

## Root Cause
When salmon quantification was removed from the pipeline (commit `9f1979f`), the salmon index channel was not accounted for in the subworkflow call. However, the subworkflow definition still requires this parameter (even if unused) for strandedness detection compatibility.

## Solution
Added `Channel.empty()` as the 5th parameter to represent the unused salmon index channel:

```groovy
FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
    ch_fastq,
    ch_fasta,
    ch_transcript_fasta,
    ch_gtf,
    Channel.empty(),  // ch_salmon_index - not used, salmon quantification removed
    ch_sortmerna_index,
    ch_bbsplit_index,
    // ... remaining parameters
)
```

## Verification
✅ Parameter count now matches: 24 passed = 24 expected
✅ Nextflow config validation passes
✅ All parameters correctly mapped to subworkflow expectations

## Parameter Mapping (After Fix)
1. ch_fastq → ch_reads
2. ch_fasta → ch_fasta
3. ch_transcript_fasta → ch_transcript_fasta
4. ch_gtf → ch_gtf
5. **Channel.empty() → ch_salmon_index** ⬅️ **FIXED**
6. ch_sortmerna_index → ch_sortmerna_index
7. ch_bbsplit_index → ch_bbsplit_index
8. ch_ribo_db → ch_rrna_fastas
9-24. (remaining parameters correctly mapped)

## Files Modified
- `workflows/rnaseq/main.nf` (line 239)

## Testing Recommendation
Test with minimal dataset to verify:
1. Pipeline initialization succeeds
2. QC/trim/filter subworkflow executes without parameter errors
3. Strandedness detection operates correctly without salmon
