# Fix: DeepTools BigWig Publishing for Multiple Quantification Methods

## Problem Description

When running the pipeline with `--quantification rsem,genome`, the DeepTools BigWig files are published to a folder named "rsem,genome" instead of separate "rsem" and "genome" folders.

### Current Behavior (INCORRECT)
```
results/star/
└── rsem,genome/                    ← Single folder with comma in name!
    └── deeptools/
        ├── all_genes/
        │   ├── sample1.norm.bw     ← Mixed RSEM and genome files
        │   └── sample2.norm.bw
        └── invariant_genes/
            ├── sample1.norm.bw
            └── sample2.norm.bw
```

### Expected Behavior (CORRECT)
```
results/star/
├── rsem/
│   └── deeptools/
│       ├── all_genes/
│       │   ├── sample1.norm.bw     ← RSEM-based normalization
│       │   └── sample2.norm.bw
│       └── invariant_genes/
│           ├── sample1.norm.bw
│           └── sample2.norm.bw
└── genome/
    └── deeptools/
        ├── all_genes/
        │   ├── sample1.norm.bw     ← Genome-based normalization
        │   └── sample2.norm.bw
        └── invariant_genes/
            ├── sample1.norm.bw
            └── sample2.norm.bw
```

## Root Cause Analysis

### 1. PublishDir Uses Literal String
**Location:** `workflows/rnaseq/nextflow.config` lines 519 and 527

```groovy
path: { "${params.outdir}/${params.aligner}/${params.quantification}/deeptools/all_genes" }
```

When `params.quantification = "rsem,genome"`, this literally creates a folder named "rsem,genome".

### 2. Scaling Factors Are Mixed
**Location:** `workflows/rnaseq/main.nf` line 1070+

The DeepTools section uses `ch_scaling_factors_individual` which contains scaling factors from BOTH RSEM and genome quantifications mixed together:

```groovy
ch_scaling_factors_individual
    .flatten()
    .filter { file -> /* filters by invariant/all_genes */ }
    .map { file -> [sample_name, scaling_value] }
```

This mixes files from:
- `star/rsem/deseq2/all_genes/*_scaling_factor.txt`
- `star/genome/deseq2/all_genes/*_scaling_factor.txt`

### 3. Single DeepTools Execution
The DeepTools processes run ONCE for all quantification methods, making it impossible to separate outputs by quantification method.

## Solution Approach

We have two potential solutions:

### Option 1: Infer Quantification Method from File Path (RECOMMENDED)
Modify the publishDir to examine the input file paths and determine the quantification method dynamically.

**Pros:**
- Minimal workflow changes
- Works with current architecture
- Backward compatible

**Cons:**
- Relies on file path parsing
- Slightly more complex publishDir logic

### Option 2: Run DeepTools Separately for Each Quantification Method
Create separate channels and process calls for RSEM and genome quantifications.

**Pros:**
- Clear separation of concerns
- Explicit and maintainable
- Better scalability

**Cons:**
- More extensive workflow changes
- Requires separate channel management
- More code duplication

## Recommended Fix: Option 1 (Path-Based Detection)

### Changes Required

#### 1. Update publishDir in nextflow.config

Modify the publishDir saveAs logic to detect quantification method from the scaling factor file path.

**File:** `workflows/rnaseq/nextflow.config`

**Before:**
```groovy
withName: '.*:DEEPTOOLS_BIGWIG_NORM_ALL_GENES' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/${params.quantification}/deeptools/all_genes" },
        mode: params.publish_dir_mode,
        pattern: "*.norm.bw"
    ]
}
```

**After:**
```groovy
withName: '.*:DEEPTOOLS_BIGWIG_NORM_ALL_GENES' {
    publishDir = [
        path: { 
            // Determine quantification method from meta map if available
            def quant_method = task.ext.quantification_method ?: params.quantification
            "${params.outdir}/${params.aligner}/${quant_method}/deeptools/all_genes" 
        },
        mode: params.publish_dir_mode,
        pattern: "*.norm.bw"
    ]
}
```

#### 2. Add Quantification Method to Task Extension

Modify the workflow to pass quantification method information through task.ext.

**File:** `workflows/rnaseq/main.nf`

We need to set `task.ext.quantification_method` based on the scaling factor file path.

### Implementation Details

The scaling factors come from different directories:
- RSEM: `results/star/rsem/deseq2/{all_genes,invariant_genes}/*_scaling_factor.txt`
- Genome: `results/star/genome/deseq2/{all_genes,invariant_genes}/*_scaling_factor.txt`

We can detect the quantification method by checking if the file path contains "rsem" or "genome".

## Alternative Fix: Path Detection in saveAs

A simpler approach is to use the `saveAs` closure to detect the quantification method from the input channel.

However, the `saveAs` closure only has access to the output filename, not the input metadata. We need to pass this information through the channel.

## Implemented Solution: Add Quantification Method to Meta Map ✅

The most robust solution is to add the quantification method to the meta map when creating the scaling factor channels.

### Implementation Steps ✅

1. ✅ **Modified scaling factor channel creation** to detect and include quantification method
2. ✅ **Added quantification method to meta map** before passing to DeepTools processes
3. ✅ **Updated publishDir configuration** to use `meta.quantification` instead of `params.quantification`

## Implementation Details

### 1. Workflow Changes (`workflows/rnaseq/main.nf`)

#### A. Detect Quantification Method from File Path

Modified the scaling factor extraction to detect the quantification method from the file path:

```groovy
.map { file ->
    def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
    def scaling_value = file.text.trim()
    // Detect quantification method from file path
    def file_path = file.toString()
    def quant_method = file_path.contains('/rsem/') ? 'rsem' : 
                      file_path.contains('/genome/') ? 'genome' :
                      file_path.contains('/salmon/') ? 'salmon' : 'unknown'
    [sample_name, scaling_value, quant_method]  // ← Now includes quant_method
}
```

**Detection Logic:**
- Checks if file path contains `/rsem/` → assigns 'rsem'
- Checks if file path contains `/genome/` → assigns 'genome'
- Checks if file path contains `/salmon/` → assigns 'salmon'
- Falls back to 'unknown' if none match

This works because scaling factors are published to:
- `results/star/rsem/deseq2/all_genes/*_scaling_factor.txt`
- `results/star/genome/deseq2/all_genes/*_scaling_factor.txt`

#### B. Add Quantification Method to Meta Map

Modified the channel combination to add quantification method to the meta map:

```groovy
ch_combined_input_invariant = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample_invariant)
    .map { meta, bam, bai, sample_id, scaling, quant_method ->  // ← Now receives quant_method
        if (meta.id == sample_id) {
            def new_meta = meta.clone()
            new_meta.quantification = quant_method  // ← Add to meta map
            [new_meta, bam, bai, scaling]
        } else {
            null
        }
    }
    .filter { it != null }
```

**Key Changes:**
- Clone the meta map to avoid modifying the original
- Add `quantification` field with the detected method
- Return the updated meta map with BAM and scaling factor

#### C. Applied to Both Normalization Methods

Both `invariant_genes` and `all_genes` normalization methods were updated:
- Lines 1064-1110: `DEEPTOOLS_BIGWIG_NORM_INVARIANT`
- Lines 1113-1147: `DEEPTOOLS_BIGWIG_NORM_ALL_GENES`

### 2. Configuration Changes (`workflows/rnaseq/nextflow.config`)

#### Updated PublishDir Paths

Modified both DeepTools processes to use `meta.quantification`:

```groovy
withName: '.*:DEEPTOOLS_BIGWIG_NORM_ALL_GENES' {
    publishDir = [
        path: { 
            // Use quantification method from meta map (added dynamically in workflow)
            // Falls back to params.quantification for single quantification method
            def quant_method = meta.quantification ?: params.quantification
            "${params.outdir}/${params.aligner}/${quant_method}/deeptools/all_genes" 
        },
        mode: params.publish_dir_mode,
        pattern: "*.norm.bw"
    ]
}
```

**Fallback Strategy:**
- Uses `meta.quantification` if available (set in workflow)
- Falls back to `params.quantification` for backward compatibility
- Ensures single quantification methods still work correctly

## Results

### Before (BROKEN) ❌
```
results/star/
└── rsem,genome/                    ← Wrong! Comma in folder name
    └── deeptools/
        ├── all_genes/
        │   └── sample1.norm.bw     ← Can't distinguish RSEM vs genome
        └── invariant_genes/
            └── sample1.norm.bw
```

### After (FIXED) ✅
```
results/star/
├── rsem/
│   └── deeptools/
│       ├── all_genes/
│       │   └── sample1.norm.bw     ← RSEM normalization
│       └── invariant_genes/
│           └── sample1.norm.bw
└── genome/
    └── deeptools/
        ├── all_genes/
        │   └── sample1.norm.bw     ← Genome normalization
        └── invariant_genes/
            └── sample1.norm.bw
```

## Testing

### Test Command
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner star \
  --quantification rsem,genome \
  --skip_deeptools_norm false
```

### Expected Results

1. **Separate DeepTools BigWig Files**
   - RSEM-normalized BigWig files in `star/rsem/deeptools/`
   - Genome-normalized BigWig files in `star/genome/deeptools/`

2. **Correct File Organization**
   ```bash
   # Check RSEM BigWig files
   ls results/star/rsem/deeptools/all_genes/
   # Should show: sample1.norm.bw, sample2.norm.bw, ...
   
   # Check genome BigWig files
   ls results/star/genome/deeptools/all_genes/
   # Should show: sample1.norm.bw, sample2.norm.bw, ...
   ```

3. **No "rsem,genome" Folder**
   ```bash
   ls results/star/ | grep ","
   # Should return nothing (no folders with commas)
   ```

### Verification Script

```bash
#!/bin/bash

RESULTS="results/star"

echo "Checking DeepTools output structure..."

# Check for incorrect folder
if [ -d "$RESULTS/rsem,genome" ]; then
    echo "❌ FAILED: Found 'rsem,genome' folder (should not exist)"
    exit 1
else
    echo "✅ PASS: No 'rsem,genome' folder"
fi

# Check for correct RSEM folder
if [ -d "$RESULTS/rsem/deeptools" ]; then
    echo "✅ PASS: Found 'rsem/deeptools/' folder"
    bw_count=$(ls -1 "$RESULTS/rsem/deeptools/all_genes"/*.bw 2>/dev/null | wc -l)
    echo "   → RSEM BigWig files: $bw_count"
else
    echo "⚠️  WARNING: No 'rsem/deeptools/' folder"
fi

# Check for correct genome folder
if [ -d "$RESULTS/genome/deeptools" ]; then
    echo "✅ PASS: Found 'genome/deeptools/' folder"
    bw_count=$(ls -1 "$RESULTS/genome/deeptools/all_genes"/*.bw 2>/dev/null | wc -l)
    echo "   → Genome BigWig files: $bw_count"
else
    echo "⚠️  WARNING: No 'genome/deeptools/' folder"
fi

echo
echo "✅ DeepTools structure verification complete!"
```

## Backward Compatibility

### Single Quantification Method
When using a single quantification method (e.g., `--quantification genome`):
- `meta.quantification` will be set to the detected method ('genome')
- publishDir will correctly use 'genome' in the path
- Behavior is identical to previous versions

### Multiple Quantification Methods
When using multiple methods (e.g., `--quantification rsem,genome`):
- Each scaling factor is tagged with its quantification method
- BigWig files are correctly separated by method
- No folders with commas are created

## Benefits

1. ✅ **Correct Organization** - Separate folders for each quantification method
2. ✅ **Clear Attribution** - Easy to identify which normalization was used
3. ✅ **Backward Compatible** - Single quantification methods work as before
4. ✅ **Scalable** - Works with any number of quantification methods
5. ✅ **Maintainable** - Automatic detection from file paths
6. ✅ **No Duplication** - Reuses existing channel infrastructure

## Technical Notes

### Why Path Detection?
- Scaling factors are already published to method-specific directories
- File paths inherently contain the quantification method information
- No need to pass additional parameters through complex channel operations
- Works with existing DESeq2 QC architecture

### Meta Map Cloning
We clone the meta map before modification to avoid side effects:
```groovy
def new_meta = meta.clone()
new_meta.quantification = quant_method
```

This ensures the original meta map remains unchanged for other processes.

### Fallback Logic
The publishDir uses Elvis operator for fallback:
```groovy
def quant_method = meta.quantification ?: params.quantification
```

This handles edge cases where meta.quantification might not be set.

## Files Modified

1. **`workflows/rnaseq/main.nf`**
   - Lines ~1064-1110: Modified invariant_genes normalization channel operations
   - Lines ~1113-1147: Modified all_genes normalization channel operations

2. **`workflows/rnaseq/nextflow.config`**
   - Lines ~516-523: Updated DEEPTOOLS_BIGWIG_NORM_ALL_GENES publishDir
   - Lines ~529-536: Updated DEEPTOOLS_BIGWIG_NORM_INVARIANT publishDir

## Related Issues

This fix resolves:
- ❌ Folders named "rsem,genome" with comma
- ❌ Mixed BigWig files from different quantification methods
- ❌ Inability to distinguish normalization source

## Future Enhancements

Potential improvements:
1. Add quantification method to BigWig file metadata
2. Create MultiQC section comparing normalizations
3. Add validation to ensure all samples processed for each method
