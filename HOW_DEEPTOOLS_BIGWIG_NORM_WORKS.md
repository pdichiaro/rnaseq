# How DEEPTOOLS_BIGWIG_NORM Works Now

## Overview

The `DEEPTOOLS_BIGWIG_NORM` process generates normalized BigWig files from BAM alignments using pre-computed scaling factors from DESeq2. This follows the **mmrnaseq strategy** where normalization factors are computed once and reused.

---

## 🔄 Complete Workflow

### Step 1: Scaling Factor Computation (DESeq2)

**Process**: `DESEQ2_QC_SALMON` (or other DESeq2 processes)

```
Input: Read count matrices
↓
DESeq2 normalization (GeneralNormalizer with optional invariant genes)
↓
Output: Individual scaling factor files
```

**Output Files**:
```
results/deseq2/all_genes/
├── sample1_scaling_factor.txt  (contains: "0.8532")
├── sample2_scaling_factor.txt  (contains: "1.2341")
└── invariant_genes/
    ├── sample1_scaling_factor.txt  (contains: "0.7821")
    └── sample2_scaling_factor.txt  (contains: "1.1543")
```

Each file contains a **single numeric value** (the scaling factor for that sample).

---

### Step 2: Channel Preparation (Workflow)

**Location**: `workflows/rnaseq/main.nf` lines ~1122-1185

#### 2A: Extract Scaling Factors from Files

```groovy
ch_scaling_per_sample_invariant = ch_scaling_factors_individual
    .flatten()                          // Flatten collection of files
    .filter { file ->                   // Filter by directory structure
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file ->                      // Extract sample name and value
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        def scaling_value = file.text.trim()  // 🔑 READ VALUE HERE
        [sample_name, scaling_value]
    }
```

**Result**: Channel emitting `[sample_id, scaling_value]` tuples
```
['sample1', '0.7821']
['sample2', '1.1543']
...
```

#### 2B: Prepare BAM Files with Indices

```groovy
ch_bam_for_deeptools = ch_genome_bam
    .join(ch_genome_bam_index, by: [0])
```

**Result**: Channel emitting `[meta, bam, bai]` tuples
```
[[id:'sample1', ...], sample1.bam, sample1.bam.bai]
[[id:'sample2', ...], sample2.bam, sample2.bam.bai]
...
```

#### 2C: Combine BAM with Scaling Factors (mmrnaseq Strategy)

```groovy
ch_combined_input = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample)     // Cartesian product
    .map { meta, bam, bai, sample_id, scaling -> 
        meta.id == sample_id ? [meta, bam, bai, scaling] : null
    }
    .filter { it != null }               // Keep only matches
```

**How it works**:
1. `.combine()` creates all possible combinations of BAM files and scaling factors
2. `.map()` checks if sample IDs match: `meta.id == sample_id`
3. If match → emit `[meta, bam, bai, scaling]`
4. If no match → emit `null`
5. `.filter()` removes all `null` values

**Result**: Channel emitting `[meta, bam, bai, scaling]` tuples
```
[[id:'sample1', ...], sample1.bam, sample1.bam.bai, '0.7821']
[[id:'sample2', ...], sample2.bam, sample2.bam.bai, '1.1543']
...
```

---

### Step 3: BigWig Generation (Process)

**Location**: `modules/local/deeptools_bw_norm/main.nf`

#### Process Signature

```groovy
process DEEPTOOLS_BIGWIG_NORM {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(bai), val(scaling)  // 🔑 scaling as VALUE

    output:
    path "*.unstranded.norm.bw" , emit: unstranded_bw
    path "*.fwd.norm.bw"        , optional:true, emit: fw_bw
    path "*.rev.norm.bw"        , optional:true, emit: rev_bw
    path "versions.yml"         , emit: versions
```

#### Key Points:
1. **`val(scaling)`**: Receives scaling factor as a **value** (not a file)
2. **Direct usage**: Can use `$scaling` directly in script
3. **No file I/O**: No need to read from a file during execution

---

## 🎯 Script Logic

The process generates different BigWig files based on:
1. **Strandedness**: unstranded, forward, reverse
2. **Read type**: single-end vs paired-end

### Example 1: Unstranded Library

```bash
# Input: meta.strandedness = 'unstranded'
# Output: Only unstranded BigWig

bamCoverage \
    --numberOfProcessors $task.cpus \
    --binSize 1 \
    --scaleFactor $scaling \    # 🔑 Use scaling value directly
    --bam $bam \
    -o ${prefix}.unstranded.norm.bw
```

### Example 2: Single-End Forward Library

```bash
# Input: meta.strandedness = 'forward', meta.single_end = true
# Output: Forward and Reverse BigWigs

# Forward strand (exclude reverse-mapped reads)
bamCoverage -b $bam \
    --scaleFactor $scaling \
    --samFlagExclude 16 \
    -o ${prefix}.fwd.norm.bw

# Reverse strand (include only reverse-mapped reads)
bamCoverage -b $bam \
    --scaleFactor $scaling \
    --samFlagInclude 16 \
    -o ${prefix}.rev.norm.bw
```

### Example 3: Paired-End Reverse Library

```bash
# Input: meta.strandedness = 'reverse', meta.single_end = false
# Output: Forward, Reverse, and Unstranded BigWigs

# 1. Generate unstranded BigWig
bamCoverage --bam $bam --scaleFactor $scaling -o ${prefix}.unstranded.norm.bw

# 2. Create FORWARD strand BigWig (FLIPPED for reverse library)
# For reverse library: reads mapped as 'reverse' represent forward transcripts
samtools view -b -f 144 $bam > ${prefix}.fwd1.bam  # R2 mapped to reverse
samtools view -b -f 64 -F 16 $bam > ${prefix}.fwd2.bam  # R1 mapped to forward
samtools merge -f ${prefix}.fwd.bam ${prefix}.fwd1.bam ${prefix}.fwd2.bam
samtools index ${prefix}.fwd.bam
bamCoverage --bam ${prefix}.fwd.bam --scaleFactor $scaling -o ${prefix}.fwd.norm.bw

# 3. Create REVERSE strand BigWig (FLIPPED for reverse library)
# For reverse library: reads mapped as 'forward' represent reverse transcripts
samtools view -b -f 128 -F 16 $bam > ${prefix}.rev1.bam  # R2 mapped to forward
samtools view -b -f 80 $bam > ${prefix}.rev2.bam  # R1 mapped to reverse
samtools merge -f ${prefix}.rev.bam ${prefix}.rev1.bam ${prefix}.rev2.bam
samtools index ${prefix}.rev.bam
bamCoverage --bam ${prefix}.rev.bam --scaleFactor $scaling -o ${prefix}.rev.norm.bw
```

---

## 🔍 Key Technical Details

### Scaling Factor Application

**What bamCoverage does with `--scaleFactor`**:
```
normalized_coverage = (raw_coverage / total_reads) * scaleFactor
```

**Example**:
- Sample has 10M reads
- Gene region has 1000 reads
- Scaling factor: 0.8532

```
Raw coverage = 1000
Coverage per million = 1000 / 10M = 0.0001
Normalized = 0.0001 * 0.8532 = 0.00008532
```

### SAM Flags for Strand Separation

#### Single-End Reads:
- **Flag 16**: Read mapped to reverse strand
- **Flag 0**: Read mapped to forward strand

#### Paired-End Reads:
- **Flag 64**: First in pair (R1)
- **Flag 128**: Second in pair (R2)
- **Flag 16**: Mapped to reverse strand

**Combinations**:
- `64 + 16 = 80`: R1 mapped to reverse
- `128 + 16 = 144`: R2 mapped to reverse
- `-f 64 -F 16`: R1 mapped to forward
- `-f 128 -F 16`: R2 mapped to forward

### Library Strandedness Logic

| Library Type | R1 Orientation | R2 Orientation | Forward BigWig | Reverse BigWig |
|--------------|----------------|----------------|----------------|----------------|
| Unstranded   | N/A            | N/A            | Not generated  | Not generated  |
| Forward (PE) | Same as transcript | Opposite of transcript | R2 fwd + R1 rev | R2 rev + R1 fwd |
| Reverse (PE) | Opposite of transcript | Same as transcript | R2 rev + R1 fwd | R2 fwd + R1 rev |

**Why the flip for reverse libraries?**
- In reverse-stranded libraries, R2 (second in pair) has the same orientation as the transcript
- But when mapping, "forward" and "reverse" refer to the genome strand
- So we need to flip the logic to get the correct transcript strand

---

## 📊 Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: DESeq2 Normalization                               │
│  ────────────────────────────────────────────────────────── │
│  Count Matrix → DESeq2 → Scaling Factors                    │
│                                                              │
│  Output: sample1_scaling_factor.txt (contains "0.8532")     │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 2A: Extract Scaling Factors                           │
│  ────────────────────────────────────────────────────────── │
│  ch_scaling_factors_individual                              │
│    .flatten()                                               │
│    .filter { /* by directory */ }                           │
│    .map { file ->                                           │
│        def sample_name = file.name - '_scaling_factor.txt'  │
│        def scaling_value = file.text.trim()  ← READ HERE    │
│        [sample_name, scaling_value]                         │
│    }                                                         │
│                                                              │
│  Output: ['sample1', '0.8532']                              │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ├──────────────────┐
                       │                  │
                       ▼                  ▼
┌────────────────────────────────┐  ┌──────────────────────────┐
│  STEP 2B: BAM Channel          │  │  Scaling Channel         │
│  ──────────────────────────    │  │  ──────────────────────  │
│  [meta, bam, bai]              │  │  [sample_id, scaling]    │
└────────────────┬───────────────┘  └──────────┬───────────────┘
                 │                              │
                 └──────────┬───────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 2C: Combine (mmrnaseq strategy)                       │
│  ────────────────────────────────────────────────────────── │
│  .combine(ch_bam, ch_scaling)  → Cartesian product          │
│  .map { meta, bam, bai, id, scale ->                        │
│      meta.id == id ? [meta, bam, bai, scale] : null         │
│  }                                                           │
│  .filter { it != null }  → Keep only matches                │
│                                                              │
│  Output: [meta, bam, bai, '0.8532']                         │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: DEEPTOOLS_BIGWIG_NORM Process                      │
│  ────────────────────────────────────────────────────────── │
│  Input: tuple val(meta), path(bam), path(bai), val(scaling) │
│                                                              │
│  Script:                                                     │
│    echo "Scaling factor: $scaling"                          │
│    bamCoverage \                                            │
│      --scaleFactor $scaling \   ← USE VALUE DIRECTLY        │
│      --bam $bam \                                           │
│      -o ${prefix}.unstranded.norm.bw                        │
│                                                              │
│  Output: Normalized BigWig files                            │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  FINAL OUTPUT                                               │
│  ────────────────────────────────────────────────────────── │
│  results/bigwig/deeptools_norm/                             │
│  ├── invariant_genes/                                       │
│  │   ├── sample1.unstranded.norm.bw                        │
│  │   ├── sample1.fwd.norm.bw                               │
│  │   └── sample1.rev.norm.bw                               │
│  └── all_genes/                                             │
│      ├── sample1.unstranded.norm.bw                        │
│      ├── sample1.fwd.norm.bw                               │
│      └── sample1.rev.norm.bw                               │
└─────────────────────────────────────────────────────────────┘
```

---

## 🔑 Key Improvements Over Previous Implementation

### Before (File-based approach):
```groovy
// OLD: Add file path to meta
ch_combined_input = ch_bam
    .map { meta, bam, bai -> 
        [meta + [scaling_factor_file: scaling_file], bam, bai] 
    }

// Process had to read file
script:
"""
SCALING_FACTOR=\$(cat "${meta.scaling_factor_file}")
bamCoverage --scaleFactor \$SCALING_FACTOR ...
"""
```

**Problems**:
- ❌ File I/O during task execution
- ❌ More complex meta map
- ❌ Extra step to read file in script
- ❌ Harder to debug

### After (Value-based approach, mmrnaseq):
```groovy
// NEW: Read value in channel operation
ch_scaling_per_sample = ch_scaling_factors
    .map { file ->
        [sample_name, file.text.trim()]  // Read once here
    }

ch_combined_input = ch_bam
    .combine(ch_scaling_per_sample)
    .map { meta, bam, bai, id, scaling -> 
        meta.id == id ? [meta, bam, bai, scaling] : null
    }
    .filter { it != null }

// Process receives value directly
input:
tuple val(meta), path(bam), path(bai), val(scaling)

script:
"""
bamCoverage --scaleFactor $scaling ...
"""
```

**Benefits**:
- ✅ Read file once during channel preparation
- ✅ Simpler meta map (no added fields)
- ✅ Direct value usage in script
- ✅ Easier to debug (value visible in logs)
- ✅ Matches mmrnaseq implementation
- ✅ More efficient execution

---

## 🎯 Multiple Normalization Methods

The workflow supports running both normalization methods simultaneously:

```bash
nextflow run nf-core/rnaseq \
    --normalization_method 'invariant_genes,all_genes'
```

**What happens**:
1. Two separate channel preparation blocks run in parallel
2. Different filtering criteria for scaling factor files
3. Two separate process calls: `DEEPTOOLS_BIGWIG_NORM_INVARIANT` and `DEEPTOOLS_BIGWIG_NORM_ALL_GENES`
4. Outputs go to different directories

**File Filtering**:

For **invariant genes**:
```groovy
.filter { file ->
    def parent_dir = file.getParent()?.getName() ?: ""
    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
    parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
}
```

For **all genes**:
```groovy
.filter { file ->
    def parent_dir = file.getParent()?.getName() ?: ""
    def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
    !(parent_dir.contains('invariant') || grandparent_dir.contains('invariant'))
}
```

---

## 🔬 Example: Complete Flow for One Sample

### Sample: `HBR_Rep1`
### Normalization: `invariant_genes`

**1. DESeq2 creates scaling factor:**
```
File: results/deseq2/all_genes/invariant_genes/HBR_Rep1_scaling_factor.txt
Content: 0.8532
```

**2. Channel extracts value:**
```groovy
['HBR_Rep1', '0.8532']
```

**3. BAM channel provides:**
```groovy
[[id:'HBR_Rep1', single_end:false, strandedness:'reverse'], 
 HBR_Rep1.bam, 
 HBR_Rep1.bam.bai]
```

**4. Combine matches by ID:**
```groovy
[[id:'HBR_Rep1', single_end:false, strandedness:'reverse'], 
 HBR_Rep1.bam, 
 HBR_Rep1.bam.bai, 
 '0.8532']
```

**5. Process generates BigWigs:**
```bash
# Unstranded
bamCoverage --bam HBR_Rep1.bam --scaleFactor 0.8532 \
    -o HBR_Rep1.unstranded.norm.bw

# Forward (flipped for reverse library)
samtools view -b -f 144 HBR_Rep1.bam > HBR_Rep1.fwd1.bam
samtools view -b -f 64 -F 16 HBR_Rep1.bam > HBR_Rep1.fwd2.bam
samtools merge HBR_Rep1.fwd.bam HBR_Rep1.fwd1.bam HBR_Rep1.fwd2.bam
bamCoverage --bam HBR_Rep1.fwd.bam --scaleFactor 0.8532 \
    -o HBR_Rep1.fwd.norm.bw

# Reverse (flipped for reverse library)
[similar process]
bamCoverage --bam HBR_Rep1.rev.bam --scaleFactor 0.8532 \
    -o HBR_Rep1.rev.norm.bw
```

**6. Output files:**
```
results/bigwig/deeptools_norm/invariant_genes/
├── HBR_Rep1.unstranded.norm.bw
├── HBR_Rep1.fwd.norm.bw
└── HBR_Rep1.rev.norm.bw
```

---

## 📝 Summary

**DEEPTOOLS_BIGWIG_NORM now works as follows**:

1. **Receives pre-computed scaling factors** from DESeq2 as values (not files)
2. **Matches samples** using the mmrnaseq `.combine()` strategy
3. **Applies scaling factors** directly to bamCoverage without file I/O
4. **Generates strand-specific BigWigs** based on library type
5. **Supports multiple normalization methods** (invariant genes, all genes)
6. **More efficient** by reading scaling factors once during channel preparation

This implementation is **identical to mmrnaseq** in strategy while adding support for multiple normalization methods simultaneously.
