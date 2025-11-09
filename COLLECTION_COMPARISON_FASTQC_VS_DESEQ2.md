# Collection Pattern Comparison: FastQC vs DESeq2

## 🔍 Key Finding

**FastQC and DESeq2 use COMPLETELY DIFFERENT collection patterns:**

- **FastQC**: Outputs are collected **WITHIN the subworkflow**, then emitted as a flattened channel
- **DESeq2**: Outputs are mixed **DIRECTLY** into the main workflow's `ch_multiqc_files` channel

This explains why DESeq2 needed `.flatten()` before the final `.collect()` while FastQC didn't!

---

## 📊 Detailed Comparison

### FastQC Collection Pattern

#### Location
`subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf`

#### Pattern
```groovy
workflow FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS {
    take:
    ch_reads
    ...
    
    main:
    ch_multiqc_files = Channel.empty()
    
    // Files are mixed INTERNALLY within the subworkflow
    if (trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE(...)
        
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .mix(ch_multiqc_files)
    }
    
    if (trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP(...)
        
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
            .mix(ch_multiqc_files)
    }
    
    emit:
    // ⭐ KEY: Files are TRANSPOSED and MAPPED before emit
    multiqc_files = ch_multiqc_files.transpose().map { it[1] }
    ...
}
```

#### In Main Workflow (workflows/rnaseq/main.nf)
```groovy
FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS(...)

// ⭐ KEY: Entire subworkflow output mixed as ONE operation
ch_multiqc_files = ch_multiqc_files.mix(
    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files
)
```

#### Key Operations
1. **`.transpose()`** - Converts `[ [meta1, file1], [meta2, file2] ]` → `[ meta1, file1, meta2, file2 ]`
2. **`.map { it[1] }`** - Extracts just the files: `[ file1, file2, file3, ... ]`
3. Files are **already flattened** before leaving the subworkflow!

---

### DESeq2 Collection Pattern

#### Location
`workflows/rnaseq/main.nf` (main workflow - NOT in a subworkflow)

#### Pattern
```groovy
// Run DESeq2 normalization processes
NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO(...)

// ⭐ KEY: Each output is mixed INDIVIDUALLY into the main channel
ch_multiqc_files = ch_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.sample_distances_txt
)
ch_multiqc_files = ch_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.pca_all_genes_txt
)
ch_multiqc_files = ch_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.pca_top_genes_txt
)
ch_multiqc_files = ch_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.read_dist_norm_txt
)

// Repeated for each DESeq2 process (invariant_genes + all_genes) × (Kallisto + Salmon + alignment) = 6 processes
// Each process emits 4 files
// Total: 6 processes × 4 files = 24 individual mix operations
```

#### Key Differences
1. **NO `.transpose()`** - Files come directly from process outputs
2. **NO `.map { it[1] }`** - Full output structure is mixed
3. Each process emits **4 files as separate outputs**, not grouped
4. Files are mixed **INDIVIDUALLY** into `ch_multiqc_files`

---

## 🧩 Why .flatten() Was Needed for DESeq2 But Not FastQC

### FastQC Channel Structure at MultiQC Input
```groovy
FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files:
[
  file1.fastqc.zip,      ← Already individual items
  file2.fastqc.zip,
  file3.fastp.json,
  file4.trim.log,
  ...
]

// When mixed into ch_multiqc_files
ch_multiqc_files = [
  file1.fastqc.zip,
  file2.fastqc.zip,
  file3.fastp.json,
  ...
]

// .collect() works perfectly!
```

### DESeq2 Channel Structure at MultiQC Input (BEFORE .flatten())
```groovy
After multiple .mix() operations with DESeq2 outputs:

ch_multiqc_files = [
  file1.fastqc.zip,                                          ← From FastQC
  file2.fastqc.zip,
  [                                                          ← From DESeq2 process
    sample1.deseq2.distances.txt,
    sample1.deseq2.pca_all.txt,
    sample1.deseq2.pca_top.txt,
    sample1.deseq2.read_dist.txt
  ],
  file3.log,
  [                                                          ← From another DESeq2 process
    sample2.deseq2.distances.txt,
    sample2.deseq2.pca_all.txt,
    ...
  ],
  ...
]

// .collect() creates: [ [file1, file2, [nested], file3, [nested], ...] ]
// MultiQC receives NESTED structure → Can't find files! ❌
```

### DESeq2 Channel Structure at MultiQC Input (AFTER .flatten())
```groovy
ch_multiqc_files.flatten():
[
  file1.fastqc.zip,
  file2.fastqc.zip,
  sample1.deseq2.distances.txt,      ← Now individual!
  sample1.deseq2.pca_all.txt,
  sample1.deseq2.pca_top.txt,
  sample1.deseq2.read_dist.txt,
  file3.log,
  sample2.deseq2.distances.txt,
  sample2.deseq2.pca_all.txt,
  ...
]

// .collect() creates: [ file1, file2, sample1.dist, sample1.pca, ... ]
// MultiQC receives FLAT structure → Finds all files! ✅
```

---

## 🔬 Technical Analysis

### Why Does This Nesting Happen?

#### Process Output Behavior
When a process emits multiple separate outputs:

```groovy
process NORMALIZE_DESEQ2_QC {
    output:
    path "*.sample.dists.txt"     , emit: sample_distances_txt
    path "*.pca.vals.txt"         , emit: pca_all_genes_txt
    path "*.pca.top.vals.txt"     , emit: pca_top_genes_txt
    path "*.read_dist_norm.txt"   , emit: read_dist_norm_txt
}
```

Each emit produces:
- `sample_distances_txt` → `[meta, [file1.txt]]` or potentially `[file1.txt, file2.txt, ...]`
- `pca_all_genes_txt` → `[meta, [file1.txt]]`
- etc.

#### Channel Mixing Behavior
When you mix channels repeatedly:
```groovy
ch_multiqc_files = ch_multiqc_files.mix(output1)
ch_multiqc_files = ch_multiqc_files.mix(output2)
ch_multiqc_files = ch_multiqc_files.mix(output3)
```

Nextflow doesn't automatically flatten the items being mixed. If `output1` contains multiple files (emitted together), they may remain grouped.

#### The Subworkflow Advantage
FastQC's subworkflow pattern:
```groovy
emit:
    multiqc_files = ch_multiqc_files.transpose().map { it[1] }
```

The `.transpose().map { it[1] }` explicitly:
1. Expands nested tuples: `[ [meta, file] ]` → `[ meta, file ]`
2. Extracts files only: `[ meta, file ]` → `[ file ]`
3. Creates a **guaranteed flat** channel of individual files

---

## 📋 Collection Patterns Summary

| Aspect | FastQC | DESeq2 |
|--------|--------|--------|
| **Location** | Inside subworkflow | Main workflow |
| **Collection Point** | Within subworkflow | Main workflow channel |
| **Pre-processing** | `.transpose().map { it[1] }` | None |
| **Structure at Mix** | Individual files | Potentially nested |
| **Need .flatten()?** | ❌ No (already flat) | ✅ Yes (can be nested) |
| **Emits per Process** | Multiple single-file outputs | 4 outputs per process |
| **Total Files** | Varies | 24 (6 processes × 4 files) |

---

## 🎯 Best Practice Lessons

### Pattern 1: Subworkflow with Transpose/Map (FastQC style)
✅ **Best for:** Modular, reusable components

```groovy
workflow MY_SUBWORKFLOW {
    main:
    ch_output = Channel.empty()
    
    PROCESS_A(...)
    ch_output = ch_output.mix(PROCESS_A.out.files)
    
    PROCESS_B(...)
    ch_output = ch_output.mix(PROCESS_B.out.files)
    
    emit:
    files = ch_output.transpose().map { it[1] }  // ← Flatten here!
}

// In main workflow
ch_multiqc_files = ch_multiqc_files.mix(MY_SUBWORKFLOW.out.files)
```

**Advantage:** Flattening happens ONCE in the subworkflow

---

### Pattern 2: Direct Mix with Final Flatten (DESeq2 style)
✅ **Best for:** Quick additions to existing workflow

```groovy
// Mix individual outputs
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_A.out.file1)
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_A.out.file2)
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_B.out.file1)

// Flatten before final collection
MULTIQC(
    ch_multiqc_files.flatten().collect(),  // ← Flatten here!
    ...
)
```

**Advantage:** Simple, no need to refactor into subworkflow

---

### Pattern 3: Flatten Individual Outputs (Alternative)
✅ **Alternative:** Flatten each output as it's mixed

```groovy
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_A.out.file1.flatten())
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_A.out.file2.flatten())
ch_multiqc_files = ch_multiqc_files.mix(PROCESS_B.out.file1.flatten())

MULTIQC(
    ch_multiqc_files.collect(),  // ← Already flat, no need to flatten again
    ...
)
```

**Advantage:** Makes flattening explicit at each step

**Disadvantage:** More verbose, repetitive

---

## 🐛 Why the Bug Wasn't Caught Earlier

### 1. FastQC Worked Fine
The subworkflow pattern with `.transpose().map { it[1] }` was already correct, so no issues appeared.

### 2. Other Tools Collected Within Processes
Many tools (STAR, RSEM, RSeQC, etc.) have their outputs collected with `.collect{it[1]}` immediately:

```groovy
ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt.collect{it[1]})
```

This collects and extracts files immediately, preventing nesting.

### 3. DESeq2 Was Different
DESeq2 QC was added WITHOUT the `.collect{it[1]}` pattern:

```groovy
// No collection at mix point!
ch_multiqc_files = ch_multiqc_files.mix(
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO.out.sample_distances_txt
)
```

This left the possibility of nested structures, which was fine... until we removed the **premature** `.collect()` that was happening later.

---

## 🔧 The Complete Fix Timeline

### Original Code (BROKEN)
```groovy
// Files mixed without collection
ch_multiqc_files = ch_multiqc_files.mix(PROCESS.out.files)

// Premature collection (WRONG PLACE)
ch_multiqc_files = ch_multiqc_files.collect()

// Later...
MULTIQC(ch_multiqc_files, ...)
```
**Problem:** Premature `.collect()` created wrong structure

---

### After Fix #1 (STILL BROKEN)
```groovy
// Files mixed without collection
ch_multiqc_files = ch_multiqc_files.mix(PROCESS.out.files)

// NO collection here (removed premature .collect())

// Later...
MULTIQC(ch_multiqc_files.collect(), ...)
```
**Problem:** Nested structures not flattened

---

### After Fix #2 (WORKING!) ✅
```groovy
// Files mixed without collection
ch_multiqc_files = ch_multiqc_files.mix(PROCESS.out.files)

// NO collection here

// Later...
MULTIQC(
    ch_multiqc_files.flatten().collect(),  // ← Added .flatten()!
    ...
)
```
**Solution:** Flatten ensures all items are individual files before collection

---

## 📊 Verification

### Check FastQC Files
```bash
# FastQC files should work as before
find results/fastqc -name "*.zip" -o -name "*.html"
```

### Check DESeq2 Files
```bash
# DESeq2 files should now appear in MultiQC
find results -name "*.deseq2.*.txt" -type f

# Check MultiQC found them
grep "deseq2-kallisto" results/multiqc/multiqc_data/multiqc.log
```

### Expected MultiQC Output
Both FastQC and DESeq2 sections should appear with all plots:
- ✅ FastQC: Per-sample quality reports
- ✅ DESeq2: 8 plots per quantifier (distances, PCA, read distributions)

---

## 💡 Key Takeaways

1. **Different collection patterns serve different purposes**
   - Subworkflows: Flatten outputs before emit
   - Direct mixing: Flatten before final collection

2. **`.transpose().map { it[1] }` is powerful**
   - Extracts files from `[meta, file]` tuples
   - Ensures flat channel structure
   - Used in well-designed subworkflows

3. **`.flatten()` is essential when**
   - Multiple outputs per process
   - Direct mixing into main channel
   - No intermediate flattening step

4. **Both patterns are valid!**
   - FastQC style: Better for reusable components
   - DESeq2 style: Simpler for one-off additions
   - Just need `.flatten()` in the right place!

---

## 🎉 Conclusion

The difference between FastQC and DESeq2 collection patterns reveals important Nextflow channel handling principles:

- **FastQC succeeds** because it flattens outputs in the subworkflow
- **DESeq2 needed help** because it mixed outputs directly without flattening
- **Both patterns work** with the right approach:
  - Subworkflow: Flatten at emit
  - Direct mix: Flatten before final collect

The fix (adding `.flatten()`) makes the DESeq2 pattern consistent with the FastQC approach's end result: a flat list of individual files ready for MultiQC! 🚀

---

**All fixes committed to: https://github.com/pdichiaro/rnaseq**
