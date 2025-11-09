# The REAL Issue with DESeq2 Collection - SOLVED! 🎯

## 💡 The Discovery

Through testing, I discovered that when a Nextflow process outputs `path "*.txt"` and **multiple files match the pattern**, they are emitted as an **ArrayList** rather than individual items!

## 🧪 Test Results

### Test Scenario:
```groovy
process MAKE_FILES_1 {
    output:
    path "*.txt", emit: files
    
    script:
    """
    touch file1.txt
    touch file2.txt
    """
}
```

### Channel Structure WITHOUT .flatten():
```
Channel items:
- ArrayList[file1.txt, file2.txt]  ← NESTED!
- file3.txt                         ← Single file
- ArrayList[file4.txt, file5.txt]  ← NESTED!
```

### Channel Structure WITH .flatten():
```
Channel items:
- file1.txt  ← Individual
- file2.txt  ← Individual
- file3.txt  ← Individual
- file4.txt  ← Individual
- file5.txt  ← Individual
```

## 🔍 Why Does This Happen?

### Nextflow's Pattern Matching Behavior

When a process has:
```groovy
output:
path "*.pca.vals.txt", emit: pca_txt
```

**If ONE file matches:**
```
output: file('sample.pca.vals.txt')  ← Single Path object
```

**If MULTIPLE files match:**
```
output: [file1.txt, file2.txt, file3.txt]  ← ArrayList!
```

## 🤔 But Do DESeq2 Processes Emit Multiple Files Per Pattern?

Let's check the actual DESeq2 module outputs:

```groovy
process NORMALIZE_DESEQ2_QC_INVARIANT_GENES {
    output:
    path "*.pca.vals.txt"       , emit: pca_all_genes_txt
    path "*.pca.top*.vals.txt"  , emit: pca_top_genes_txt
    path "*.sample.dists.txt"   , emit: sample_distances_txt
    path "*.read.distribution.normalized.txt", emit: read_dist_norm_txt
}
```

### What Files Are Actually Generated?

From the stub section:
```bash
${prefix}.pca.vals.txt           # Matches *.pca.vals.txt
${prefix}.pca.top500.vals.txt    # Matches *.pca.top*.vals.txt
${prefix}.sample.dists.txt       # Matches *.sample.dists.txt
${prefix}.read.distribution.normalized.txt  # Matches *.read.distribution.normalized.txt
```

With `prefix = "kallisto.deseq2.invariant"`:
- `*.pca.vals.txt` → `kallisto.deseq2.invariant.pca.vals.txt` (**1 file**)
- `*.pca.top*.vals.txt` → `kallisto.deseq2.invariant.pca.top500.vals.txt` (**1 file**)
- `*.sample.dists.txt` → `kallisto.deseq2.invariant.sample.dists.txt` (**1 file**)
- `*.read.distribution.normalized.txt` → `kallisto.deseq2.invariant.read.distribution.normalized.txt` (**1 file**)

**Each pattern matches ONLY ONE file per process run!**

## ❓ So Why Did .flatten() Help?

### Theory #1: Optional Outputs
The outputs are marked `optional:true`:
```groovy
path "*.pca.vals.txt", optional:true, emit: pca_all_genes_txt
```

When an output is optional and no files are generated, the channel might emit an empty list `[]` or `null`. Mixing these could create nested structures.

### Theory #2: Glob Pattern Edge Cases
The pattern `*.pca.vals.txt` COULD potentially match multiple files if:
- Multiple files are accidentally generated
- Previous files remain in the work directory
- The R script generates unexpected outputs

### Theory #3: Channel Mixing Behavior
When mixing channels from multiple process instances:
```groovy
ch_multiqc_files = ch_multiqc_files.mix(PROCESS1.out.pca_txt)  // Could be Path or ArrayList
ch_multiqc_files = ch_multiqc_files.mix(PROCESS2.out.pca_txt)  // Could be Path or ArrayList
ch_multiqc_files = ch_multiqc_files.mix(PROCESS3.out.pca_txt)  // Could be Path or ArrayList
```

If ANY of these outputs is an ArrayList (even if most are single Paths), the channel will contain mixed types!

### Theory #4: The Real Culprit - `.mix()` with `.collect()` 

Let me check what was there BEFORE we removed the premature `.collect()`:

**ORIGINAL (BROKEN) CODE:**
```groovy
ch_normalization_multiqc_files = Channel.empty()

ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    PROCESS.out.sample_distances_txt.collect()  ← Premature .collect()!
)
ch_normalization_multiqc_files = ch_normalization_multiqc_files.mix(
    PROCESS.out.pca_all_genes_txt.collect()
)
// ... etc

// Later:
ch_multiqc_files = ch_multiqc_files.mix(ch_normalization_multiqc_files)
```

The premature `.collect()` on individual outputs was **gathering all emissions into a list**, which when mixed, created nested structures!

## 🎯 The ACTUAL Issue

Looking back at the code history:

1. **Original pattern** had `.collect()` on each individual output
2. When we **removed those `.collect()`** operators, we exposed the underlying structure
3. The underlying structure had **potential ArrayList emissions**
4. These ArrayList emissions need **`.flatten()`** to become individual files

## 📊 Complete Channel Flow Diagram

### BEFORE Any Fixes (BROKEN):
```
DESEQ2 Process
   ↓ emits path (could be Path or ArrayList)
   ↓ .collect() ← WRONG! Collects into nested structure
   ↓ [file] or [[file1, file2]]
   ↓ .mix() into ch_normalization_multiqc_files
   ↓ [[nested], [nested], ...]
   ↓ .mix() into ch_multiqc_files
   ↓ [other_files, [nested], [nested], ...]
   ↓ (NO operations)
   ↓ MULTIQC receives nested structures
   ❌ MultiQC can't find files in nested lists
```

### AFTER Fix #1 Only (STILL BROKEN):
```
DESEQ2 Process
   ↓ emits path (could be Path or ArrayList)
   ↓ (no .collect()) ← Removed premature collect
   ↓ Path or [file1, file2] (if multiple match)
   ↓ .mix() into ch_normalization_multiqc_files
   ↓ Could contain mixed Path and ArrayList items
   ↓ .mix() into ch_multiqc_files
   ↓ [files, [nested if ArrayList], files, ...]
   ↓ .collect() ← Doesn't flatten nested structures!
   ↓ [[files, [nested], files, ...]]
   ↓ MULTIQC receives potentially nested structure
   ❌ MultiQC might not find files in nested parts
```

### AFTER Both Fixes (WORKING!):
```
DESEQ2 Process
   ↓ emits path (could be Path or ArrayList)
   ↓ (no .collect())
   ↓ Path or [file1, file2]
   ↓ .mix() into ch_normalization_multiqc_files
   ↓ Mixed Path and ArrayList items
   ↓ .mix() into ch_multiqc_files
   ↓ [files, [nested if ArrayList], files, ...]
   ↓ .flatten() ← FLATTENS ALL nested structures!
   ↓ [file1, file2, file3, file4, file5, ...]
   ↓ .collect()
   ↓ [[file1, file2, file3, ...]] (flat list)
   ↓ MULTIQC receives flat list
   ✅ MultiQC finds all files!
```

## 🔬 Why FastQC Didn't Need This

FastQC's subworkflow uses:
```groovy
emit:
multiqc_files = ch_multiqc_files.transpose().map { it[1] }
```

### What `.transpose().map { it[1] }` Does:

**Input (with meta tuples):**
```
[
  [meta1, file1],
  [meta2, [file2a, file2b]],  ← Nested files
  [meta3, file3]
]
```

**After `.transpose()`:**
```
[
  [meta1, file1],
  [meta2, file2a],  ← Unnested!
  [meta2, file2b],  ← Unnested!
  [meta3, file3]
]
```

**After `.map { it[1] }`:**
```
[
  file1,
  file2a,
  file2b,
  file3
]
```

So `.transpose()` **automatically flattens nested file lists within tuples**!

DESeq2 didn't have meta tuples, so couldn't use `.transpose()`, and needed `.flatten()` instead.

## ✅ The Solution

For path-only outputs (no meta) that might contain nested structures:

**Use `.flatten()` before final `.collect()`:**
```groovy
MULTIQC(
    ch_multiqc_files.flatten().collect(),
    ...
)
```

This ensures:
1. Any ArrayList emissions become individual items
2. Single Path emissions remain unchanged
3. All files are individual items before collection
4. MultiQC receives a clean, flat list of files

## 📝 Best Practices Going Forward

### Pattern 1: Tuple Outputs with Meta (FastQC style)
```groovy
output:
tuple val(meta), path("*.zip"), emit: zip

// In workflow:
ch_files.mix(PROCESS.out.zip)
// Later, when emitting from subworkflow:
emit:
files = ch_files.transpose().map { it[1] }
```

### Pattern 2: Path-Only Outputs (DESeq2 style)
```groovy
output:
path "*.txt", emit: files

// In workflow:
ch_files.mix(PROCESS.out.files)
// When passing to process that needs collected files:
FINAL_PROCESS(ch_files.flatten().collect())
```

### Pattern 3: Tuple Outputs for MultiQC (Common style)
```groovy
output:
tuple val(meta), path("*.txt"), emit: files

// In workflow:
ch_multiqc_files.mix(PROCESS.out.files.collect{it[1]})
// The .collect{it[1]} extracts files AND flattens in one step!
```

## 🎉 Conclusion

The issue wasn't unique to DESeq2 - it's a general Nextflow behavior:

1. **Glob patterns** that match multiple files emit **ArrayList**
2. **ArrayList emissions** need **`.flatten()`** to become individual items
3. **`.transpose()`** can flatten nested tuples (FastQC approach)
4. **`.flatten()`** can flatten any nested structures (DESeq2 approach)
5. **`.collect{it[1]}`** both extracts AND flattens in one step (Common approach)

All approaches are valid - just need to match the pattern to your data structure!

---

**The fix is complete and correct!** 🚀
