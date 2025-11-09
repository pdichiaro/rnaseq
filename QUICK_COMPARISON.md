# Quick Visual Comparison: FastQC vs DESeq2 Collection

## 🎯 The Key Difference

```
┌─────────────────────────────────────────────────────────────────┐
│                        FastQC Pattern                           │
└─────────────────────────────────────────────────────────────────┘

Subworkflow:
┌────────────────────────────────────────┐
│  FASTQC processes                      │
│  ↓                                     │
│  .mix() internally                     │
│  ↓                                     │
│  .transpose().map { it[1] }            │  ← Flattening happens HERE
│  ↓                                     │
│  emit: multiqc_files (already flat)    │
└────────────────────────────────────────┘
         ↓
Main Workflow:
ch_multiqc_files.mix(SUBWORKFLOW.out.multiqc_files)  ← Already flat
         ↓
MULTIQC(ch_multiqc_files.collect(), ...)             ← Works perfectly!


┌─────────────────────────────────────────────────────────────────┐
│                        DESeq2 Pattern                           │
└─────────────────────────────────────────────────────────────────┘

Main Workflow:
DESeq2 processes (emit 4 files each)
         ↓
ch_multiqc_files.mix(PROCESS.out.file1)              ← Potentially nested
ch_multiqc_files.mix(PROCESS.out.file2)
ch_multiqc_files.mix(PROCESS.out.file3)
ch_multiqc_files.mix(PROCESS.out.file4)
         ↓                                            ← No flattening yet!
MULTIQC(ch_multiqc_files.flatten().collect(), ...)   ← Flattening happens HERE
```

---

## 📊 Side-by-Side Comparison

| Aspect | FastQC | DESeq2 |
|--------|--------|--------|
| **Where are files mixed?** | Inside subworkflow | Main workflow |
| **When is flattening done?** | At subworkflow emit | Before MultiQC call |
| **Flattening method** | `.transpose().map { it[1] }` | `.flatten()` |
| **Pattern complexity** | Higher (subworkflow) | Lower (direct mix) |
| **Reusability** | ✅ High | ⚠️ Lower |
| **When it works** | Always (built-in) | After adding `.flatten()` |

---

## 🔍 Channel Structure Examples

### FastQC (at MultiQC input)
```groovy
ch_multiqc_files = [
  sample1_R1.fastqc.zip,      ← Individual files
  sample1_R2.fastqc.zip,
  sample2_R1.fastqc.zip,
  sample2_R2.fastqc.zip,
  sample1.fastp.json,
  sample2.fastp.json
]
```
✅ Already flat! `.collect()` works fine.

---

### DESeq2 BEFORE .flatten() ❌
```groovy
ch_multiqc_files = [
  sample1_R1.fastqc.zip,
  sample1_R2.fastqc.zip,
  [                                    ← Nested structure!
    kallisto.deseq2.distances.txt,
    kallisto.deseq2.pca_all.txt,
    kallisto.deseq2.pca_top.txt,
    kallisto.deseq2.read_dist.txt
  ],
  sample2_R1.fastqc.zip,
  [                                    ← Another nested structure!
    salmon.deseq2.distances.txt,
    salmon.deseq2.pca_all.txt,
    ...
  ]
]
```
❌ Nested! MultiQC can't find files in nested lists.

---

### DESeq2 AFTER .flatten() ✅
```groovy
ch_multiqc_files.flatten() = [
  sample1_R1.fastqc.zip,
  sample1_R2.fastqc.zip,
  kallisto.deseq2.distances.txt,      ← Now individual!
  kallisto.deseq2.pca_all.txt,
  kallisto.deseq2.pca_top.txt,
  kallisto.deseq2.read_dist.txt,
  sample2_R1.fastqc.zip,
  salmon.deseq2.distances.txt,        ← Now individual!
  salmon.deseq2.pca_all.txt,
  ...
]
```
✅ Flat! MultiQC finds all files.

---

## 💡 The Lesson

Both patterns are valid, just need the right approach:

### Pattern 1: Flatten in Subworkflow (FastQC)
```groovy
workflow FASTQC_SUBWORKFLOW {
    emit:
    multiqc_files = ch_files.transpose().map { it[1] }
}

// In main workflow
ch_multiqc_files = ch_multiqc_files.mix(FASTQC_SUBWORKFLOW.out.multiqc_files)
MULTIQC(ch_multiqc_files.collect(), ...)  // ← No .flatten() needed
```

### Pattern 2: Flatten at Collection (DESeq2)
```groovy
// Direct mixing
ch_multiqc_files = ch_multiqc_files.mix(PROCESS.out.files)

// Flatten before collection
MULTIQC(ch_multiqc_files.flatten().collect(), ...)  // ← .flatten() needed!
```

---

## 🎯 When to Use Which Pattern?

### Use Subworkflow Pattern (FastQC style) when:
- ✅ Building reusable components
- ✅ Multiple related processes
- ✅ Want to encapsulate logic
- ✅ Planning to share/publish

### Use Direct Mix Pattern (DESeq2 style) when:
- ✅ Quick additions to existing workflow
- ✅ One-off custom processes
- ✅ Don't need reusability
- ✅ Want simplicity

**Just remember:** Add `.flatten()` if you're mixing directly!

---

## 🚀 Result

Both FastQC and DESeq2 now work perfectly in MultiQC! 🎉

```
MultiQC Report:
├── FastQC
│   ├── Per base sequence quality
│   ├── Per sequence quality scores
│   └── ...
└── DESeq2 Kallisto QC
    ├── Sample Distance (All Genes)
    ├── Sample Distance (Invariant Genes)
    ├── PCA All Genes
    ├── PCA Top 500
    └── Read Distribution
```
