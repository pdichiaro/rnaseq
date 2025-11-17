# Scaling Factors Output - Complete Diagram

## 🎯 Quick Answer

**YES!** The pipeline creates `scaling_dat.txt` as a consolidated annotation file containing all scaling factors for all samples.

---

## 📦 Complete Output Structure

```
results/deseq2/
│
├── salmon/                                    (Quantification method)
│   │
│   ├── 📄 scaling_dat.txt                     ← CONSOLIDATED FILE (all samples)
│   │   ├─ Format: TSV (tab-separated)
│   │   ├─ Columns: sample, size_factor
│   │   └─ Purpose: Human-readable reference
│   │
│   ├── 📁 scaling_factors/                    ← INDIVIDUAL FILES (per sample)
│   │   ├── sample1_scaling_factor.txt         (single value: "0.8234")
│   │   ├── sample2_scaling_factor.txt         (single value: "1.0456")
│   │   ├── sample3_scaling_factor.txt         (single value: "0.9123")
│   │   └── ...
│   │   └─ Purpose: Used by DEEPTOOLS_BIGWIG_NORM
│   │
│   ├── Quality_Control/
│   │   ├── PCA plots
│   │   ├── Heatmaps
│   │   └── ...
│   │
│   ├── Read_Distribution/
│   │   └── Distribution plots
│   │
│   └── invariant_genes/                       (Optional: if enabled)
│       │
│       ├── 📄 scaling_dat.txt                 ← CONSOLIDATED FILE (invariant genes)
│       │
│       ├── 📁 scaling_factors/                ← INDIVIDUAL FILES (invariant genes)
│       │   ├── sample1_scaling_factor.txt
│       │   ├── sample2_scaling_factor.txt
│       │   └── ...
│       │
│       └── normalization/
│           ├── invariant_genes_list.txt
│           ├── gene_variances.txt
│           └── ...
│
└── star_rsem/                                 (If STAR+RSEM alignment used)
    ├── 📄 scaling_dat.txt
    ├── 📁 scaling_factors/
    │   └── ...
    └── invariant_genes/
        ├── 📄 scaling_dat.txt
        └── 📁 scaling_factors/
            └── ...
```

---

## 📋 File Format Examples

### 1. Consolidated File: `scaling_dat.txt`

**Location**: `results/deseq2/salmon/scaling_dat.txt`

**Format**: Tab-separated values (TSV)

**Content**:
```tsv
sample	size_factor
HBR_Rep1	0.8234567
HBR_Rep2	0.9876543
HBR_Rep3	0.7654321
UHR_Rep1	1.2345678
UHR_Rep2	1.1234567
UHR_Rep3	1.3456789
```

**Purpose**:
- ✅ Single file with ALL samples
- ✅ Easy to view and compare
- ✅ Import into Excel, R, Python
- ✅ Quality control reference
- ✅ Documentation of normalization

**How to use**:
```bash
# View in terminal
cat results/deseq2/salmon/scaling_dat.txt

# View as table
column -t -s $'\t' results/deseq2/salmon/scaling_dat.txt

# Find outliers (sort by size factor)
sort -t$'\t' -k2 -n results/deseq2/salmon/scaling_dat.txt
```

---

### 2. Individual Files: `scaling_factors/*.txt`

**Location**: `results/deseq2/salmon/scaling_factors/`

**Files Created**:
```
HBR_Rep1_scaling_factor.txt
HBR_Rep2_scaling_factor.txt
HBR_Rep3_scaling_factor.txt
UHR_Rep1_scaling_factor.txt
UHR_Rep2_scaling_factor.txt
UHR_Rep3_scaling_factor.txt
```

**Content of each file** (single numeric value):

`HBR_Rep1_scaling_factor.txt`:
```
0.8234567
```

`UHR_Rep1_scaling_factor.txt`:
```
1.2345678
```

**Purpose**:
- ✅ One file per sample
- ✅ Easy to parse programmatically
- ✅ Used by pipeline processes (DEEPTOOLS_BIGWIG_NORM)
- ✅ Robust design (no parsing required)

**How pipeline uses them**:
```groovy
// Read individual files in channel operation
ch_scaling_per_sample = ch_scaling_factors_individual
    .flatten()
    .map { file ->
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        def scaling_value = file.text.trim()  // Read the single value
        [sample_name, scaling_value]
    }
    // Result: ['HBR_Rep1', '0.8234567']
```

---

## 🔄 Data Flow: From DESeq2 to BigWig

```
┌─────────────────────────────────────────────────────────────────┐
│  Step 1: DESeq2 Normalization                                   │
│  ─────────────────────────────────────────────────────────────  │
│                                                                  │
│  Input: Count matrix (all samples)                              │
│  Process: NORMALIZE_DESEQ2_QC_ALL_GENES                         │
│  Method: DESeq2 estimateSizeFactors()                           │
│                                                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 2: Generate Output Files                                  │
│  ─────────────────────────────────────────────────────────────  │
│                                                                  │
│  TWO types of outputs created simultaneously:                   │
│                                                                  │
│  1️⃣  CONSOLIDATED FILE                                          │
│      └─ scaling_dat.txt (all samples in one table)             │
│                                                                  │
│  2️⃣  INDIVIDUAL FILES                                           │
│      └─ scaling_factors/sample_name_scaling_factor.txt         │
│                                                                  │
└───────────────────┬─────────────────────┬───────────────────────┘
                    │                     │
                    │                     │
         ┌──────────▼──────────┐  ┌───────▼──────────┐
         │  For Humans         │  │  For Pipeline    │
         │  ─────────────────  │  │  ──────────────  │
         │  • View in terminal │  │  • Read in       │
         │  • Import to R/Excel│  │    channels      │
         │  • Quality control  │  │  • Match with    │
         │  • Documentation    │  │    samples       │
         │  • Manual analysis  │  │  • Pass to       │
         │                     │  │    processes     │
         └─────────────────────┘  └────────┬─────────┘
                                           │
                                           ▼
                              ┌────────────────────────┐
                              │  Step 3: Channel Ops   │
                              │  ──────────────────    │
                              │  .flatten()            │
                              │  .map { extract }      │
                              │  .combine(BAM)         │
                              │  .filter { match }     │
                              └────────┬───────────────┘
                                       │
                                       ▼
                              ┌────────────────────────┐
                              │  Step 4: Process       │
                              │  ──────────────────    │
                              │  DEEPTOOLS_BIGWIG_NORM │
                              │  Input: meta, bam,     │
                              │         bai, scaling   │
                              └────────┬───────────────┘
                                       │
                                       ▼
                              ┌────────────────────────┐
                              │  Step 5: Output        │
                              │  ──────────────────    │
                              │  Normalized BigWig     │
                              │  files                 │
                              └────────────────────────┘
```

---

## 🔍 Example: 6-Sample Experiment

### Experiment Setup
- **3 HBR samples**: High Brain RNA (replicates 1-3)
- **3 UHR samples**: Universal Human Reference (replicates 1-3)
- **Quantification**: Salmon
- **Normalization**: Both methods (all_genes + invariant_genes)

---

### Output Directory Structure

```
results/deseq2/salmon/
│
├── 📄 scaling_dat.txt                          ← All genes method (consolidated)
│   ┌──────────────────────────────────────┐
│   │ sample       size_factor             │
│   │ HBR_Rep1     0.8234567               │
│   │ HBR_Rep2     0.9876543               │
│   │ HBR_Rep3     0.7654321               │
│   │ UHR_Rep1     1.2345678               │
│   │ UHR_Rep2     1.1234567               │
│   │ UHR_Rep3     1.3456789               │
│   └──────────────────────────────────────┘
│
├── 📁 scaling_factors/                         ← All genes method (individual)
│   ├── HBR_Rep1_scaling_factor.txt  → "0.8234567"
│   ├── HBR_Rep2_scaling_factor.txt  → "0.9876543"
│   ├── HBR_Rep3_scaling_factor.txt  → "0.7654321"
│   ├── UHR_Rep1_scaling_factor.txt  → "1.2345678"
│   ├── UHR_Rep2_scaling_factor.txt  → "1.1234567"
│   └── UHR_Rep3_scaling_factor.txt  → "1.3456789"
│
├── HBR_Rep1_normalized_counts.txt
├── HBR_Rep2_normalized_counts.txt
├── ... (other QC files)
│
└── invariant_genes/
    │
    ├── 📄 scaling_dat.txt                      ← Invariant genes method (consolidated)
    │   ┌──────────────────────────────────────┐
    │   │ sample       size_factor             │
    │   │ HBR_Rep1     0.7821234               │
    │   │ HBR_Rep2     0.9543210               │
    │   │ HBR_Rep3     0.7123456               │
    │   │ UHR_Rep1     1.3123456               │
    │   │ UHR_Rep2     1.1876543               │
    │   │ UHR_Rep3     1.4234567               │
    │   └──────────────────────────────────────┘
    │
    ├── 📁 scaling_factors/                     ← Invariant genes method (individual)
    │   ├── HBR_Rep1_scaling_factor.txt  → "0.7821234"
    │   ├── HBR_Rep2_scaling_factor.txt  → "0.9543210"
    │   ├── HBR_Rep3_scaling_factor.txt  → "0.7123456"
    │   ├── UHR_Rep1_scaling_factor.txt  → "1.3123456"
    │   ├── UHR_Rep2_scaling_factor.txt  → "1.1876543"
    │   └── UHR_Rep3_scaling_factor.txt  → "1.4234567"
    │
    └── normalization/
        ├── invariant_genes_list.txt
        └── gene_variances.txt
```

---

### How Each File Type Is Used

#### 📄 Consolidated File (`scaling_dat.txt`)

**View all factors at once**:
```bash
$ cat results/deseq2/salmon/scaling_dat.txt
sample      size_factor
HBR_Rep1    0.8234567
HBR_Rep2    0.9876543
HBR_Rep3    0.7654321
UHR_Rep1    1.2345678
UHR_Rep2    1.1234567
UHR_Rep3    1.3456789
```

**Import into R for visualization**:
```r
scaling <- read.table("results/deseq2/salmon/scaling_dat.txt", 
                      header=TRUE, sep="\t")

# Plot size factors
barplot(scaling$size_factor, 
        names.arg=scaling$sample, 
        las=2, 
        ylab="Size Factor",
        main="DESeq2 Size Factors - All Genes")
abline(h=1, col="red", lty=2)  # Reference line at 1.0
```

**Compare normalization methods**:
```bash
# All genes method
cat results/deseq2/salmon/scaling_dat.txt

# Invariant genes method
cat results/deseq2/salmon/invariant_genes/scaling_dat.txt

# Side-by-side comparison
paste results/deseq2/salmon/scaling_dat.txt \
      results/deseq2/salmon/invariant_genes/scaling_dat.txt
```

---

#### 📁 Individual Files (`scaling_factors/*.txt`)

**Used by pipeline**:
```groovy
// Workflow reads these individual files
ch_scaling_factors_individual = NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors_individual

// Extract sample names and values
ch_scaling_per_sample = ch_scaling_factors_individual
    .flatten()  // Separate each file
    .map { file ->
        // Sample name from filename
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        // Read the single value in the file
        def scaling_value = file.text.trim()
        // Emit tuple: [sample_name, scaling_value]
        [sample_name, scaling_value]
    }

// Result channel:
// ['HBR_Rep1', '0.8234567']
// ['HBR_Rep2', '0.9876543']
// ['HBR_Rep3', '0.7654321']
// ['UHR_Rep1', '1.2345678']
// ['UHR_Rep2', '1.1234567']
// ['UHR_Rep3', '1.3456789']
```

**Match with BAM files**:
```groovy
ch_combined = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample)
    .map { meta, bam, bai, sample_id, scaling -> 
        meta.id == sample_id ? [meta, bam, bai, scaling] : null
    }
    .filter { it != null }

// Result for HBR_Rep1:
// [[id:'HBR_Rep1', ...], HBR_Rep1.bam, HBR_Rep1.bam.bai, '0.8234567']
```

**Send to process**:
```groovy
DEEPTOOLS_BIGWIG_NORM (
    ch_combined
)

// Process receives:
// tuple val(meta), path(bam), path(bai), val(scaling)
// 
// Can use directly in script:
// bamCoverage --scaleFactor $scaling ...
```

---

## 🆚 Comparison: Two File Types

| Aspect | scaling_dat.txt | Individual files |
|--------|----------------|------------------|
| **Format** | Tab-separated table | Single value per file |
| **Content** | All samples in one file | One file per sample |
| **Purpose** | Human reference, QC | Pipeline processing |
| **Used by** | Manual inspection, R/Python | DEEPTOOLS_BIGWIG_NORM |
| **Advantages** | Easy comparison, complete view | Robust parsing, clean channels |
| **When to use** | Quality control, documentation | Automated workflows |

---

## ✅ Summary

**Two complementary outputs**:

1. **`scaling_dat.txt`** (consolidated)
   - ✅ All samples in one place
   - ✅ Tab-separated format
   - ✅ Human-readable
   - ✅ Perfect for QC and analysis

2. **Individual `*_scaling_factor.txt` files**
   - ✅ One file per sample
   - ✅ Single value (easy to parse)
   - ✅ Used by pipeline
   - ✅ Robust design

**Both are created automatically** by the DESeq2 normalization processes, providing flexibility for different use cases!

---

## 📚 See Also

- **SCALING_DAT_FILE_DOCUMENTATION.md**: Complete technical documentation (456 lines)
- **HOW_DEEPTOOLS_BIGWIG_NORM_WORKS.md**: How scaling factors are used in normalization
- **IMPLEMENTATION_SUMMARY.md**: Overall implementation details
