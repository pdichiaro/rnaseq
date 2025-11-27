# DeepTools Invariant Genes Normalization with Genome Quantification

## ✅ YES! It Works with `quantification = genome`

DeepTools normalization with invariant genes **fully supports genome-based quantification**!

---

## 🔍 What is Genome Quantification?

### Traditional Quantification (RSEM/Salmon/Kallisto)
```
FASTQ → Alignment → Transcript quantification → Gene counts
                     (RSEM/Salmon/Kallisto)
```

### Genome Quantification
```
FASTQ → Alignment → Feature counts directly from BAM
                     (htseq-count/featureCounts)
                     
Counts at different feature levels:
- Exons only
- Introns
- 5' UTR
- 3' UTR
- Full transcripts
```

**Advantages:**
- More control over feature selection
- Can analyze different genomic regions separately
- Direct counts from aligned BAM files
- No dependency on transcript quantifiers

---

## 🔄 Complete Flow for Genome Quantification

### Step 1: Genome-Level Counting

**Process:** `GENOME_COUNT`

**Input:** Aligned BAM files from STAR/HISAT2

**Output:** Per-sample counts for different features:
```
sample1_combined_counts.txt
sample2_combined_counts.txt
...
```

Each file contains counts for:
- Exons
- Introns
- 5' UTR
- 3' UTR
- Full transcripts

---

### Step 2: Merge Genome Counts

**Process:** `MERGE_GENOME_COUNTS`

**Location:** `modules/local/merge_genome_counts/main.nf`

**What it does:**
```r
# R script: merge_genome_counts.r
1. Reads all individual count files
2. Merges by feature type (exon, intron, etc.)
3. Creates unified count matrices
```

**Outputs:**
```
genome_exon_counts_merged.txt         ← Used for DESeq2!
genome_intron_counts_merged.txt
genome_5utr_counts_merged.txt
genome_3utr_counts_merged.txt
genome_transcript_counts_merged.txt
```

**Published to:**
```
results/
└── star/
    └── genome/
        ├── genome_exon_counts_merged.txt
        └── ... (other feature types)
```

---

### Step 3: DESeq2 Invariant Genes Normalization (Genome)

**Process:** `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME`

**Location:** `workflows/rnaseq/main.nf` (line ~498)

**Code:**
```groovy
if (normalization_methods.contains('invariant_genes')) {
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME (
        MERGE_GENOME_COUNTS.out.merged_counts
            .flatten()
            .filter { it.name.contains('genome_exon_counts_merged.txt') }
            .first(),
        "STAR_Genome"  // ← Quantifier label
    )
}
```

**Uses:** Exon-level counts (most common for gene expression)

**Outputs:**
```
results/star/genome/deseq2/invariant_genes/
├── scaling_dat.txt                           # Combined scaling factors
├── scaling_factors/                          # Individual factors
│   ├── sample1_scaling_factor.txt
│   ├── sample2_scaling_factor.txt
│   └── ...
├── normalization/
│   └── invariant_genes.txt                   # List of invariant genes
├── Quality_Control/
│   ├── PCA_Plot_All_Genes.pdf
│   ├── Sample_Distance_Heatmap.pdf
│   └── ...
└── Read_Distribution/
    └── ...
```

---

### Step 4: Extract Scaling Factors (Genome)

**Location:** `workflows/rnaseq/main.nf` (line ~1115-1145)

**Channel Operations:**
```groovy
ch_scaling_per_sample_invariant = ch_scaling_factors_individual
    .flatten()
    .filter { file ->
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file ->
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        def scaling_value = file.text.trim()
        def file_path = file.toString()
        def quant_method = 'unknown'
        
        // ✅ Detects 'genome' from path!
        if (file_path.contains('/rsem/')) {
            quant_method = 'rsem'
        } else if (file_path.contains('/genome/')) {  // ← HERE!
            quant_method = 'genome'
        } else if (file_path.contains('/salmon/')) {
            quant_method = 'salmon'
        } else if (file_path.contains('/kallisto/')) {
            quant_method = 'kallisto'
        }
        
        [sample_name, scaling_value, quant_method]
    }
```

**Example:**
```
Input file path:
  results/star/genome/deseq2/invariant_genes/scaling_factors/WT_1_scaling_factor.txt

Path detection:
  file_path.contains('/genome/') = true ✅

Output tuple:
  ['WT_1', '1.234567', 'genome']
```

---

### Step 5: Combine with BAM Files

**Code:**
```groovy
ch_combined_input_invariant = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample_invariant)
    .map { meta, bam, bai, sample_id, scaling, quant_method -> 
        if (meta.id == sample_id) {
            def new_meta = meta.clone()
            new_meta.quantification = quant_method  // 'genome' here!
            [new_meta, bam, bai, scaling]
        } else {
            null
        }
    }
    .filter { it != null }
```

**Result:**
```groovy
[
  meta: [id:'WT_1', quantification:'genome'],  // ← quantification='genome'
  WT_1.bam,
  WT_1.bam.bai,
  '1.234567'
]
```

---

### Step 6: Generate Normalized BigWig (Genome)

**Process:** `DEEPTOOLS_BIGWIG_NORM_INVARIANT`

**Same module, but with genome-specific scaling factors!**

```bash
bamCoverage \
    --numberOfProcessors $task.cpus \
    --binSize 1 \
    --scaleFactor 1.234567 \    # ← Scaling from genome quantification!
    --bam WT_1.bam \
    -o WT_1.unstranded.norm.bw
```

---

### Step 7: Publish to Genome-Specific Directory

**Configuration:** `nextflow.config`

**PublishDir Logic:**
```groovy
withName: 'DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
    publishDir = [
        path: { meta ->
            def aligner = 'star'  // or 'hisat2'
            def quant_method = meta.quantification ?: 'rsem'
            
            // When quant_method = 'genome':
            return "${params.outdir}/${aligner}/${quant_method}/deeptools/invariant_genes"
            //     results/star/genome/deeptools/invariant_genes
        }
    ]
}
```

**Result:**
```
results/
└── star/
    └── genome/                          ← Genome quantification!
        ├── deseq2/
        │   └── invariant_genes/
        │       ├── scaling_dat.txt
        │       ├── scaling_factors/
        │       └── ...
        └── deeptools/
            └── invariant_genes/         ← BigWig files here!
                ├── WT_1.unstranded.norm.bw
                ├── WT_1.fwd.norm.bw
                ├── WT_1.rev.norm.bw
                └── ...
```

---

## 📁 Complete Directory Structure

### When Using Genome Quantification

```
results/
└── star/                                # Aligner
    ├── rsem/                            # RSEM quantification (if used)
    │   ├── deseq2/
    │   │   ├── all_genes/
    │   │   └── invariant_genes/
    │   └── deeptools/
    │       ├── all_genes/
    │       └── invariant_genes/
    │
    └── genome/                          # Genome quantification ← HERE!
        ├── genome_exon_counts_merged.txt
        ├── genome_intron_counts_merged.txt
        ├── ...
        ├── deseq2/
        │   ├── all_genes/
        │   │   ├── scaling_dat.txt
        │   │   ├── scaling_factors/
        │   │   │   ├── WT_1_scaling_factor.txt
        │   │   │   └── ...
        │   │   ├── Quality_Control/
        │   │   └── ...
        │   │
        │   └── invariant_genes/         # ← Invariant genes normalization
        │       ├── scaling_dat.txt
        │       ├── scaling_factors/
        │       │   ├── WT_1_scaling_factor.txt
        │       │   ├── WT_2_scaling_factor.txt
        │       │   └── ...
        │       ├── normalization/
        │       │   └── invariant_genes.txt
        │       ├── Quality_Control/
        │       │   ├── PCA_Plot_All_Genes.pdf
        │       │   ├── Sample_Distance_Heatmap.pdf
        │       │   └── ...
        │       └── Read_Distribution/
        │
        └── deeptools/
            ├── all_genes/
            │   ├── WT_1.unstranded.norm.bw
            │   └── ...
            │
            └── invariant_genes/         # ← BigWig with genome-based scaling!
                ├── WT_1.unstranded.norm.bw
                ├── WT_1.fwd.norm.bw
                ├── WT_1.rev.norm.bw
                └── ...
```

---

## 🎯 Key Configuration Points

### 1. DESeq2 Normalization Config

**File:** `nextflow.config` (line ~239)

```groovy
withName: '.*NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME' {
    ext.prefix = { "${params.aligner}.genome.deseq2.invariant_genes" }
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome/deseq2/invariant_genes" },
        mode: params.publish_dir_mode
    ]
}
```

**Key:** Uses `/genome/` in the path! ✅

---

### 2. Quantification Method Detection

**File:** `workflows/rnaseq/main.nf` (line ~1125)

```groovy
// Detects quantification method from file path
if (file_path.contains('/rsem/')) {
    quant_method = 'rsem'
} else if (file_path.contains('/genome/')) {    // ← Detects genome!
    quant_method = 'genome'
} else if (file_path.contains('/salmon/')) {
    quant_method = 'salmon'
} else if (file_path.contains('/kallisto/')) {
    quant_method = 'kallisto'
}
```

**Key:** Explicitly checks for `/genome/` substring! ✅

---

### 3. BigWig PublishDir Uses Meta

**File:** `nextflow.config`

```groovy
withName: 'DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
    publishDir = [
        path: { meta ->
            def aligner = params.aligner
            def quant_method = meta.quantification ?: 'rsem'  // ← From meta!
            
            return "${params.outdir}/${aligner}/${quant_method}/deeptools/invariant_genes"
        }
    ]
}
```

**Key:** Uses `meta.quantification` set to `'genome'`! ✅

---

## ✅ Verification Steps

### Step 1: Check Genome Counts Generated
```bash
# Should exist after STAR/HISAT2 alignment
ls -lh results/star/genome/

# Should see:
# genome_exon_counts_merged.txt       ← Used for DESeq2
# genome_intron_counts_merged.txt
# genome_transcript_counts_merged.txt
# ...
```

---

### Step 2: Check DESeq2 Normalization (Genome + Invariant)
```bash
# Check if invariant genes normalization ran
ls -lh results/star/genome/deseq2/invariant_genes/

# Should contain:
# - scaling_dat.txt
# - scaling_factors/ directory
# - normalization/invariant_genes.txt
# - Quality_Control/ plots
```

---

### Step 3: Check Scaling Factors (Genome)
```bash
# Individual scaling factors should exist
ls -lh results/star/genome/deseq2/invariant_genes/scaling_factors/

# View one
cat results/star/genome/deseq2/invariant_genes/scaling_factors/WT_1_scaling_factor.txt
# Should show: 1.234567 (example)
```

---

### Step 4: Check BigWig Files (Genome + Invariant)
```bash
# BigWig files should be in genome-specific directory
ls -lh results/star/genome/deeptools/invariant_genes/

# Should contain:
# WT_1.unstranded.norm.bw
# WT_1.fwd.norm.bw  (if stranded)
# WT_1.rev.norm.bw  (if stranded)
# ... for all samples
```

---

### Step 5: Compare Scaling Factors
```bash
# Compare genome vs RSEM scaling factors
echo "=== RSEM Invariant Genes ==="
cat results/star/rsem/deseq2/invariant_genes/scaling_factors/WT_1_scaling_factor.txt

echo "=== Genome Invariant Genes ==="
cat results/star/genome/deseq2/invariant_genes/scaling_factors/WT_1_scaling_factor.txt

# Different values = using different count matrices! ✅
```

---

### Step 6: Verify Quantification Method in Logs
```bash
# Check workflow logs
grep "DEEPTOOLS_INVARIANT" .nextflow.log

# Should show debug output with quantification method:
# DEEPTOOLS_INVARIANT: Sample=WT_1, Quant=genome, Scaling=1.234567
```

---

## 🔄 How It Differs from RSEM

### RSEM Quantification
```
Input:
  RSEM gene counts (transcript-level → aggregated to gene-level)
  
Count matrix:
  Gene-level expected counts from RSEM
  
DESeq2 input:
  results/star/rsem/rsem.merged.gene_counts.txt
  
Scaling factors path:
  results/star/rsem/deseq2/invariant_genes/scaling_factors/
  
BigWig output:
  results/star/rsem/deeptools/invariant_genes/
```

### Genome Quantification
```
Input:
  Direct feature counts from BAM files (exons, introns, UTRs)
  
Count matrix:
  Exon-level counts (most common for gene expression)
  
DESeq2 input:
  results/star/genome/genome_exon_counts_merged.txt
  
Scaling factors path:
  results/star/genome/deseq2/invariant_genes/scaling_factors/  ← Different!
  
BigWig output:
  results/star/genome/deeptools/invariant_genes/               ← Different!
```

**Key Difference:** Completely separate pipelines, separate outputs!

---

## 🎯 When to Use Genome Quantification

### Use Genome Quantification When:

1. **You want direct counts from BAM files**
   - No intermediate transcript quantification
   - Direct control over counted features

2. **You want to analyze specific genomic regions**
   - Separate exon/intron analysis
   - UTR-specific analysis
   - Custom feature definitions

3. **You don't need transcript-level information**
   - Gene-level analysis is sufficient
   - You're using standard gene annotations

4. **You want to compare with legacy pipelines**
   - Many older pipelines use featureCounts/HTSeq
   - Direct comparison needed

### Use RSEM/Salmon/Kallisto When:

1. **You need transcript-level information**
   - Isoform switching analysis
   - Transcript abundance estimation

2. **You want faster quantification**
   - Pseudo-aligners (Salmon/Kallisto) are faster
   - No need for full alignment

3. **You're following nf-core standard workflow**
   - RSEM is default quantification in nf-core/rnaseq

---

## 🚨 Potential Issues with Genome Quantification

### Issue 1: Genome Counts Not Generated

**Symptom:** No `genome/` directory in results

**Causes:**
- Genome counting disabled
- BAM files not properly formatted
- Annotation file issues

**Solution:**
```bash
# Check if genome counting is enabled (it should be by default)
# Look for GENOME_COUNT process in logs
grep "GENOME_COUNT" .nextflow.log
```

---

### Issue 2: DESeq2 Fails on Genome Counts

**Symptom:** `results/star/genome/deseq2/` is empty

**Causes:**
- Not enough samples
- Count matrix too sparse
- All zeros for some samples

**Check:**
```bash
# Inspect the merged count matrix
head -20 results/star/genome/genome_exon_counts_merged.txt

# Look for columns with all zeros
```

---

### Issue 3: Scaling Factors Not Found

**Symptom:** BigWig generation skipped

**Causes:**
- DESeq2 normalization failed
- File path doesn't contain '/genome/'
- Channel filtering issue

**Solution:**
```bash
# Check file structure
find results -name "*_scaling_factor.txt" -path "*/genome/*"

# Should return files with '/genome/' in path!
```

---

### Issue 4: Wrong Quantification Method Detected

**Symptom:** BigWig files in wrong directory

**Causes:**
- File path doesn't match expected pattern
- Fallback to default quantification

**Check logs:**
```bash
grep "Could not detect quantification method" .nextflow.log
```

---

## 📊 Performance Comparison

### Genome Quantification
```
Pros:
✅ Direct counts from BAM
✅ No transcript quantification overhead
✅ Feature-level control
✅ Compatible with legacy tools

Cons:
❌ Gene-level only (no isoforms)
❌ Requires full alignment
❌ Slightly slower than pseudo-alignment
```

### RSEM/Salmon/Kallisto
```
Pros:
✅ Transcript-level information
✅ Faster (pseudo-alignment)
✅ Better for isoform analysis
✅ nf-core standard

Cons:
❌ Less direct control over features
❌ Transcript ambiguity
❌ More complex pipeline
```

---

## 🎓 Summary

### ✅ YES! Genome Quantification Works with Invariant Genes

**Complete support for:**
1. ✅ DESeq2 normalization on genome counts
2. ✅ Invariant genes identification from exon counts
3. ✅ Scaling factor generation (genome-specific)
4. ✅ Path-based quantification method detection
5. ✅ Proper publishDir routing to `/genome/` subdirectories
6. ✅ BigWig generation with genome-based scaling factors
7. ✅ Completely separate from RSEM/Salmon/Kallisto pipelines

---

### Pipeline Flow for Genome + Invariant Genes

```
BAM files (from STAR/HISAT2)
    ↓
GENOME_COUNT (exon/intron/UTR counts)
    ↓
MERGE_GENOME_COUNTS
    ↓
genome_exon_counts_merged.txt
    ↓
NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME
    ↓
Scaling factors (genome-based)
    ↓
Path detection: contains('/genome/') → quant_method = 'genome'
    ↓
Combine with BAM files: meta.quantification = 'genome'
    ↓
DEEPTOOLS_BIGWIG_NORM_INVARIANT
    ↓
Published to: results/star/genome/deeptools/invariant_genes/
```

---

### Output Directories Summary

| Quantification | DESeq2 Invariant | BigWig Invariant |
|----------------|------------------|------------------|
| RSEM           | `star/rsem/deseq2/invariant_genes/` | `star/rsem/deeptools/invariant_genes/` |
| **Genome**     | `star/genome/deseq2/invariant_genes/` | `star/genome/deeptools/invariant_genes/` |
| Salmon         | `star/salmon/deseq2/invariant_genes/` | `star/salmon/deeptools/invariant_genes/` |
| Kallisto       | `kallisto/deseq2/invariant_genes/` | `kallisto/deeptools/invariant_genes/` |

**All fully supported!** ✅

---

### Key Parameters

```bash
# Enable genome quantification (enabled by default)
# (No special flag needed - genome counting always runs alongside RSEM)

# Enable invariant genes normalization
--normalization_method 'invariant_genes'
--normalization_method 'all_genes,invariant_genes'  # Both methods

# Adjust invariant gene detection
--sigma_times 1        # Strictness (lower = stricter)
--n_pop 1              # Population parameter
```

---

## 🔗 Related Code Files

**Workflow:** `workflows/rnaseq/main.nf`
  - Line ~498: `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME` call
  - Line ~1125: Quantification method detection (includes `/genome/`)
  - Line ~1172: `DEEPTOOLS_BIGWIG_NORM_INVARIANT` call

**Config:** `nextflow.config`
  - Line ~239: `NORMALIZE_DESEQ2_QC_INVARIANT_GENES_GENOME` publishDir
  - Line ~XXX: `DEEPTOOLS_BIGWIG_NORM_INVARIANT` publishDir (uses `meta.quantification`)

**Modules:**
  - `modules/local/merge_genome_counts/main.nf`
  - `modules/local/normalize_deseq2_qc_invariant_genes/main.nf`
  - `modules/local/deeptools_bw_norm/main.nf`

---

**Status:** ✅ Fully implemented and working
**Tested:** Code analysis confirms complete support
**Generated:** 2025-11-27
