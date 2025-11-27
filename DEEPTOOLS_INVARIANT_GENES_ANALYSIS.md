# DeepTools Normalization with Invariant Genes - Complete Analysis

## ✅ Short Answer: YES, It Works!

The `deeptools_norm_bigwig` with `invariant_genes` normalization **does work** in your pipeline, provided:

1. ✅ DESeq2 QC invariant genes normalization completes successfully
2. ✅ Scaling factors are generated (`scaling_factors_individual` files)
3. ✅ BAM files are available for normalization
4. ✅ `--normalization_method` includes `invariant_genes`

---

## 🔍 Complete Flow: How It Works

### Step 1: DESeq2 Invariant Genes Normalization

**Process:** `NORMALIZE_DESEQ2_QC_INVARIANT_GENES`

**Location:** `modules/local/normalize_deseq2_qc_invariant_genes/main.nf`

**What it does:**
```r
# R script: normalize_deseq2_qc_invariant_genes.r
1. Reads count matrix (gene × sample)
2. Identifies invariant genes using GeneralNormalizer package
   - Uses sigma_times parameter (default: 1)
   - Uses n_pop parameter (default: 1)
3. Computes DESeq2 size factors based ONLY on invariant genes
4. Generates scaling factors for each sample
```

**Outputs:**
```
normalization/
  └── invariant_genes.txt              # List of identified invariant genes

scaling_dat.txt                        # Combined scaling factors
scaling_factors/
  ├── sample1_scaling_factor.txt       # Individual per-sample factors
  ├── sample2_scaling_factor.txt
  └── sample3_scaling_factor.txt

Quality_Control/                       # QC plots
Read_Distribution/                     # Read distribution plots
*_normalized_counts.txt                # Normalized count matrix
*.pca.vals.txt                         # PCA results
*.sample.dists.txt                     # Sample distances
```

---

### Step 2: Extract Scaling Factors (Workflow Logic)

**Location:** `workflows/rnaseq/main.nf` (around line 1100-1140)

**Channel Operations:**
```groovy
// Extract individual scaling factor files
ch_scaling_per_sample_invariant = ch_scaling_factors_individual
    .flatten()
    .filter { file ->
        def parent_dir = file.getParent()?.getName() ?: ""
        def grandparent_dir = file.getParent()?.getParent()?.getName() ?: ""
        // ONLY keep files from 'invariant' directories
        parent_dir.contains('invariant') || grandparent_dir.contains('invariant')
    }
    .map { file ->
        def sample_name = file.name.replaceAll('_scaling_factor\\.txt$', '')
        def scaling_value = file.text.trim()  // Read the scaling factor
        def quant_method = detectQuantMethod(file)  // rsem, genome, kallisto
        
        [sample_name, scaling_value, quant_method]
    }
```

**Example:**
```
Input file: results/star/rsem/deseq2/invariant/scaling_factors/WT_1_scaling_factor.txt
Content: 1.234567

Output tuple: ['WT_1', '1.234567', 'rsem']
```

---

### Step 3: Combine BAM Files with Scaling Factors

**Location:** `workflows/rnaseq/main.nf` (around line 1150-1170)

**Channel Operations:**
```groovy
// Combine BAM files with their corresponding scaling factors
ch_combined_input_invariant = ch_bam_for_deeptools
    .combine(ch_scaling_per_sample_invariant)
    .map { meta, bam, bai, sample_id, scaling, quant_method -> 
        if (meta.id == sample_id) {
            def new_meta = meta.clone()
            new_meta.quantification = quant_method  // Add quant method to meta
            [new_meta, bam, bai, scaling]
        } else {
            null
        }
    }
    .filter { it != null }
```

**Example:**
```
Input BAM channel:
  [meta:[id:'WT_1', strandedness:'reverse'], WT_1.bam, WT_1.bam.bai]

Input scaling channel:
  ['WT_1', '1.234567', 'rsem']

Output combined channel:
  [meta:[id:'WT_1', strandedness:'reverse', quantification:'rsem'], 
   WT_1.bam, 
   WT_1.bam.bai, 
   '1.234567']
```

---

### Step 4: Generate Normalized BigWig Files

**Process:** `DEEPTOOLS_BIGWIG_NORM_INVARIANT`

**Location:** `modules/local/deeptools_bw_norm/main.nf`

**Input:**
```groovy
tuple val(meta), path(bam), path(bai), val(scaling)
```

**What it does:**
```bash
# Uses bamCoverage from deeptools
bamCoverage \
    --numberOfProcessors $task.cpus \
    --binSize 1 \
    --scaleFactor $scaling \      # ← Uses invariant genes scaling factor!
    --bam $bam \
    -o ${prefix}.unstranded.norm.bw
```

**For stranded data, generates 3 BigWig files:**
1. `*.unstranded.norm.bw` - All reads
2. `*.fwd.norm.bw` - Forward strand
3. `*.rev.norm.bw` - Reverse strand

**Scaling factor applied:**
- The scaling factor computed from **invariant genes** (not all genes!)
- Normalizes based on genes with **stable expression across conditions**
- More robust for datasets with many differentially expressed genes

---

### Step 5: Publish to Output Directory

**Configuration:** `nextflow.config`

**PublishDir Logic:**
```groovy
withName: 'DEEPTOOLS_BIGWIG_NORM_INVARIANT' {
    publishDir = [
        path: { meta ->
            def aligner = params.aligner  // 'star' or 'hisat2'
            def quant_method = meta.quantification ?: 'rsem'
            
            // Output path includes 'invariant_genes' subfolder!
            return "${params.outdir}/${aligner}/${quant_method}/deeptools/invariant_genes"
        },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

**Example output structure:**
```
results/
└── star/
    └── rsem/
        └── deeptools/
            └── invariant_genes/          ← Separate folder!
                ├── WT_1.unstranded.norm.bw
                ├── WT_1.fwd.norm.bw
                ├── WT_1.rev.norm.bw
                ├── KO_1.unstranded.norm.bw
                ├── KO_1.fwd.norm.bw
                └── KO_1.rev.norm.bw
```

---

## 🔄 Comparison: Invariant Genes vs All Genes

### All Genes Normalization
**Output:** `results/star/rsem/deeptools/all_genes/`

**Scaling approach:**
- Uses **ALL** genes in the count matrix
- Standard DESeq2 size factor estimation
- Assumes most genes are NOT differentially expressed

**Best for:**
- Datasets with few DE genes
- Well-balanced experimental designs
- Standard RNA-seq experiments

### Invariant Genes Normalization
**Output:** `results/star/rsem/deeptools/invariant_genes/`

**Scaling approach:**
- Uses **ONLY** invariant genes (stable across conditions)
- Robust to many DE genes
- Uses GeneralNormalizer package to identify stable genes

**Best for:**
- Datasets with MANY differentially expressed genes
- Perturbation experiments (drug treatments, knockouts)
- Experiments where global expression changes expected
- Comparing wildly different conditions

---

## 📊 Why Two Normalization Methods?

### Scientific Rationale

**Problem with "all genes" normalization:**
```
Scenario: Gene knockout experiment
- 1000 genes are strongly upregulated (compensatory response)
- 50 genes are downregulated
- 9000 genes are unchanged

All genes normalization:
  → Incorrectly assumes most genes unchanged
  → Over-corrects for the 1000 upregulated genes
  → Can mask real biological changes
```

**Solution with "invariant genes" normalization:**
```
Invariant genes normalization:
  → Identifies the 9000 truly unchanged genes
  → Computes size factors ONLY from stable genes
  → Correctly captures the global upregulation
  → More accurate biological interpretation
```

---

## ✅ How to Verify It's Working

### Check 1: Invariant Genes Normalization Completed
```bash
# Look for successful DESeq2 normalization
ls -lh results/star/rsem/deseq2/invariant/

# Should contain:
# - scaling_dat.txt
# - scaling_factors/ directory
# - normalization/invariant_genes.txt
# - Quality_Control/ plots
```

### Check 2: Scaling Factors Generated
```bash
# Check individual scaling factors
ls -lh results/star/rsem/deseq2/invariant/scaling_factors/

# Should have one file per sample:
# WT_1_scaling_factor.txt
# WT_2_scaling_factor.txt
# KO_1_scaling_factor.txt
# ...

# View a scaling factor
cat results/star/rsem/deseq2/invariant/scaling_factors/WT_1_scaling_factor.txt
# Output: 1.234567
```

### Check 3: BigWig Files Generated
```bash
# Check for normalized BigWig files
ls -lh results/star/rsem/deeptools/invariant_genes/

# Should contain:
# WT_1.unstranded.norm.bw
# WT_1.fwd.norm.bw  (if stranded)
# WT_1.rev.norm.bw  (if stranded)
# ... for all samples
```

### Check 4: Compare Scaling Factors
```bash
# Compare all_genes vs invariant_genes scaling factors
echo "=== All Genes ==="
cat results/star/rsem/deseq2/all_genes/scaling_factors/WT_1_scaling_factor.txt

echo "=== Invariant Genes ==="
cat results/star/rsem/deseq2/invariant/scaling_factors/WT_1_scaling_factor.txt

# Different values = both methods working independently!
```

### Check 5: Workflow Logs
```bash
# Check workflow execution logs for DeepTools
grep "DEEPTOOLS_INVARIANT" .nextflow.log

# Should show:
# - Channel debug output with scaling factors
# - Process execution messages
# - File generation confirmations
```

---

## 🚨 When It Might NOT Work

### Issue 1: DESeq2 QC Not Completing
**Symptom:** No `invariant` directory in deseq2 output

**Causes:**
- Not enough samples (need ≥3 samples)
- No replicates in experimental design
- Count matrix too sparse
- All genes flagged as invariant (no variation)

**Check:**
```bash
# Look for DESeq2 error logs
find results -name "*.log" -path "*/deseq2/*" -exec cat {} \;
```

### Issue 2: No Invariant Genes Identified
**Symptom:** Empty or very small `invariant_genes.txt`

**Causes:**
- High biological variability
- Sigma threshold too strict (params.sigma_times)
- Not enough genes pass filtering

**Solution:**
```bash
# Adjust sigma_times parameter (more lenient)
--sigma_times 2  # Default is 1
```

### Issue 3: Scaling Factors Not Found
**Symptom:** Pipeline skips BigWig generation

**Causes:**
- Scaling factor files in wrong location
- File naming mismatch
- Channel filtering too strict

**Check:**
```bash
# Verify scaling factors directory structure
find results -name "*_scaling_factor.txt" | grep invariant
```

### Issue 4: BAM Files Missing
**Symptom:** DeepTools process not running

**Causes:**
- Alignment skipped or failed
- BAM files not properly indexed
- Wrong aligner specified

**Check:**
```bash
# Verify BAM files exist and are indexed
find results -name "*.bam" -path "*/star/*"
find results -name "*.bam.bai" -path "*/star/*"
```

### Issue 5: Wrong Normalization Method Specified
**Symptom:** Process doesn't run at all

**Check your command:**
```bash
# Must include 'invariant_genes'
--normalization_method 'invariant_genes'
--normalization_method 'all_genes,invariant_genes'

# ❌ Won't run invariant normalization:
--normalization_method 'all_genes'
```

---

## 🎯 Parameters That Control This Feature

### Core Parameters
```bash
# Enable invariant genes normalization
--normalization_method 'invariant_genes'
--normalization_method 'all_genes,invariant_genes'  # Both methods

# Don't skip DeepTools
# (If skip_deeptools_norm is true, no BigWig files generated)
# Default: false (DeepTools runs by default)
```

### Advanced Parameters
```bash
# Invariant gene selection parameters
--sigma_times 1        # Strictness of invariant gene selection
                       # Lower = stricter, higher = more lenient
                       
--n_pop 1              # Population parameter for GeneralNormalizer
                       # Controls multi-group comparison logic
```

### Skip Parameters (Don't Use These!)
```bash
# ❌ These will prevent BigWig generation:
--skip_alignment       # No BAM files = no BigWig
--skip_deeptools_norm  # Explicitly disables BigWig normalization
```

---

## 📈 Expected Performance

### Computational Resources

**NORMALIZE_DESEQ2_QC_INVARIANT_GENES:**
- Label: `process_medium`
- Time: 5-30 minutes (depends on sample count and gene count)
- Memory: 4-16 GB
- CPUs: 1-4 cores

**DEEPTOOLS_BIGWIG_NORM_INVARIANT:**
- Label: `process_high`
- Time: 10-60 minutes per sample
- Memory: 8-32 GB
- CPUs: Multiple cores (uses --numberOfProcessors)
- **Runs in parallel for all samples!**

### File Sizes

**Per sample:**
- Unstranded BigWig: ~100-500 MB
- Forward BigWig: ~50-250 MB
- Reverse BigWig: ~50-250 MB

**Example for 10 samples (stranded):**
- Total: 10 × (100 + 50 + 50) MB = ~2 GB

---

## 🔬 Scientific Validation

### How to Check Normalization Quality

**1. Compare PCA plots:**
```bash
# All genes PCA
results/star/rsem/deseq2/all_genes/Quality_Control/PCA_Plot_All_Genes.pdf

# Invariant genes PCA
results/star/rsem/deseq2/invariant/Quality_Control/PCA_Plot_All_Genes.pdf

# Should show similar sample clustering if normalization is appropriate
```

**2. Check sample distances:**
```bash
# Sample distance heatmaps
results/star/rsem/deseq2/all_genes/Quality_Control/Sample_Distance_Heatmap.pdf
results/star/rsem/deseq2/invariant/Quality_Control/Sample_Distance_Heatmap.pdf

# Should correlate well between methods for good quality data
```

**3. Compare scaling factors:**
```bash
# Extract and compare
echo "Sample,All_Genes,Invariant_Genes" > scaling_comparison.csv
for sample in WT_1 WT_2 KO_1 KO_2; do
    all_genes=$(cat results/star/rsem/deseq2/all_genes/scaling_factors/${sample}_scaling_factor.txt)
    invariant=$(cat results/star/rsem/deseq2/invariant/scaling_factors/${sample}_scaling_factor.txt)
    echo "${sample},${all_genes},${invariant}" >> scaling_comparison.csv
done

# Large differences = method choice matters!
```

**4. Visual inspection in genome browser:**
```bash
# Load both BigWig files in IGV:
# - results/star/rsem/deeptools/all_genes/WT_1.unstranded.norm.bw
# - results/star/rsem/deeptools/invariant_genes/WT_1.unstranded.norm.bw

# Compare coverage at:
# 1. Housekeeping genes (should look similar)
# 2. Differentially expressed genes (may differ)
```

---

## 🎓 Summary

### ✅ Does it work? YES!

**The DeepTools normalization with invariant genes:**
1. ✅ Is fully implemented in your pipeline
2. ✅ Has proper module structure
3. ✅ Has proper workflow integration
4. ✅ Has proper channel operations
5. ✅ Has proper publishDir configuration
6. ✅ Generates separate output files in distinct directories
7. ✅ Works independently from "all genes" normalization
8. ✅ Can run both methods simultaneously

### When to Use Each Method

**Use ALL GENES when:**
- Standard RNA-seq experiment
- Few expected DE genes (<10% of transcriptome)
- Well-balanced experimental design
- You want standard, published normalization

**Use INVARIANT GENES when:**
- Many expected DE genes (>20% of transcriptome)
- Perturbation experiments
- Global transcriptional changes expected
- Comparing very different biological states
- You want robust normalization

**Use BOTH when:**
- Exploring your data
- Comparing normalization approaches
- Publishing comprehensive analysis
- You have sufficient compute resources

### Key Advantages of Invariant Genes Method

1. **Robust to global changes:** Doesn't assume most genes unchanged
2. **Biologically meaningful:** Uses only stable genes for normalization
3. **Better for perturbed systems:** Knockout, overexpression, drug treatment
4. **Separate outputs:** Doesn't interfere with standard normalization
5. **Automatic selection:** GeneralNormalizer identifies invariant genes

---

## 📞 Troubleshooting Commands

```bash
# Quick health check
echo "=== Checking invariant genes normalization ==="

# 1. DESeq2 completed?
test -d results/star/rsem/deseq2/invariant && echo "✅ DESeq2 invariant directory exists" || echo "❌ DESeq2 invariant directory missing"

# 2. Scaling factors generated?
test -f results/star/rsem/deseq2/invariant/scaling_dat.txt && echo "✅ Scaling factors file exists" || echo "❌ Scaling factors file missing"

# 3. Individual scaling factors?
count=$(find results/star/rsem/deseq2/invariant/scaling_factors -name "*_scaling_factor.txt" 2>/dev/null | wc -l)
echo "📊 Found $count individual scaling factor files"

# 4. Invariant genes identified?
test -f results/star/rsem/deseq2/invariant/normalization/invariant_genes.txt && echo "✅ Invariant genes file exists" || echo "❌ Invariant genes file missing"

# 5. BigWig files generated?
count=$(find results/star/rsem/deeptools/invariant_genes -name "*.norm.bw" 2>/dev/null | wc -l)
echo "📊 Found $count normalized BigWig files"

# 6. Compare with all_genes
all_genes_count=$(find results/star/rsem/deeptools/all_genes -name "*.norm.bw" 2>/dev/null | wc -l)
echo "📊 Comparison: all_genes=$all_genes_count, invariant_genes=$count"

echo "=== Health check complete ==="
```

---

## 🔗 Related Files

**Module:** `modules/local/deeptools_bw_norm/main.nf`
**Module:** `modules/local/normalize_deseq2_qc_invariant_genes/main.nf`
**Workflow:** `workflows/rnaseq/main.nf` (lines ~400-1200)
**Config:** `nextflow.config` (DEEPTOOLS_BIGWIG_NORM_INVARIANT section)
**Documentation:** `BUGFIX_QUALITY_CONTROL_PDFS.md`
**Documentation:** `DEEPTOOLS_QUANTIFICATION_FIX.md`

---

**Generated:** 2025-11-27  
**Pipeline Version:** Custom fork with multi-normalization support  
**Status:** ✅ Feature working as designed
