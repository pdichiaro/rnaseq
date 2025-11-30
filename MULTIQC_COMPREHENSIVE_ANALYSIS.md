# Comprehensive MultiQC Workflow Analysis

## Executive Summary

After thorough analysis of the entire RNA-seq workflow, I found **2 critical bugs** and **1 potential issue** with MultiQC:

1. ✅ **FIXED**: Invalid aligner conditional check
2. ✅ **FIXED**: Missing publishDir directive  
3. ⚠️ **POTENTIAL ISSUE**: No fallback MultiQC if no aligner/pseudo-aligner is used

---

## Complete MultiQC Flow

### 1. MultiQC Modules in Use

The pipeline uses **3 different MultiQC process definitions**:

| Module | Location | Usage | Status |
|--------|----------|-------|--------|
| `MULTIQC` | `modules/nf-core/multiqc/main.nf` | Included but **NOT USED** | ⚠️ Unused |
| `MULTIQC_WITH_SUBFOLDERS` | `modules/local/multiqc_with_subfolders/main.nf` | Used for all aligner-specific reports | ✅ Fixed |
| `MULTIQC_CUSTOM_BIOTYPE` | `modules/local/multiqc_custom_biotype/main.nf` | Generates biotype count tables | ✅ Working |

### 2. MultiQC Instances Created

The workflow creates **separate MultiQC reports** for each aligner/pseudo-aligner:

```groovy
// Line 77-79: Three instances of MULTIQC_WITH_SUBFOLDERS
MULTIQC_STAR     → when params.aligner == 'star'
MULTIQC_HISAT2   → when params.aligner == 'hisat2'
MULTIQC_KALLISTO → when params.pseudo_aligner == 'kallisto'
```

---

## Channel Flow Analysis

### Input Channels to MultiQC

MultiQC collects files from multiple channels throughout the workflow:

#### **Shared QC Files** (ch_multiqc_shared)
Goes into ALL MultiQC reports:
- FastQC (raw and trimmed) - Line 291-292
- Trimming stats (Fastp/Cutadapt) - Line 291-292
- Strandedness inference - Line 291-292
- Sample status checks - Line 771
- Workflow summary - Line 1284-1288
- Software versions - Line 1284-1288
- Methods description - Line 1284-1288

#### **STAR-Specific Files** (ch_multiqc_star_files + ch_multiqc_star_qc + ch_multiqc_star_quant)
- STAR alignment logs - Line 339-340
- BAM statistics (if not using UMI) - Line 368-370
- UMI deduplication stats (if using UMI) - Line 363-364
- RSEM quantification stats - Line 390-391
- DESeq2 QC plots - Line 564-576

#### **HISAT2-Specific Files** (ch_multiqc_hisat2_files + ch_multiqc_hisat2_qc + ch_multiqc_hisat2_quant)
- HISAT2 alignment summary - Line 602-603
- BAM statistics - Line 630-633
- DESeq2 QC plots - Line 716-728

#### **Kallisto-Specific Files** (ch_multiqc_kallisto_files)
- Kallisto quantification logs - Line 1055
- DESeq2 QC plots - Line 1128-1138

#### **Common Post-Alignment QC** (ch_multiqc_files)
- Preseq library complexity - Line 793
- Picard MarkDuplicates metrics - Line 807-810
- Biotype counts - Line 850
- Qualimap RNA-seq metrics - Line 910
- dupRadar duplication analysis - Line 919
- RSeQC comprehensive QC - Line 937-944
- Strandedness comparison - Line 1001
- Kraken2/Bracken contamination - Line 1015, 1022

---

## Critical Issues Found and Fixed

### Bug #1: Invalid Aligner Conditional ✅ FIXED

**Location:** `workflows/rnaseq/main.nf` line 1317

**Problem:**
```groovy
if (params.aligner && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem' || params.aligner == 'star'))
```

**Why Wrong:**
- Valid aligners: `['star', 'hisat2']` (line 144)
- `'star_salmon'` and `'star_rsem'` are **not valid** aligner values
- These are legacy from old API where aligner+quantifier were combined
- Modern API: separate `params.aligner` and `params.quantification`

**Fix Applied:**
```groovy
if (params.aligner == 'star')
```

**Commit:** `0f9b338`

---

### Bug #2: Missing publishDir Directive ✅ FIXED

**Location:** `modules/local/multiqc_with_subfolders/main.nf`

**Problem:**
The MULTIQC_WITH_SUBFOLDERS process had **NO publishDir directive**, so outputs remained in work directory only.

**Impact:**
- Process executed successfully
- HTML report generated
- But files were **never published** to results directory
- Users couldn't find their MultiQC reports

**Fix Applied:**
```groovy
publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode, 
    saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
```

**Commit:** `5dcceec`

---

## Potential Issue #3: No Fallback MultiQC ⚠️

### The Problem

**Current Logic:**
```groovy
if (!params.skip_multiqc) {
    // ... setup channels ...
    
    if (params.aligner == 'star') {
        MULTIQC_STAR(...)
    }
    if (params.aligner == 'hisat2') {
        MULTIQC_HISAT2(...)
    }
    if (params.pseudo_aligner == 'kallisto') {
        MULTIQC_KALLISTO(...)
    }
}
```

**Issue:** If someone runs with:
- `--skip_alignment` (no aligner)
- `--skip_pseudo_alignment` (no pseudo-aligner)
- Or invalid aligner somehow bypasses validation

Then **NO MultiQC report will be generated** even though QC files exist (FastQC, trimming stats, etc.)!

### Recommendation

Add a fallback generic MultiQC for shared QC when no aligner is used:

```groovy
if (!params.skip_multiqc) {
    // ... existing setup ...
    
    ch_multiqc_report = Channel.empty()
    
    // STAR MultiQC Report
    if (params.aligner == 'star') {
        MULTIQC_STAR(...)
        ch_multiqc_report = ch_multiqc_report.mix(MULTIQC_STAR.out.report)
    }
    // HISAT2 MultiQC Report  
    else if (params.aligner == 'hisat2') {
        MULTIQC_HISAT2(...)
        ch_multiqc_report = ch_multiqc_report.mix(MULTIQC_HISAT2.out.report)
    }
    // Kallisto MultiQC Report
    else if (params.pseudo_aligner == 'kallisto') {
        MULTIQC_KALLISTO(...)
        ch_multiqc_report = ch_multiqc_report.mix(MULTIQC_KALLISTO.out.report)
    }
    // Fallback: Generate MultiQC with shared QC only
    else if (params.skip_alignment || (!params.aligner && !params.pseudo_aligner)) {
        MULTIQC(  // Use the nf-core/multiqc module
            ch_multiqc_shared.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = ch_multiqc_report.mix(MULTIQC.out.report)
    }
}
```

**Note:** This would require adding publishDir to `modules/nf-core/multiqc/main.nf` as well.

---

## Configuration Analysis

### Parameters (nextflow.config)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_multiqc` | `false` | Skip MultiQC report generation |
| `multiqc_config` | `null` | Path to custom MultiQC config |
| `multiqc_title` | `null` | Custom title for report |
| `multiqc_logo` | `null` | Path to custom logo |
| `multiqc_methods_description` | `null` | Custom methods description |

### MultiQC Config (workflows/rnaseq/assets/multiqc/multiqc_config.yml)

**Key Settings:**
- Report sections ordered by priority (5000 = highest)
- Separate sections for each quantifier (deseq2-star-genome-qc, etc.)
- Module search patterns optimized for speed
- Specific modules enabled: fastqc, cutadapt, fastp, star, hisat2, rsem, kallisto, samtools, picard, preseq, rseqc, qualimap, kraken
- FastQC split into raw and trimmed sections
- % Dups hidden (shown by Picard instead)

---

## Output Structure

With fixes applied, MultiQC outputs will be published to:

```
${params.outdir}/
└── multiqc/
    ├── multiqc_report.html       # Main HTML report
    ├── multiqc_data/              # Parsed data in JSON/TSV
    │   ├── multiqc_general_stats.txt
    │   ├── multiqc_sources.txt
    │   ├── multiqc_data.json
    │   └── ...
    └── multiqc_plots/             # Optional plot data
        └── ...
```

**Note:** Report name can be customized with `--multiqc_title` parameter.

---

## Verification Steps

### 1. Check if MultiQC Will Run

```bash
# Check parameters in your run
params.skip_multiqc == false  # Should be false (default)
params.aligner == 'star'      # Or 'hisat2'
# OR
params.pseudo_aligner == 'kallisto'
```

### 2. Monitor Execution

Look for these process names in logs:
```
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:MULTIQC_STAR (1)
# OR
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:MULTIQC_HISAT2 (1)
# OR
[xx/xxxxxx] process > NFCORE_RNASEQ:RNASEQ:MULTIQC_KALLISTO (1)
```

### 3. Verify Output

```bash
# After completion
ls -lh ${OUTDIR}/multiqc/
# Should show:
# - multiqc_report.html
# - multiqc_data/
```

### 4. Check MultiQC Work Directory

If report isn't published, check work directory:
```bash
find work -name "multiqc_report.html" -type f
```

If found in work but not in results → publishDir issue
If not found anywhere → MultiQC didn't run

---

## Summary

### ✅ Fixed Issues
1. Conditional check now correctly uses `params.aligner == 'star'`
2. publishDir directive added to MULTIQC_WITH_SUBFOLDERS module
3. Both fixes committed and pushed to main branch

### ⚠️ Potential Enhancement
- Consider adding fallback MultiQC for workflows without alignment
- Would require enabling the unused `MULTIQC` module with publishDir

### 📊 Current Status
- MultiQC will execute for: STAR, HISAT2, Kallisto
- MultiQC reports will publish to: `${params.outdir}/multiqc/`
- All QC data properly collected and integrated

---

## Testing Recommendations

1. **Run with -resume:**
   ```bash
   nextflow run . -resume --outdir results
   ```
   Cached processes will be reused, only MultiQC will re-run with fixes.

2. **Check logs:**
   ```bash
   grep -i "MULTIQC" .nextflow.log
   ```

3. **Verify output:**
   ```bash
   ls -lh results/multiqc/
   cat results/multiqc/multiqc_report.html  # Should exist
   ```

4. **Test different aligners:**
   ```bash
   # Test STAR (default)
   nextflow run . --aligner star
   
   # Test HISAT2
   nextflow run . --aligner hisat2
   
   # Test Kallisto
   nextflow run . --pseudo_aligner kallisto --skip_alignment
   ```

---

## Date
2025-11-30

## Files Modified
- `workflows/rnaseq/main.nf` (line 1317) - Fixed conditional
- `modules/local/multiqc_with_subfolders/main.nf` (line 9) - Added publishDir
- `BUGFIX_MULTIQC.md` - Documentation
- `MULTIQC_COMPREHENSIVE_ANALYSIS.md` - This file
