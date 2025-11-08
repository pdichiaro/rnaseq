# Repository Structure Overview

## Current Date/Time
Generated: 2025-11-08 at 21:35:46

## Repository: pdichiaro/rnaseq
**Branch:** main  
**Last Commit:** 85da866 - Fix DESeq2 MultiQC integration for pseudo-alignment quantifiers

---

## Root Directory Files

### Core Nextflow Files
- `main.nf` - Main entry point for the pipeline
- `nextflow.config` - Global pipeline configuration
- `nextflow_schema.json` - Parameter validation schema
- `test.config` - Test-specific configuration

### Documentation Files
- **`DESEQ2_MULTIQC_FIX.md`** ✅ - Comprehensive fix documentation (270 lines)
- `DESEQ2_STRUCTURE.txt` - Analysis of DESeq2 workflow structure
- `config_params.txt` - Parameter documentation
- `ro-crate-metadata.json` - Research Object metadata
- `tower.yml` - Seqera Platform configuration
- `modules.json` - nf-core modules tracking

---

## Directory Structure

```
rnaseq/
├── assets/                    # Static pipeline assets
├── bin/                       # Executable scripts used by processes
├── conf/                      # Additional configuration files
├── docs/                      # Pipeline documentation
├── modules/                   # nf-core modules
│   └── local/                 # Custom local modules
│   └── nf-core/              # Community modules
├── subworkflows/             # Reusable workflow components
│   └── local/
│   └── nf-core/
├── test_multiqc/             # MultiQC testing resources
└── workflows/
    └── rnaseq/               # Main RNA-seq workflow
        ├── main.nf           # ✅ MODIFIED - Contains the fix
        ├── nextflow.config
        └── assets/
```

---

## Key Modified Files

### 1. `workflows/rnaseq/main.nf` (Lines 1000-1040)

**Modified:** 2025-11-08 21:09:00  
**Purpose:** Fix DESeq2 MultiQC integration

**Change Details:**
```groovy
// Line ~1004: NEW - Explicit quantifier determination
def pseudo_quantifier = params.pseudo_aligner ? params.pseudo_aligner.capitalize() : "Salmon"

// Lines ~1008, 1023: UPDATED - Use pseudo_quantifier instead of params.pseudo_aligner
NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO (
    ch_counts_gene,
    pseudo_quantifier  // ← Changed from: params.pseudo_aligner ?: "Salmon"
)

NORMALIZE_DESEQ2_QC_ALL_GENES_PSEUDO (
    ch_counts_gene,
    pseudo_quantifier  // ← Changed from: params.pseudo_aligner ?: "Salmon"
)
```

**Impact:**
- Ensures quantifier is properly capitalized ("Kallisto" or "Salmon")
- Provides consistent value to both normalization processes
- Fixes file naming mismatch that prevented MultiQC detection

---

## Git History (Recent)

```
85da866 (HEAD -> main, origin/main) Fix DESeq2 MultiQC integration for pseudo-alignment quantifiers
6a0b7ea Fix DESeq2 QC MultiQC integration and remove non-existent output channels
c3cd377 Fix MultiQC DESeq2 section configuration
b429ae3 Fix duplicate custom_data declaration in MultiQC config
a14569b Remove unused dendrogram generation scripts
```

**Commit 85da866 Changes:**
- ✅ Added: `DESEQ2_MULTIQC_FIX.md` (270 lines)
- ✅ Modified: `workflows/rnaseq/main.nf` (8 lines changed)

---

## Clean Status

✅ **No temporary files remaining**
- Removed: CRITICAL_VERIFICATION.md
- Removed: HOW_TO_APPLY_FIX.md
- Removed: DIAGNOSTIC_COMMANDS.sh
- Removed: QUICK_CHECK.sh
- Removed: run_with_fix.sh
- Removed: RUN_INSTRUCTIONS.txt

✅ **Git working directory clean**
- All changes committed and pushed
- No untracked files
- No staged changes

---

## Documentation Summary

### DESEQ2_MULTIQC_FIX.md Contents

1. **Problem Summary** - Root cause analysis
2. **Solution Applied** - Code changes with before/after comparison
3. **Technical Details** - How the fix works
4. **File Naming Convention** - Expected output patterns
5. **Verification Steps** - How to confirm the fix works
6. **Testing Recommendations** - Best practices for validation
7. **Additional Notes** - Edge cases and considerations

**Total:** 270 lines of comprehensive documentation

---

## Expected Output Files (After Fix)

When running with Kallisto (`--pseudo_aligner 'kallisto'`):

```
results/deseq2_qc/
├── kallisto.deseq2.all_genes.pca.vals.txt
├── kallisto.deseq2.all_genes.sample.dists.txt
├── kallisto.deseq2.all_genes.read.distribution.normalized.txt
├── kallisto.deseq2.all_genes.size_factors.txt
├── kallisto.deseq2.all_genes.pca.top500.vals.txt
└── kallisto.deseq2.all_genes.rds
```

When running with Salmon (`--pseudo_aligner 'salmon'` or default):

```
results/deseq2_qc/
├── salmon.deseq2.all_genes.pca.vals.txt
├── salmon.deseq2.all_genes.sample.dists.txt
├── salmon.deseq2.all_genes.read.distribution.normalized.txt
├── salmon.deseq2.all_genes.size_factors.txt
├── salmon.deseq2.all_genes.pca.top500.vals.txt
└── salmon.deseq2.all_genes.rds
```

---

## MultiQC Detection Patterns

From MultiQC configuration, these are the expected file patterns:

```yaml
custom_data:
  deseq2_kallisto_qc:
    file_format: 'tsv'
    section_name: 'DESeq2 Kallisto QC'
    pconfig:
      - 'kallisto.deseq2.*.pca.vals.txt'
      - 'kallisto.deseq2.*.sample.dists.txt'
      - 'kallisto.deseq2.*.read.distribution.normalized.txt'
  
  deseq2_salmon_qc:
    file_format: 'tsv'
    section_name: 'DESeq2 Salmon QC'
    pconfig:
      - 'salmon.deseq2.*.pca.vals.txt'
      - 'salmon.deseq2.*.sample.dists.txt'
      - 'salmon.deseq2.*.read.distribution.normalized.txt'
```

**The fix ensures files are named to match these patterns exactly.**

---

## Next Steps for Users

### 1. Update Pipeline
```bash
nextflow pull pdichiaro/rnaseq -r main
```

### 2. Run Pipeline
```bash
nextflow run pdichiaro/rnaseq \
  -r main \
  -resume \
  --pseudo_aligner 'kallisto' \
  -profile singularity \
  --input samplesheet.csv \
  --outdir results \
  --fasta genome.fa \
  --gtf annotation.gtf \
  --transcript_fasta transcripts.fa \
  --normalization_method 'all_genes'
```

### 3. Verify Fix
```bash
# Check MultiQC report
firefox results/multiqc/multiqc_report.html

# Check DESeq2 output files
ls -la results/deseq2_qc/kallisto.deseq2.*
```

**Expected Result:** DESeq2 section appears in MultiQC report with all plots

---

## Repository URLs

- **GitHub:** https://github.com/pdichiaro/rnaseq
- **Branch:** main
- **Commit:** 85da866bdae3ed2678f7eab8b60b7eb89ad9d6cd

---

**Generated by:** Seqera AI  
**Workspace:** showcase (ID: 40230138858677)  
**Organization:** community
