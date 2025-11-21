# Improvement: Organize Per-Sample Genome Counts in Sample Folders

## Overview
Reorganized per-sample genome count files into individual sample folders to match RSEM's organizational structure and improve output clarity.

## Problem Description

### Previous Structure (Less Organized)
Per-sample count files from genome quantification were published directly at the root of the `genome/` directory:

```
results/
└── star/
    └── genome/
        ├── sample1_transcript_counts.txt         ❌ All files at root
        ├── sample1_intron_counts.txt
        ├── sample1_exon_counts.txt
        ├── sample1_5utr_counts.txt
        ├── sample1_3utr_counts.txt
        ├── sample1_combined_counts.txt
        ├── sample1_summary.txt
        ├── sample2_transcript_counts.txt
        ├── sample2_intron_counts.txt
        ├── sample2_exon_counts.txt
        ├── ... (many files at root level)
        ├── genome_exon_counts_merged.txt
        └── genome_transcript_counts_merged.txt
```

**Issues:**
- **Cluttered root directory** with many per-sample files
- **Inconsistent with RSEM** which organizes per-sample results in folders
- **Difficult to navigate** when working with many samples
- **Harder to identify** merged results vs per-sample results

### RSEM Structure (Well Organized)
RSEM already uses sample folders:

```
results/
└── star/
    └── rsem/
        ├── sample1/                              ✅ Per-sample folder
        │   ├── sample1.genes.results
        │   ├── sample1.isoforms.results
        │   ├── sample1.stat/
        │   └── ...
        ├── sample2/                              ✅ Per-sample folder
        │   ├── sample2.genes.results
        │   ├── sample2.isoforms.results
        │   └── ...
        ├── merged_gene_counts.txt                ✅ Merged at root
        └── merged_transcript_counts.txt
```

## Solution

### Changed GENOME_COUNT publishDir Configuration

**Before:**
```groovy
withName: '.*:GENOME_COUNT' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome" },  // ❌ All at root
        mode: params.publish_dir_mode,
        saveAs: { filename -> 
            if (filename.equals('versions.yml') || filename.endsWith('_counts.rds')) {
                return null
            }
            return filename
        }
    ]
}
```

**After:**
```groovy
withName: '.*:GENOME_COUNT' {
    publishDir = [
        path: { "${params.outdir}/${params.aligner}/genome/${meta.id}" },  // ✅ Per-sample folders
        mode: params.publish_dir_mode,
        saveAs: { filename -> 
            if (filename.equals('versions.yml') || filename.endsWith('_counts.rds')) {
                return null
            }
            return filename
        }
    ]
}
```

**Key Change:** Added `/${meta.id}` to the path, creating individual sample folders.

## Expected Output Structure

### New Genome Quantification Structure (Clean & Organized)

```
results/
└── star/
    └── genome/
        ├── ──────────────────────────────────────
        │  PER-SAMPLE FOLDERS (GENOME_COUNT outputs)
        ├── ──────────────────────────────────────
        ├── sample1/                              ✅ Individual folder
        │   ├── sample1_transcript_counts.txt
        │   ├── sample1_intron_counts.txt
        │   ├── sample1_exon_counts.txt
        │   ├── sample1_5utr_counts.txt
        │   ├── sample1_3utr_counts.txt
        │   ├── sample1_combined_counts.txt
        │   └── sample1_summary.txt
        ├── sample2/                              ✅ Individual folder
        │   ├── sample2_transcript_counts.txt
        │   ├── sample2_intron_counts.txt
        │   ├── sample2_exon_counts.txt
        │   ├── sample2_5utr_counts.txt
        │   ├── sample2_3utr_counts.txt
        │   ├── sample2_combined_counts.txt
        │   └── sample2_summary.txt
        ├── sample3/                              ✅ Individual folder
        │   └── ... (same structure)
        │
        ├── ──────────────────────────────────────
        │  MERGED COUNTS AT ROOT (MERGE_GENOME_COUNTS)
        ├── ──────────────────────────────────────
        ├── genome_transcript_counts_merged.txt   ✅ Easy to find
        ├── genome_intron_counts_merged.txt
        ├── genome_exon_counts_merged.txt
        ├── genome_5utr_counts_merged.txt
        ├── genome_3utr_counts_merged.txt
        ├── genome_transcript_merge_summary.txt
        ├── genome_intron_merge_summary.txt
        ├── genome_exon_merge_summary.txt
        ├── genome_5utr_merge_summary.txt
        ├── genome_3utr_merge_summary.txt
        │
        ├── ──────────────────────────────────────
        │  DESEQ2 NORMALIZATION & QC
        ├── ──────────────────────────────────────
        ├── deseq2/
        │   ├── all_genes/
        │   │   ├── Quality_Control/
        │   │   └── Read_Distribution/
        │   └── invariant_genes/
        │       ├── Quality_Control/
        │       └── Read_Distribution/
        │
        ├── ──────────────────────────────────────
        │  DEEPTOOLS NORMALIZED BIGWIG
        ├── ──────────────────────────────────────
        └── deeptools/
            ├── all_genes/
            │   ├── sample1.norm.bw
            │   ├── sample2.norm.bw
            │   └── sample3.norm.bw
            └── invariant_genes/
                ├── sample1.norm.bw
                ├── sample2.norm.bw
                └── sample3.norm.bw
```

### HISAT2 Structure (Same Organization)
```
results/
└── hisat2/
    └── genome/
        ├── sample1/                              ✅ Individual folders
        ├── sample2/
        ├── sample3/
        ├── genome_*_counts_merged.txt           ✅ Merged at root
        ├── deseq2/
        └── deeptools/
```

## Comparison: Before vs After

### Before (Cluttered)
```
star/genome/
├── sample1_exon_counts.txt           ← 7 files per sample
├── sample1_transcript_counts.txt
├── sample1_intron_counts.txt
├── sample1_5utr_counts.txt
├── sample1_3utr_counts.txt
├── sample1_combined_counts.txt
├── sample1_summary.txt
├── sample2_exon_counts.txt           ← Another 7 files
├── sample2_transcript_counts.txt
├── ... (35 files for 5 samples!)
├── genome_exon_counts_merged.txt     ← Hard to find
└── genome_transcript_counts_merged.txt
```

**With 10 samples:** 70+ files at root level! 😱

### After (Organized)
```
star/genome/
├── sample1/                          ← Organized in folder
│   └── 7 count files
├── sample2/                          ← Organized in folder
│   └── 7 count files
├── sample3/
├── ... (only 10 folders)
├── genome_exon_counts_merged.txt     ← Easy to find!
└── genome_transcript_counts_merged.txt
```

**With 10 samples:** 10 folders + merged files at root! 🎉

## Consistency with RSEM

### RSEM Structure
```
star/rsem/
├── sample1/          ← Per-sample folder
├── sample2/          ← Per-sample folder
├── merged_*.txt      ← Merged at root
├── deseq2/
└── deeptools/
```

### Genome Structure (Now Matches!)
```
star/genome/
├── sample1/          ← Per-sample folder ✅
├── sample2/          ← Per-sample folder ✅
├── genome_*.txt      ← Merged at root ✅
├── deseq2/           ✅
└── deeptools/        ✅
```

**Perfect consistency!** Both quantification methods follow the same organizational pattern.

## Benefits

### 1. Clean Root Directory ✨
- Only merged files and subdirectories at root level
- Easy to identify key analysis files
- Professional and organized output structure

### 2. Sample-Level Organization 📁
- All per-sample results grouped together
- Easy to inspect individual sample quality
- Easy to share sample-specific results

### 3. Consistency Across Pipeline 🔄
- RSEM uses sample folders ✅
- Genome now uses sample folders ✅
- DESeq2 uses method folders ✅
- DeepTools uses method folders ✅

### 4. Scalability 📊
- Structure remains clean with 100+ samples
- Easy to navigate regardless of sample count
- Scripts can easily iterate over sample folders

### 5. Clarity for Users 👥
- Clear separation: per-sample vs merged results
- Intuitive organization matches user expectations
- Reduces confusion about which files to use

## Files Modified

### Modified
- **Workflow Configuration:** `workflows/rnaseq/nextflow.config`
  - Changed `GENOME_COUNT` publishDir path from `${params.aligner}/genome` to `${params.aligner}/genome/${meta.id}`
  - Single line change with significant organizational improvement

## Impact Assessment

### No Breaking Changes ✅
- Configuration change only affects output organization
- No changes to file content or naming
- No changes to workflow logic
- Downstream processes (DESeq2, DeepTools) unaffected (they use merged counts)

### Improved User Experience ✅
- Cleaner, more professional output structure
- Easier to navigate and understand results
- Consistent with established bioinformatics conventions

### Minor Migration Note
Users upgrading from previous versions will see a new folder structure:
- **Old:** `star/genome/sample1_exon_counts.txt`
- **New:** `star/genome/sample1/sample1_exon_counts.txt`

**Migration is automatic** - just re-run the pipeline and results will be organized in the new structure.

## Testing Recommendations

### Test 1: Verify Sample Folders Created
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner star \
  --quantification genome

# Verify sample folders exist
ls -lh results/star/genome/
# Should show: sample1/ sample2/ sample3/ genome_*.txt

# Verify per-sample files are in folders
ls -lh results/star/genome/sample1/
# Should show: sample1_exon_counts.txt, sample1_transcript_counts.txt, etc.
```

### Test 2: Verify Merged Files at Root
```bash
# Merged files should be at root level, not in sample folders
ls results/star/genome/*.txt
# Should show: genome_exon_counts_merged.txt, genome_transcript_counts_merged.txt, etc.
```

### Test 3: Verify DESeq2 Still Works
```bash
# DESeq2 uses merged exon counts from root
ls results/star/genome/deseq2/all_genes/Quality_Control/*.pdf
# Should show: heatmaps, PCA plots, etc.
```

### Test 4: Verify HISAT2 Structure
```bash
nextflow run . \
  --input samplesheet.csv \
  --aligner hisat2 \
  --quantification genome

# Verify same organization for HISAT2
ls -lh results/hisat2/genome/
# Should show sample folders + merged files
```

## Technical Details

### How meta.id Works
The `meta.id` variable contains the sample identifier from the samplesheet. The publishDir path interpolation `${meta.id}` creates a folder with the sample name.

**Example:**
- Sample ID: `WT_REP1`
- Output path: `results/star/genome/WT_REP1/`
- Files: `WT_REP1_exon_counts.txt`, etc.

### Workflow Compatibility
The `MERGE_GENOME_COUNTS` process collects all per-sample files regardless of folder structure:
```groovy
MERGE_GENOME_COUNTS (
    GENOME_COUNT.out.combined_counts.map { meta, counts -> counts }.collect(),
    ...
)
```

The `.collect()` operator gathers files from all sample folders automatically, so the merge process is unaffected by the folder organization.

## Related Changes

This improvement builds on:
1. **Unified genome quantification structure** - Both STAR and HISAT2 support
2. **Quality Control PDFs** - Published in organized subdirectories
3. **DeepTools BigWig** - Published in organized subdirectories

All improvements work together to create a clean, professional, and consistent output structure.

## Summary

✅ **Single-line configuration change**
✅ **Significant organizational improvement**
✅ **Consistent with RSEM structure**
✅ **Scales well with sample count**
✅ **No breaking changes to analysis**
✅ **Improved user experience**

The pipeline now follows best practices for output organization:
- Per-sample results in sample folders
- Merged/summarized results at root level
- Analysis results in organized subdirectories
