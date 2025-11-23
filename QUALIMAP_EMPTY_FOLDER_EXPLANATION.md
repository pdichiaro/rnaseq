# Qualimap Empty Folder Issue - Explanation

## What is Qualimap?

**Qualimap** is a quality control tool for RNA-seq data that analyzes BAM alignment files. It provides:

1. **Mapping Quality Metrics**:
   - Coverage statistics across genes
   - Read distribution across genomic features (exons, introns, intergenic)
   - Junction analysis
   - 5'-3' bias assessment

2. **Output Files**:
   - `qualimapReport.html` - Interactive HTML report
   - `rnaseq_qc_results.txt` - Text summary of QC metrics
   - Coverage plots and statistics
   - Genomic origin plots (showing % reads in exons, introns, etc.)

3. **Purpose**: 
   - Detect biases in sequencing/mapping
   - Assess overall data quality
   - Verify appropriate coverage of gene features
   - Identify problems before downstream analysis

## Why is the Qualimap Folder Empty?

The Qualimap folder is empty because of **sample filtering based on mapping percentage**. Here's the workflow logic:

### Pipeline Flow for Qualimap:

```
1. STAR Alignment generates BAM files
   ↓
2. Calculate percent mapped reads for each sample
   ↓
3. FILTER: Keep only samples with ≥ min_mapped_reads % mapped
   ↓
4. Qualimap runs on FILTERED BAM files (ch_genome_bam)
   ↓
5. If NO samples pass filter → ch_genome_bam is EMPTY → Qualimap produces nothing
```

### The Filtering Logic (workflows/rnaseq/main.nf, lines 655-699):

```groovy
ch_genome_bam_bai_mapping = ch_genome_bam
    .join(ch_genome_bam_index)
    .join(ch_percent_mapped, remainder: true)
    .map{ row ->
        def (meta, bam, index) = row[0..2]
        def percent_mapped = row.size() == 4 ? row[3] : null
        // ⚠️ KEY LINE: Filter samples below threshold
        def pass = percent_mapped != null ? percent_mapped >= params.min_mapped_reads.toFloat() : null
        return [ meta, bam, index, percent_mapped, pass ]
    }

// ... later ...

// ⚠️ CRITICAL: Only samples that PASS the filter are kept
map_filtered_genome_bam_bai = ch_genome_bam_bai_mapping.bam
    .filter { meta, bam, index, pass -> pass || pass == null }

ch_genome_bam = map_filtered_genome_bam_bai.bam

// ... later ...

// ⚠️ Qualimap uses the FILTERED channel
if (!params.skip_qualimap) {
    QUALIMAP_RNASEQ (
        ch_genome_bam,  // ← If this is empty, no Qualimap output
        ch_gtf.map { [ [:], it ] }
    )
}
```

## Common Reasons for Empty Qualimap Folder

### 1. **Low Mapping Percentage** (Most Common)
   - **Default threshold**: `min_mapped_reads = 5` (5% mapped reads)
   - If all samples have < 5% uniquely mapped reads, they're filtered out
   - **Check**: Look at STAR log files in `star/log/` for mapping percentages

### 2. **Skip Flag Enabled**
   - `--skip_qualimap true` in your run command
   - `--skip_qc true` disables all QC including Qualimap

### 3. **No Genome BAM Files Generated**
   - Using kallisto pseudoaligner without STAR
   - Qualimap requires genome-aligned BAM files

### 4. **Alignment Issues**
   - STAR alignment failed for all samples
   - Genome index mismatch
   - Reference genome/GTF issues

## How to Diagnose

### Step 1: Check Your Run Parameters
```bash
# Check if Qualimap is disabled
cat .nextflow.log | grep "skip_qualimap"

# Check min_mapped_reads threshold
cat .nextflow.log | grep "min_mapped_reads"
```

### Step 2: Check STAR Mapping Statistics
```bash
# Look at STAR final logs
cat results/star/log/*Log.final.out | grep "Uniquely mapped reads %"

# Example output:
# Uniquely mapped reads %: 2.34%  ← TOO LOW! Below 5% threshold
# Uniquely mapped reads %: 85.67% ← GOOD! Above threshold
```

### Step 3: Check Failed Samples Report
```bash
# MultiQC should report failed samples
cat results/multiqc/multiqc_data/fail_mapped_samples_mqc.tsv

# Or check the work directory
find work -name "fail_mapped_samples_mqc.tsv" -exec cat {} \;
```

### Step 4: Verify ch_genome_bam Contents
```bash
# Check if BAM files were generated
ls -lh results/star/

# Check if BAM files exist but were filtered
ls -lh work/*/*/*.bam | wc -l
```

## Solutions

### Solution 1: Lower the Mapping Threshold (if using test data)
```bash
nextflow run workflows/rnaseq/main.nf \
  --min_mapped_reads 1 \    # ← Lower from 5 to 1
  ...other params...
```

**Use when**: Testing with small/subsampled data that naturally has low mapping rates

### Solution 2: Investigate Low Mapping Rates
If all real samples have < 5% mapping:
- **Check reference genome**: Is it the correct species?
- **Check GTF annotation**: Does it match the genome version?
- **Check data quality**: Run FastQC - are reads very short/low quality?
- **Check for contamination**: Reads might be from a different organism
- **Check library prep**: Was it stranded? Did you specify strandedness correctly?

### Solution 3: Run Qualimap Manually (After Pipeline)
If you want to run Qualimap on filtered samples anyway:

```bash
# Find your BAM files
BAM_FILE="results/star/sample1.Aligned.sortedByCoord.out.bam"
GTF_FILE="path/to/genes.gtf"

# Run Qualimap manually
qualimap rnaseq \
  -bam $BAM_FILE \
  -gtf $GTF_FILE \
  -outdir qualimap_manual_output \
  -p strand-specific-reverse \  # Adjust based on your library
  -pe  # Include for paired-end data
```

### Solution 4: Disable Sample Filtering (Advanced)
Modify `workflows/rnaseq/main.nf` line 699 to NOT filter:

```groovy
# BEFORE (filters samples):
map_filtered_genome_bam_bai = ch_genome_bam_bai_mapping.bam
    .filter { meta, bam, index, pass -> pass || pass == null }

# AFTER (keeps all samples):
map_filtered_genome_bam_bai = ch_genome_bam_bai_mapping.bam
    .map { meta, bam, index, pass -> 
        bam: [meta, bam]
        index: [meta, index]
    }
```

⚠️ **Warning**: This will include low-quality samples in downstream analysis!

## Expected Qualimap Output Structure

When working correctly:
```
results/star/qualimap/
├── sample1/
│   ├── qualimapReport.html          # Main HTML report
│   ├── rnaseq_qc_results.txt        # Text summary
│   ├── images_qualimapReport/       # Plot images
│   ├── raw_data_qualimapReport/     # Raw data files
│   └── css/                         # Styling
├── sample2/
│   └── ...
└── sample3/
    └── ...
```

## Key MultiQC Plots from Qualimap

When Qualimap runs successfully, MultiQC includes:

1. **Gene Coverage Profile**: Shows 5' to 3' coverage across gene bodies
   - Detects 5' or 3' bias
   - Indicates RNA degradation issues

2. **Genomic Origin of Reads**: Bar chart showing read distribution
   - % in exons (should be highest for RNA-seq)
   - % in introns
   - % in intergenic regions
   - % in ribosomal RNA

## Configuration Parameters

Relevant parameters in `nextflow.config`:
```groovy
// Line 98: Enable/disable Qualimap
skip_qualimap = false

// Line 103: Minimum mapping threshold
min_mapped_reads = 5  // Samples below 5% mapped are filtered

// Line 94: Disable all QC (including Qualimap)
skip_qc = false
```

## Summary

**Qualimap folder is empty when**:
1. ✅ All samples fail the `min_mapped_reads` threshold (< 5% mapped)
2. ✅ Qualimap is disabled (`--skip_qualimap true`)
3. ✅ No genome BAM files are generated

**Most likely cause**: Your samples have low mapping percentages and are filtered out before Qualimap runs.

**Quick fix for testing**: Use `--min_mapped_reads 1` or check STAR logs to verify mapping rates.
