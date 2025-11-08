# Suggestion: Enhance Kallisto n_unique Reporting in MultiQC

## Current Status

✅ **Already Enabled!** The MultiQC configuration already includes `n_unique` visibility:

```yaml
table_columns_visible:
  "Kallisto":
    n_unique: true
```

**Location:** `workflows/rnaseq/assets/multiqc/multiqc_config.yml` (line 113)

---

## What is `n_unique`?

The **`n_unique`** metric from Kallisto represents:
- **Number of uniquely mapped reads** to transcripts
- **Important QC metric** for assessing quantification confidence
- **Complementary to total processed reads**

### Why It Matters:
1. **Multi-mapping assessment** - Low n_unique suggests many reads map to multiple transcripts
2. **Ambiguity detection** - Helps identify samples with high isoform ambiguity
3. **Quantification confidence** - Higher n_unique generally means more reliable quantification
4. **Library complexity** - Can indicate library preparation issues

---

## Current Implementation

### ✅ What's Already Working

**1. General Statistics Table:**
- `n_unique` is **visible by default** in MultiQC General Statistics
- Shows up alongside other Kallisto metrics (processed reads, pseudoaligned %)
- Available from Kallisto's `run_info.json` output

**2. Data Source:**
```json
// From: results/kallisto/{sample}/run_info.json
{
  "n_targets": 246605,
  "n_bootstraps": 0,
  "n_processed": 50000000,
  "n_pseudoaligned": 43500000,
  "n_unique": 38500000,    // ← This value
  "p_pseudoaligned": 87.0,
  "p_unique": 77.0,
  "kallisto_version": "0.48.0",
  "start_time": "..."
}
```

---

## Enhancement Suggestions

### 🎯 **Option 1: Add Dedicated n_unique Visualization Section**

**Why:** Create a focused visualization specifically for n_unique metrics

**Implementation:**

Add to `multiqc_config.yml` in the `custom_data` section:

```yaml
custom_data:
  kallisto-n-unique-stats:
    section_name: "Kallisto Uniquely Mapped Reads Analysis"
    description: >
      Analysis of uniquely mapped reads from Kallisto quantification.
      The n_unique metric represents reads that map unambiguously to a single transcript.
      Higher percentages indicate better quantification confidence and lower isoform ambiguity.
    section_anchor: kallisto-unique-analysis
    file_format: tsv
    plot_type: bar_graph
    pconfig:
      title: "Kallisto Unique Read Mapping Rates"
      ylab: "Percentage of Reads (%)"
      xlab: "Sample"
      id: kallisto_n_unique_percentage
      cpswitch_counts_label: "Number of Reads"
      data_labels:
        - name: "Unique Reads (%)"
          ylab: "Percentage"
        - name: "Multi-mapped Reads (%)"
          ylab: "Percentage"
```

**Required:** Create a custom script to generate this data from `run_info.json`

---

### 🎯 **Option 2: Enhanced Kallisto Section with n_unique Plots**

**Why:** Expand the existing Kallisto section with additional metrics

**Benefits:**
- Compare n_unique across samples
- Visualize unique vs multi-mapped reads
- Track n_unique as QC threshold

**Add to Report Section Order:**
```yaml
report_section_order:
  kallisto:
    order: 2000
  kallisto-unique-analysis:
    order: 1999  # Right after main Kallisto section
```

---

### 🎯 **Option 3: Add n_unique to DESeq2 QC Context**

**Why:** Show n_unique alongside DESeq2 normalization metrics

**Rationale:**
- Samples with low n_unique may have different normalization behavior
- Useful for identifying outliers in PCA/sample distance plots
- Can explain DESeq2 QC anomalies

**Implementation:**

Create a summary table showing:
```
Sample | n_unique | n_unique (%) | DESeq2 Size Factor | PCA Outlier
-------|----------|--------------|-------------------|-------------
S1     | 38.5M    | 77.0%        | 1.02              | No
S2     | 25.3M    | 50.6%        | 0.87              | Yes  ← Low n_unique!
S3     | 40.1M    | 80.2%        | 1.05              | No
```

---

## Recommended Implementation Plan

### ✅ **Phase 1: Verify Current Functionality** (DONE)
- [x] Confirm `n_unique: true` in multiqc_config.yml
- [x] Check Kallisto module outputs `run_info.json`
- [x] Verify MultiQC detects n_unique automatically

### 🔧 **Phase 2: Add Custom n_unique Visualization** (RECOMMENDED)

**2.1. Create n_unique Summary Script**

Location: `bin/kallisto_unique_summary.py`

```python
#!/usr/bin/env python3
"""
Generate Kallisto n_unique summary for MultiQC
"""
import json
import sys
from pathlib import Path

def parse_kallisto_stats(results_dir):
    """Parse Kallisto run_info.json files"""
    stats = {}
    
    for json_file in Path(results_dir).rglob("run_info.json"):
        sample = json_file.parent.name
        
        with open(json_file) as f:
            data = json.load(f)
        
        n_processed = data.get("n_processed", 0)
        n_unique = data.get("n_unique", 0)
        n_pseudoaligned = data.get("n_pseudoaligned", 0)
        
        pct_unique = (n_unique / n_processed * 100) if n_processed > 0 else 0
        pct_multi = ((n_pseudoaligned - n_unique) / n_processed * 100) if n_processed > 0 else 0
        
        stats[sample] = {
            "Unique_reads": n_unique,
            "Multi_mapped_reads": n_pseudoaligned - n_unique,
            "Unique_pct": pct_unique,
            "Multi_mapped_pct": pct_multi
        }
    
    return stats

def write_multiqc_file(stats, output_file):
    """Write stats in MultiQC-compatible format"""
    with open(output_file, 'w') as f:
        # Header
        f.write("Sample\tUnique_reads\tMulti_mapped_reads\tUnique_pct\tMulti_mapped_pct\n")
        
        # Data
        for sample, data in sorted(stats.items()):
            f.write(f"{sample}\t{data['Unique_reads']}\t{data['Multi_mapped_reads']}\t"
                   f"{data['Unique_pct']:.2f}\t{data['Multi_mapped_pct']:.2f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: kallisto_unique_summary.py <results_dir> <output_file>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    stats = parse_kallisto_stats(results_dir)
    write_multiqc_file(stats, output_file)
    print(f"Generated {output_file} with {len(stats)} samples")
```

**2.2. Add Process to Workflow**

Location: `workflows/rnaseq/main.nf` (after Kallisto quantification)

```groovy
/*
 * Generate Kallisto unique read summary for MultiQC
 */
process KALLISTO_UNIQUE_SUMMARY {
    label 'process_low'
    
    input:
    path(kallisto_results)
    
    output:
    path("kallisto_unique_summary.txt"), emit: summary
    
    script:
    """
    kallisto_unique_summary.py . kallisto_unique_summary.txt
    """
}

// In workflow section, after KALLISTO_QUANT:
if (!params.skip_pseudo_alignment && params.pseudo_aligner == 'kallisto') {
    KALLISTO_UNIQUE_SUMMARY(
        KALLISTO_QUANT.out.results.collect()
    )
    
    ch_multiqc_files = ch_multiqc_files.mix(
        KALLISTO_UNIQUE_SUMMARY.out.summary.collect().ifEmpty([])
    )
}
```

**2.3. Add MultiQC Configuration**

Add to `workflows/rnaseq/assets/multiqc/multiqc_config.yml`:

```yaml
custom_data:
  kallisto-unique-summary:
    section_name: "Kallisto Unique vs Multi-mapped Reads"
    description: >
      Comparison of uniquely mapped reads (n_unique) versus multi-mapped reads 
      in Kallisto quantification. Higher unique mapping rates indicate better 
      quantification confidence. Low n_unique values may suggest:
      <ul>
        <li>High isoform similarity in the transcriptome</li>
        <li>Repetitive sequences or gene families</li>
        <li>Library complexity issues</li>
        <li>Potential genomic DNA contamination</li>
      </ul>
    section_anchor: kallisto-unique-analysis
    file_format: tsv
    plot_type: bargraph
    pconfig:
      title: "Kallisto: Unique vs Multi-mapped Reads"
      id: kallisto_unique_multimapped
      ylab: "Number of Reads"
      cpswitch_counts_label: "Number of Reads"
      cpswitch_percent_label: "Percentage of Reads"
      yDecimals: false
      tt_percentages: true
      use_legend: true
      data_labels:
        - name: "Read Mapping Classification"
          ylab: "Reads"
    fn: "kallisto_unique_summary.txt"

report_section_order:
  kallisto-unique-analysis:
    order: 1999  # Right after kallisto (2000)
```

---

### 🎯 **Phase 3: Add QC Thresholds** (OPTIONAL)

**3.1. Add Parameter to Control Minimum n_unique**

Location: `nextflow.config`

```groovy
params {
    // Kallisto QC thresholds
    min_unique_reads         = 1000000  // Minimum unique reads
    min_unique_pct           = 50.0     // Minimum % unique of pseudoaligned
}
```

**3.2. Implement QC Check**

```groovy
process KALLISTO_QC_CHECK {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(run_info)
    
    output:
    tuple val(meta), path(run_info), emit: passed
    path("*.fail.txt"), optional: true, emit: failed
    
    script:
    """
    #!/usr/bin/env python3
    import json
    
    with open("${run_info}") as f:
        data = json.load(f)
    
    n_unique = data["n_unique"]
    n_pseudoaligned = data["n_pseudoaligned"]
    pct_unique = (n_unique / n_pseudoaligned * 100) if n_pseudoaligned > 0 else 0
    
    min_reads = ${params.min_unique_reads}
    min_pct = ${params.min_unique_pct}
    
    if n_unique < min_reads or pct_unique < min_pct:
        with open("${meta.id}.fail.txt", "w") as f:
            f.write(f"Sample: ${meta.id}\\n")
            f.write(f"n_unique: {n_unique} (threshold: {min_reads})\\n")
            f.write(f"% unique: {pct_unique:.2f} (threshold: {min_pct})\\n")
            f.write("FAILED QC\\n")
        exit(1)
    """
}
```

---

## Benefits of Implementation

### 📊 **Enhanced QC Reporting**
1. **Visual comparison** of n_unique across all samples
2. **Quick identification** of low-confidence samples
3. **Better understanding** of multi-mapping issues

### 🔍 **Improved Troubleshooting**
1. **Correlate** n_unique with PCA outliers
2. **Identify** samples needing deeper investigation
3. **Track** consistency across batches/runs

### 📈 **Publication-Ready Metrics**
1. **Standard metric** for RNA-seq QC reporting
2. **Transparent** quantification quality assessment
3. **Reproducible** QC thresholds

---

## Example Output

### General Statistics Table (Current)
```
Sample  | Processed | Pseudoaligned | n_unique  | % Pseudoaligned | % Unique
--------|-----------|---------------|-----------|-----------------|----------
Sample1 | 50.0M     | 43.5M (87%)   | 38.5M     | 87.0%           | 77.0%
Sample2 | 48.5M     | 40.2M (83%)   | 25.3M     | 82.9%           | 52.2% ⚠️
Sample3 | 52.1M     | 45.8M (88%)   | 40.1M     | 87.9%           | 77.0%
```

### New Section (Proposed)
**"Kallisto Unique vs Multi-mapped Reads"**
- Bar graph showing unique (green) vs multi-mapped (orange) reads
- Percentage view showing proportion of each
- Sortable by % unique
- Highlight samples below thresholds

---

## Limitations & Considerations

### ⚠️ **Known Challenges**

1. **Multi-isoform Genes**
   - Low n_unique is **expected** for genes with many isoforms
   - Not always a QC failure

2. **Transcriptome Complexity**
   - Some organisms naturally have more ambiguous mappings
   - Thresholds should be organism-specific

3. **Library Type Dependency**
   - Poly-A selected: expect higher n_unique
   - Total RNA: expect lower n_unique (rRNA multi-mapping)

### 💡 **Best Practices**

1. **Set appropriate thresholds** based on organism and library type
2. **Compare within batches** rather than absolute values
3. **Investigate outliers** individually
4. **Cross-reference** with other QC metrics (PCA, sample distance)

---

## Implementation Checklist

- [ ] Phase 1: Verify current n_unique visibility (ALREADY DONE ✅)
- [ ] Phase 2: Decide if enhanced visualization is needed
  - [ ] Create `bin/kallisto_unique_summary.py` script
  - [ ] Add `KALLISTO_UNIQUE_SUMMARY` process to workflow
  - [ ] Update MultiQC config with new custom_data section
  - [ ] Test with sample data
- [ ] Phase 3 (Optional): Implement QC thresholds
  - [ ] Add parameters to `nextflow.config`
  - [ ] Create QC check process
  - [ ] Add to fail_samples reporting

---

## Alternative: Quick Win

### **Minimal Enhancement** (No code changes!)

**Just emphasize existing n_unique in the report:**

Add to `multiqc_config.yml` introduction text:

```yaml
report_comment: >
  <h4>Key Kallisto QC Metric: n_unique</h4>
  <p>Pay special attention to the <b>n_unique</b> column in the General Statistics table. 
  This represents the number of reads that map uniquely to transcripts. Low values may indicate:
  <ul>
    <li>High isoform ambiguity</li>
    <li>Multi-mapping issues</li>
    <li>Library quality concerns</li>
  </ul>
  A good rule of thumb: aim for >70% of pseudoaligned reads being unique.
  </p>
```

---

## Conclusion

✅ **Current Status:** `n_unique` is already visible in General Statistics  
🎯 **Recommendation:** Implement Phase 2 for enhanced visualization  
📊 **Expected Impact:** Better QC insights and easier troubleshooting  
⏱️ **Estimated Effort:** ~2-4 hours for Phase 2 implementation  

The enhanced visualization would make `n_unique` much more prominent and useful for QC decisions, especially when comparing multiple samples or identifying outliers.

---

**Questions? Suggestions?**
Feel free to reach out or open an issue on GitHub!

**Generated by:** Seqera AI  
**Date:** 2025-11-08  
**Workspace:** showcase (ID: 40230138858677)
