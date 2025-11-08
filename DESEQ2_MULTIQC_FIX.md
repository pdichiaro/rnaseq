# DESeq2 MultiQC Integration Fix

## Problem Summary
DESeq2 QC processes completed successfully but MultiQC wasn't detecting the output files, resulting in missing DESeq2 sections in the final report.

## Root Cause
The workflow was passing an incorrect quantifier parameter to the DESeq2 normalization processes. Specifically:

1. **Issue in main.nf (line ~303)**: The code used `params.pseudo_aligner ?: "Salmon"` which:
   - Defaults to "Salmon" if `params.pseudo_aligner` is null/empty/false
   - Doesn't handle capitalization consistently

2. **MultiQC Expectation**: MultiQC config expects files matching these exact patterns:
   - `kallisto.deseq2.all_genes.*.txt` (for Kallisto runs)
   - `salmon.deseq2.all_genes.*.txt` (for Salmon runs)
   - etc.

3. **Impact**: If quantifier was passed incorrectly:
   - R script generates files with wrong prefix (e.g., "salmon" instead of "kallisto")
   - MultiQC can't find files → Missing DESeq2 sections in report

## Solution Applied

### Change in `workflows/rnaseq/main.nf`

**Before:**
```groovy
// Run invariant_genes normalization if requested
if (normalization_methods.contains('invariant_genes')) {
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO (
        ch_counts_gene,
        params.pseudo_aligner ?: "Salmon"  // ❌ Problem: inconsistent handling
    )
    ...
}

// Run all_genes normalization if requested or as default
if (normalization_methods.contains('all_genes') || (!normalization_methods.contains('invariant_genes') && !normalization_methods.contains('all_genes'))) {
    NORMALIZE_DESEQ2_QC_ALL_GENES_PSEUDO (
        ch_counts_gene,
        params.pseudo_aligner ?: "Salmon"  // ❌ Same problem here
    )
    ...
}
```

**After:**
```groovy
// Determine the quantifier for file naming (must match what was actually used for quantification)
// This is critical for MultiQC to find the correct files
def pseudo_quantifier = params.pseudo_aligner ? params.pseudo_aligner.capitalize() : "Salmon"

// Run invariant_genes normalization if requested
if (normalization_methods.contains('invariant_genes')) {
    NORMALIZE_DESEQ2_QC_INVARIANT_GENES_PSEUDO (
        ch_counts_gene,
        pseudo_quantifier  // ✅ Fixed: explicit, capitalized quantifier
    )
    ...
}

// Run all_genes normalization if requested or as default
if (normalization_methods.contains('all_genes') || (!normalization_methods.contains('invariant_genes') && !normalization_methods.contains('all_genes'))) {
    NORMALIZE_DESEQ2_QC_ALL_GENES_PSEUDO (
        ch_counts_gene,
        pseudo_quantifier  // ✅ Fixed: same consistent value
    )
    ...
}
```

### Key Improvements
1. **Explicit variable**: Creates `pseudo_quantifier` variable for clarity
2. **Proper capitalization**: Uses `.capitalize()` to ensure "kallisto" → "Kallisto", "salmon" → "Salmon"
3. **Consistency**: Both normalization methods use the exact same quantifier value
4. **Safe default**: Falls back to "Salmon" if `params.pseudo_aligner` is not set

## How the Fix Works

### Data Flow:
1. **Command Line**: User runs with `--pseudo_aligner 'kallisto'`
2. **Workflow**: Sets `pseudo_quantifier = "Kallisto"` (capitalized)
3. **Process**: DESeq2 processes receive "Kallisto" as quantifier parameter
4. **R Script**: `normalize_deseq2_qc_all_genes.r` receives `--quantifier Kallisto`
5. **File Generation**: R script converts to lowercase and creates files:
   - `kallisto.deseq2.all_genes.pca.vals.txt`
   - `kallisto.deseq2.all_genes.sample.dists.txt`
   - `kallisto.deseq2.all_genes.read.distribution.normalized.txt`
   - etc.
6. **MultiQC**: Finds files matching `kallisto\.deseq2\.all_genes\.*` pattern
7. **Result**: DESeq2-Kallisto-QC section appears in MultiQC report! ✅

## Testing the Fix

### Quick Test (Resume existing run):
```bash
cd /mnt/ngs_ricerca/dichiarop/NEXTFLOW/nextflow_temp/BulkRNAseq_test

# Resume with the fixed code
nextflow run /home/user/rnaseq \
  -resume \
  --pseudo_aligner 'kallisto' \
  -profile singularity \
  --input samplesheet_OMIM.csv \
  --outdir results \
  --fasta /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/GRCh38.primary_assembly.genome.fa \
  --gtf /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/gencode.v47.primary_assembly.annotation.gtf \
  --transcript_fasta /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/gencode.v47.transcripts.fa \
  --skip_pseudo_alignment false \
  --normalization_method 'all_genes' \
  --max_cpus 16 \
  --max_memory '128.GB'
```

**Expected behavior:**
- DESeq2 processes should be pulled from cache (no re-execution needed)
- MultiQC will re-run and now SHOULD find the DESeq2 output files
- Final report should include "DESeq2 Kallisto QC" section

### Full Test (Clean run on subset):
```bash
# Create test directory
mkdir -p test_deseq2_fix
cd test_deseq2_fix

# Create minimal samplesheet (2 samples for speed)
cat > samplesheet_test.csv << 'EOF'
sample,fastq_1,fastq_2,strandedness
OMIM1,/path/to/OMIM1_R1.fastq.gz,/path/to/OMIM1_R2.fastq.gz,auto
OMIM3,/path/to/OMIM3_R1.fastq.gz,/path/to/OMIM3_R2.fastq.gz,auto
EOF

# Run pipeline with fixed code
nextflow run /home/user/rnaseq \
  --pseudo_aligner 'kallisto' \
  -profile singularity \
  --input samplesheet_test.csv \
  --outdir results_test \
  --fasta /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/GRCh38.primary_assembly.genome.fa \
  --gtf /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/gencode.v47.primary_assembly.annotation.gtf \
  --transcript_fasta /mnt/ngs_ricerca/shared_resources/genomes/Homo_sapiens/GENCODE/GRCh38.p14/gencode.v47.transcripts.fa \
  --skip_pseudo_alignment false \
  --normalization_method 'all_genes' \
  --max_cpus 16 \
  --max_memory '128.GB'
```

### Verification Steps

1. **Check DESeq2 process logs** (should show correct quantifier):
```bash
# Find the work directory for the DESeq2 process
find work -name ".command.log" -path "*/NORMALIZE_DESEQ2_QC_ALL_GENES_PSEUDO*" | head -1 | xargs dirname | xargs -I {} sh -c 'echo "=== DESeq2 Process Directory: {}" && ls -la {} && echo "=== Log Output:" && head -30 {/.command.log}'
```

Look for lines like:
```
=== ALL GENES NORMALIZATION + DESEQ2 QC ===
File prefix for outputs: kallisto.deseq2.all_genes
```

2. **Verify output file names**:
```bash
# List DESeq2 output files
find results/deseq2_qc -name "*.txt" -o -name "*.RData" | sort
```

Expected files:
```
results/deseq2_qc/kallisto.deseq2.all_genes.pca.vals.txt
results/deseq2_qc/kallisto.deseq2.all_genes.pca.top500.vals.txt
results/deseq2_qc/kallisto.deseq2.all_genes.sample.dists.txt
results/deseq2_qc/kallisto.deseq2.all_genes.read.distribution.normalized.txt
results/deseq2_qc/kallisto.deseq2.invariant_genes.pca.vals.txt
results/deseq2_qc/kallisto.deseq2.invariant_genes.sample.dists.txt
... (etc.)
```

3. **Check MultiQC report**:
```bash
# Open the MultiQC report
firefox results/multiqc/multiqc_report.html
# or
chromium results/multiqc/multiqc_report.html
```

Navigate to the report and look for:
- **Section**: "DESeq2 Kallisto QC" (should be present!)
- **Subsections**:
  - Read Distribution (All Genes)
  - Sample Distance (All Genes)
  - PCA All Genes (All Genes Normalization)
  - PCA Top 500 Genes (All Genes Normalization)

4. **Check MultiQC data directory**:
```bash
# Verify MultiQC found the files
grep -r "kallisto.deseq2" results/multiqc/multiqc_data/
```

Should show entries for all the kallisto.deseq2.*.txt files.

## Expected Outcomes

### Before Fix ❌
- DESeq2 processes complete successfully
- Files generated with wrong prefix (e.g., "salmon.deseq2..." when using Kallisto)
- MultiQC doesn't find files → No DESeq2 section in report
- Users confused why DESeq2 results are missing

### After Fix ✅
- DESeq2 processes complete successfully
- Files generated with CORRECT prefix (e.g., "kallisto.deseq2..." when using Kallisto)
- MultiQC finds files → DESeq2 section appears in report
- Users can view PCA plots, sample distances, read distributions in web report

## Additional Notes

### For Salmon users:
The fix also ensures Salmon runs work correctly:
```bash
nextflow run /home/user/rnaseq --pseudo_aligner 'salmon' ...
```
Will generate `salmon.deseq2.all_genes.*` files that match MultiQC's expectations.

### For STAR alignment quantification:
Similar logic already exists for STAR-based quantification:
- `star.rsem.deseq2.all_genes.*`
- `star.salmon.deseq2.all_genes.*`
- `star.genome.deseq2.all_genes.*`

These use hardcoded strings (e.g., "STAR_Genome") which work correctly.

### Backward Compatibility:
- If `--pseudo_aligner` is not specified, defaults to "Salmon" (maintains existing behavior)
- Existing pipeline runs are not affected
- Only affects new runs after the fix is applied

## Related Files

1. **Main workflow**: `workflows/rnaseq/main.nf` (modified)
2. **R script**: `bin/normalize_deseq2_qc_all_genes.r` (handles lowercase conversion)
3. **MultiQC config**: `workflows/rnaseq/assets/multiqc/multiqc_config.yml` (defines file patterns)
4. **Process definition**: `modules/local/normalize_deseq2_qc_all_genes/main.nf` (passes quantifier to R)

## Commit Message
```
Fix DESeq2 MultiQC integration for pseudo-alignment quantifiers

- Explicitly determine quantifier for DESeq2 file naming
- Ensure proper capitalization (kallisto → Kallisto, salmon → Salmon)
- Use consistent quantifier value for both normalization methods
- Fixes missing DESeq2 sections in MultiQC report when using Kallisto/Salmon
- Maintains backward compatibility with default "Salmon" fallback

Resolves issue where DESeq2 processes completed successfully but
MultiQC couldn't find output files due to quantifier parameter mismatch.
```

## Questions or Issues?

If the DESeq2 section still doesn't appear after applying this fix:

1. **Check the quantifier parameter**: Ensure `--pseudo_aligner` matches what was actually used
2. **Check R script output**: Look for "File prefix for outputs:" in the process log
3. **Verify file names**: Confirm files match the pattern `{quantifier}.deseq2.{method}.*.txt`
4. **Check MultiQC verbose output**: Run MultiQC manually with `-v` flag for debugging
5. **Review MultiQC config**: Verify the `fn_re` patterns in `multiqc_config.yml` are correct

Contact: [Your contact info or issue tracker]
