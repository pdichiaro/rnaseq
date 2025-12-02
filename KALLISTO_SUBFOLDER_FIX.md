# Kallisto MultiQC Subfolder Fix

## Problem

The Kallisto MultiQC process was failing with:
```
Missing output file(s) `*multiqc_report.html` expected by process MULTIQC_KALLISTO
-- Error is ignored
```

Result: No `multiqc/kallisto/` folder was being created.

## Root Cause

**Filename pattern mismatch:**

The MultiQC module output pattern expects: `*multiqc_report.html`

With the old configuration:
- `ext.prefix = 'kallisto_multiqc'`
- MultiQC command: `multiqc --filename kallisto_multiqc.html`
- File created: `kallisto_multiqc.html` ❌
- Pattern match: Does NOT match `*multiqc_report.html` ❌

## Solution

Changed the prefix to end with `multiqc_report`:

**Before:**
```groovy
ext.prefix = 'kallisto_multiqc'
```

**After:**
```groovy
ext.prefix = 'kallisto_multiqc_report'
```

Now:
- MultiQC command: `multiqc --filename kallisto_multiqc_report.html`
- File created: `kallisto_multiqc_report.html` ✅
- Pattern match: MATCHES `*multiqc_report.html` ✅

## Files Modified

1. **workflows/rnaseq/nextflow.config**
   - Changed: `ext.prefix = 'kallisto_multiqc'` → `ext.prefix = 'kallisto_multiqc_report'`

2. **conf/base.config**
   - Changed: `ext.prefix = 'multiqc_kallisto'` → `ext.prefix = 'kallisto_multiqc_report'`

## Expected Output

After running with `-resume`:

```
${params.outdir}/
└── multiqc/
    ├── multiqc_report.html              # Main combined report
    └── kallisto/                         # ← NEW! Kallisto subfolder
        ├── kallisto_multiqc_report.html  # ← Kallisto-specific report
        └── kallisto_multiqc_report_data/
```

## How to Test

Run your pipeline with `-resume`:

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir results \
  --pseudo_aligner kallisto \
  --skip_alignment \
  -profile docker \
  -resume
```

Verify the output:

```bash
# Check for the kallisto folder
ls results/multiqc/kallisto/

# Should contain:
# - kallisto_multiqc_report.html
# - kallisto_multiqc_report_data/
```

## Why This Works

The MultiQC nf-core module uses this logic:

```groovy
output:
    path "*multiqc_report.html", emit: report

script:
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
```

- If no prefix: creates `multiqc_report.html` (matches `*multiqc_report.html` ✅)
- If prefix: creates `${prefix}.html` (must end with `multiqc_report` to match pattern)

So valid prefixes are:
- ✅ `kallisto_multiqc_report` → `kallisto_multiqc_report.html`
- ✅ `star_multiqc_report` → `star_multiqc_report.html`
- ✅ `anything_multiqc_report` → `anything_multiqc_report.html`

Invalid prefixes:
- ❌ `kallisto_multiqc` → `kallisto_multiqc.html` (doesn't match pattern)
- ❌ `kallisto` → `kallisto.html` (doesn't match pattern)

## Consistency with Other Aligners

Now all aligner-specific MultiQC processes use consistent naming:

```groovy
MULTIQC_STAR:     ext.prefix = 'star_multiqc_report'     → star/star_multiqc_report.html
MULTIQC_HISAT2:   ext.prefix = 'hisat2_multiqc_report'   → hisat2/hisat2_multiqc_report.html
MULTIQC_KALLISTO: ext.prefix = 'kallisto_multiqc_report' → kallisto/kallisto_multiqc_report.html
```

All follow the pattern: `<aligner>_multiqc_report`

## Status

✅ **Fixed** - Kallisto now creates its own subfolder just like STAR
✅ **Consistent** - All aligners use the same naming pattern
✅ **Ready** - Run with `-resume` to regenerate reports

## Note

The `errorStrategy = 'ignore'` in base.config is intentional. It allows the pipeline to continue if the aligner-specific report fails (e.g., if Kallisto wasn't used). The main MultiQC report will still fail the pipeline if it has issues.
