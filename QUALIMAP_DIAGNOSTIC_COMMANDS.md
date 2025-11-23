# Qualimap Diagnostic Commands

If Qualimap folder is empty even with `--min_mapped_reads 0`, run these diagnostic commands:

## 1. Check if Qualimap Tasks Were Executed

```bash
# Search for Qualimap in Nextflow log
grep -i "qualimap" .nextflow.log | head -20

# Check if Qualimap process was submitted
grep "QUALIMAP_RNASEQ" .nextflow.log | wc -l

# Look for Qualimap in work directory
find work -name "qualimap*" -o -name "*.html" | grep -i qualimap
```

## 2. Check BAM Files Available for Qualimap

```bash
# Check if BAM files exist in results
ls -lh results/star/*.bam 2>/dev/null || echo "No BAM files in results/star/"
ls -lh results/hisat2/*.bam 2>/dev/null || echo "No BAM files in results/hisat2/"

# Check BAM files in work directory
find work -name "*.bam" -type f | head -10

# Check if genome BAM channel had content
grep "ch_genome_bam" .nextflow.log | head -20
```

## 3. Check Qualimap Configuration

```bash
# Verify skip_qualimap is false
grep "skip_qualimap" nextflow.config

# Check if skip_qc is disabled
grep "skip_qc" nextflow.config

# Look for any Qualimap-specific config
grep -i "qualimap" nextflow.config
```

## 4. Check Qualimap Work Directory

```bash
# Find Qualimap work directories
find work -type d -name "*QUALIMAP*"

# If found, check their contents
QUALIMAP_WORK=$(find work -type d -name "*QUALIMAP*" | head -1)
if [ -n "$QUALIMAP_WORK" ]; then
    echo "Found Qualimap work dir: $QUALIMAP_WORK"
    ls -lha "$QUALIMAP_WORK"
    cat "$QUALIMAP_WORK/.command.log" 2>/dev/null || echo "No command log found"
fi
```

## 5. Check Nextflow Execution Trace

```bash
# Look at execution trace (if enabled)
cat pipeline_info/execution_trace*.txt | grep -i qualimap

# Check timeline (if enabled)
cat pipeline_info/execution_timeline*.html | grep -i qualimap
```

## 6. Check Output Directory Structure

```bash
# List all directories in output
find results -type d | sort

# Check if qualimap folder exists but is empty
ls -la results/*/qualimap/ 2>/dev/null || echo "No qualimap folder found"
ls -la results/qualimap/ 2>/dev/null || echo "No top-level qualimap folder"
```

## 7. Check for Qualimap Errors

```bash
# Search for errors related to Qualimap
grep -i "qualimap.*error\|qualimap.*fail" .nextflow.log

# Check stderr from Qualimap processes
find work -name ".command.err" -path "*/QUALIMAP*" -exec echo "=== {} ===" \; -exec cat {} \;
```

## 8. Check Process Completion Status

```bash
# Look for Qualimap process status
grep "QUALIMAP_RNASEQ.*COMPLETED\|QUALIMAP_RNASEQ.*CACHED" .nextflow.log
grep "QUALIMAP_RNASEQ.*FAILED\|QUALIMAP_RNASEQ.*ERROR" .nextflow.log
```

## 9. Check Channel Contents (Advanced)

If you have access to Nextflow tower or reports:

```bash
# Check if reports were generated
ls -lh results/pipeline_info/

# Look at DAG to see if Qualimap was included
cat results/pipeline_info/pipeline_dag*.html 2>/dev/null
```

## 10. Manual Qualimap Test

Test if Qualimap would work manually:

```bash
# Find a BAM file
BAM=$(find results -name "*.bam" -type f | head -1)
echo "Testing with BAM: $BAM"

# Find GTF
GTF=$(find . -name "*.gtf" -o -name "*.gtf.gz" | head -1)
echo "Using GTF: $GTF"

# Try running Qualimap manually
if [ -f "$BAM" ] && [ -f "$GTF" ]; then
    qualimap rnaseq \
        -bam "$BAM" \
        -gtf "$GTF" \
        -outdir test_qualimap_manual \
        --java-mem-size=4G
fi
```

## Common Issues & Interpretations

### Scenario 1: No Qualimap Tasks Found
```
$ grep "QUALIMAP" .nextflow.log | wc -l
0
```
**Cause**: Qualimap is being skipped
**Check**: 
- Is `--skip_qualimap true` in your command?
- Is `--skip_qc true` in your command?
- Is `ch_genome_bam` empty before Qualimap runs?

### Scenario 2: Qualimap Tasks Failed
```
$ grep "QUALIMAP_RNASEQ.*FAILED" .nextflow.log
[...] QUALIMAP_RNASEQ (sample1) FAILED
```
**Cause**: Qualimap execution error
**Check**: Error logs in work directory

### Scenario 3: Qualimap Completed But No Output
```
$ grep "QUALIMAP_RNASEQ.*COMPLETED" .nextflow.log
[...] QUALIMAP_RNASEQ (sample1) COMPLETED
```
**Cause**: PublishDir issue or output in unexpected location
**Check**: Work directory for actual outputs

### Scenario 4: Qualimap Not Run Due to Empty Channel
```
$ grep "ch_genome_bam" .nextflow.log
[...] ch_genome_bam is empty
```
**Cause**: No BAM files in ch_genome_bam channel
**Reasons**:
- Using pseudoaligner only (kallisto) without genome alignment
- All BAM files filtered out (even with min_mapped_reads=0)
- Alignment failed

## Quick Diagnostic Script

Save as `diagnose_qualimap.sh` and run:

```bash
#!/bin/bash

echo "=== Qualimap Diagnostic Report ==="
echo ""

echo "1. Checking if Qualimap tasks were executed..."
QUALIMAP_COUNT=$(grep -c "QUALIMAP_RNASEQ" .nextflow.log 2>/dev/null || echo "0")
echo "   Qualimap mentions in log: $QUALIMAP_COUNT"

echo ""
echo "2. Checking for Qualimap work directories..."
WORK_DIRS=$(find work -type d -name "*QUALIMAP*" 2>/dev/null | wc -l)
echo "   Qualimap work directories found: $WORK_DIRS"

echo ""
echo "3. Checking for BAM files..."
BAM_COUNT=$(find results -name "*.bam" 2>/dev/null | wc -l)
echo "   BAM files in results: $BAM_COUNT"

echo ""
echo "4. Checking configuration..."
echo "   skip_qualimap: $(grep skip_qualimap nextflow.config 2>/dev/null || echo 'not found')"
echo "   skip_qc: $(grep skip_qc nextflow.config 2>/dev/null || echo 'not found')"

echo ""
echo "5. Checking output directories..."
if [ -d "results/qualimap" ]; then
    QUALIMAP_FILES=$(find results/qualimap -type f 2>/dev/null | wc -l)
    echo "   Qualimap folder exists with $QUALIMAP_FILES files"
elif [ -d "results/star/qualimap" ]; then
    QUALIMAP_FILES=$(find results/star/qualimap -type f 2>/dev/null | wc -l)
    echo "   Qualimap folder exists at results/star/qualimap with $QUALIMAP_FILES files"
else
    echo "   No qualimap output folder found"
fi

echo ""
echo "6. Checking for Qualimap failures..."
FAILURES=$(grep "QUALIMAP.*FAIL\|QUALIMAP.*ERROR" .nextflow.log 2>/dev/null | wc -l)
echo "   Qualimap failures: $FAILURES"

if [ $FAILURES -gt 0 ]; then
    echo ""
    echo "   Error details:"
    grep "QUALIMAP.*FAIL\|QUALIMAP.*ERROR" .nextflow.log 2>/dev/null | head -5
fi

echo ""
echo "=== End Diagnostic Report ==="
```

## Next Steps Based on Results

1. **If no Qualimap tasks executed**: Check if `ch_genome_bam` is populated
2. **If tasks failed**: Check error logs in work directory
3. **If tasks completed**: Check work directory for outputs and publishDir configuration
4. **If no BAM files**: Check alignment configuration and logs

## Key Question to Answer

**Is Qualimap running at all, or is it being skipped?**

Run:
```bash
grep -i "qualimap" .nextflow.log | grep -i "submit\|start\|complete\|cache"
```

This will show you if Qualimap tasks were attempted.
