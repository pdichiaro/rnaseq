# Understanding Channel Structure for DESeq2 Outputs

## The Core Question

When a process outputs `path "*.pca.vals.txt"`, what is the channel structure?

### Case 1: Single File Match
```groovy
output:
path "*.pca.vals.txt", emit: pca_txt
```

If the process produces: `sample1.pca.vals.txt`

**Channel structure:**
```
pca_txt = [ file('sample1.pca.vals.txt') ]
```

###Case 2: Multiple Files Match
```groovy
output:
path "*.pca.vals.txt", emit: pca_txt
```

If the process produces multiple files:
- `sample1.pca.vals.txt`
- `sample2.pca.vals.txt`

**Channel structure (CRITICAL):**
```
pca_txt = [ 
  file('sample1.pca.vals.txt'),
  file('sample2.pca.vals.txt')
]
```

**OR potentially:**
```
pca_txt = [
  [file('sample1.pca.vals.txt'), file('sample2.pca.vals.txt')]
]
```

## The Real Question: Does DESeq2 Emit Multiple Files?

Looking at the process:
```groovy
process NORMALIZE_DESEQ2_QC_INVARIANT_GENES {
    input:
    path counts
    val quantifier

    output:
    path "*.pca.vals.txt", emit: pca_all_genes_txt
```

### In the Script Section:
The R script generates files like:
- `{prefix}.pca.vals.txt` (ONE file per process run)

So each DESeq2 process run produces **ONE file per output**, not multiple!

## Comparing with Tuple Outputs

### FastQC (Tuple Output):
```groovy
process FASTQC {
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.zip"), emit: zip
}
```

**Channel structure:**
```
zip = [
  [meta1, file1.zip],
  [meta2, file2.zip],
  [meta3, file3.zip]
]
```

To extract just files: `.collect{it[1]}` gives `[file1.zip, file2.zip, file3.zip]`

### DESeq2 (Path-only Output):
```groovy
process NORMALIZE_DESEQ2_QC {
    output:
    path "*.pca.vals.txt", emit: pca_all_genes_txt
}
```

**Channel structure:**
```
pca_all_genes_txt = file('kallisto.deseq2.pca.vals.txt')
```

**No tuple, no meta, just the file!**

## The Mix Operation

### When mixing path-only outputs:
```groovy
ch_multiqc_files = Channel.empty()
ch_multiqc_files = ch_multiqc_files.mix(PROCESS1.out.file)  // Adds: file1
ch_multiqc_files = ch_multiqc_files.mix(PROCESS2.out.file)  // Adds: file2  
ch_multiqc_files = ch_multiqc_files.mix(PROCESS3.out.file)  // Adds: file3
```

**Result:**
```
ch_multiqc_files = [file1, file2, file3]
```

This should already be flat!

## So Why Did .flatten() Help?

The answer might be that some DESeq2 outputs emit **ARRAYS OF FILES** even though it's a single pattern match!

Let me check the actual glob pattern behavior...

### If the Pattern Matches Multiple Files in ONE Process Run:
```bash
# In process directory:
ls *.pca.vals.txt
# Output:
# kallisto.deseq2.invariant.pca.vals.txt
# kallisto.deseq2.allgenes.pca.vals.txt  # IF multiple files exist
```

Then the output would be:
```
pca_all_genes_txt = [
  [file1.txt, file2.txt]  # ← Nested array!
]
```

## Hypothesis: The Glob Pattern Matches Multiple Files

Looking at the output patterns:
```groovy
path "*.pca.vals.txt"       , emit: pca_all_genes_txt
path "*.pca.top*.vals.txt"  , emit: pca_top_genes_txt
path "*.sample.dists.txt"   , emit: sample_distances_txt
```

These patterns could potentially match multiple files if the R script generates multiple files!

## Testing Theory: Check R Script Output

We need to verify what files the R script actually generates. From the stub:
```groovy
stub:
touch ${prefix}.pca.vals.txt
touch ${prefix}.pca.top500.vals.txt
touch ${prefix}.sample.dists.txt
```

So with `prefix = "kallisto.deseq2.invariant"`, we get:
- `kallisto.deseq2.invariant.pca.vals.txt` ← matches `*.pca.vals.txt`
- `kallisto.deseq2.invariant.pca.top500.vals.txt` ← matches `*.pca.top*.vals.txt`
- `kallisto.deseq2.invariant.sample.dists.txt` ← matches `*.sample.dists.txt`

Each pattern should match **ONE file per process run**.

## The Real Issue Must Be Different!

If each pattern matches only ONE file, and we're mixing path-only outputs, the channel should be flat.

**UNLESS...**

What if the process is called **MULTIPLE TIMES** (for different samples), and `.collect()` was needed to gather all instances?

Let me check how the process is called...
