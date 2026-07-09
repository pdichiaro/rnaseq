# RNA-seq Alignment Filtering Summary

## STAR and HISAT2 filtering workflow

This document summarizes how aligned reads are retained or excluded after genome alignment in the RNA-seq pipeline when using **STAR** or **HISAT2**.

The key point is that the two aligners now use a more harmonized BAM filtering strategy for **genome-based quantification**:

- **STAR**: multimapping is controlled at the STAR alignment level through `--outFilterMultimapNmax`. In addition, when `params.quantification == 'genome'`, the STAR genome BAM is filtered after alignment with `samtools view`, using the same flag logic adopted for HISAT2. This filter is not applied to RSEM/transcriptome-based quantification.
- **HISAT2**: HISAT2 output is piped directly into `samtools view`, where unmapped reads, mate-unmapped reads in paired-end mode, secondary alignments and, in paired-end mode, non-proper pairs are removed before BAM sorting.

Important distinction:

```text
proper-pair filtering ≠ unique-mapping filtering
```

The `samtools` filters described here clean the BAM at the alignment-flag level. They do **not** guarantee unique-only mapping. Unique/multimapping behavior is controlled separately by the aligner-specific settings.

---

## Complete filtering workflow: STAR

```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃              STAR ALIGNMENT AND BAM PROCESSING                  ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

📁 INPUT: sample_R1.fq.gz (+ sample_R2.fq.gz if paired-end)
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 1: STAR ALIGNMENT                                           │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  STAR                                                             │
│    --quantMode TranscriptomeSAM                                   │
│    --outSAMtype BAM Unsorted                                      │
│    --outSAMattributes NH HI AS NM MD                              │
│    --readFilesCommand zcat                                        │
│    --outFilterMultimapNmax 20                                     │
│                                                                   │
│  Standard/genome mode additionally uses:                          │
│    --twopassMode Basic                                            │
│    --runRNGseed 0                                                 │
│    --alignSJDBoverhangMin 1                                       │
│    --outSAMstrandField intronMotif                                │
│                                                                   │
│  RSEM-compatible mode additionally uses:                          │
│    --outSAMunmapped Within                                        │
│    --outFilterType BySJout                                        │
│    --outFilterMismatchNmax 999                                    │
│    --outFilterMismatchNoverLmax 0.04                              │
│    --alignIntronMin 20                                            │
│    --alignIntronMax 1000000                                       │
│    --alignMatesGapMax 1000000                                     │
│    --alignSJoverhangMin 8                                         │
│    --alignSJDBoverhangMin 1                                       │
│    --sjdbScore 1                                                  │
│                                                                   │
│  Output: Aligned.out.bam / genome BAM + transcriptome BAM          │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    STAR genome BAM
    ├─ BAM is generated directly by STAR
    ├─ Multimapping is controlled by --outFilterMultimapNmax
    ├─ Default allows multimappers up to 20 loci
    └─ Optional unique-only mode can be set with --outFilterMultimapNmax 1
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 2: CONDITIONAL STAR BAM FILTERING                           │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  Applied only when:                                               │
│                                                                   │
│    params.quantification == 'genome'                              │
│                                                                   │
│  Single-end filtering:                                            │
│                                                                   │
│    samtools view -bS -F 4 -F 256 \                                │
│      sample.Aligned.out.bam > sample.Aligned.filtered.out.bam      │
│                                                                   │
│  Removes:                                                         │
│    -F 4     unmapped reads                                        │
│    -F 256   secondary alignments                                  │
│                                                                   │
│  Paired-end filtering:                                            │
│                                                                   │
│    samtools view -bS -f 2 -F 4 -F 8 -F 256 \                      │
│      sample.Aligned.out.bam > sample.Aligned.filtered.out.bam      │
│                                                                   │
│  Keeps/removes:                                                   │
│    -f 2     keep only proper pairs                                │
│    -F 4     remove unmapped reads                                 │
│    -F 8     remove reads with unmapped mate                       │
│    -F 256   remove secondary alignments                           │
│                                                                   │
│  The filtered BAM replaces the original Aligned.out.bam name       │
│  so that downstream steps receive the filtered genome BAM.         │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    filtered STAR genome BAM
    ├─ No unmapped reads
    ├─ No secondary alignments
    ├─ No mate-unmapped reads in paired-end mode
    ├─ Proper pairs only in paired-end mode
    └─ Multimappers may still be present unless STAR is run unique-only
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 3: BAM_SORT_STATS_SAMTOOLS                                  │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  samtools sort                                                    │
│  samtools index                                                   │
│  samtools stats                                                   │
│  samtools flagstat                                                │
│  samtools idxstats                                                │
│                                                                   │
│  Output: sorted/indexed genome BAM + QC metrics                    │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    sorted STAR BAM
    ├─ Sorted and indexed
    ├─ QC metrics generated
    └─ For genome quantification, this BAM derives from the filtered BAM
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 4: OPTIONAL DEDUPLICATION / DUPLICATE MARKING               │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  If params.with_umi = true:                                       │
│    BAM_DEDUP_UMI_STAR                                             │
│                                                                   │
│  Else if !params.skip_markduplicates:                             │
│    Picard MarkDuplicates                                          │
│                                                                   │
│  Note: Picard MarkDuplicates marks duplicate reads by default.     │
│  It does not remove them unless extra Picard arguments are used.   │
└───────────────────────────────────────────────────────────────────┘
```

---

## STAR filtering interpretation

### Main STAR filtering parameters

| Parameter | Meaning | Effect |
|----------|---------|--------|
| `--outFilterMultimapNmax 20` | Maximum number of loci a read can map to | Default STAR behavior in this pipeline; reads mapping to more than 20 loci are excluded by STAR |
| `--outFilterMultimapNmax 1` | Strict unique-mapping mode | Optional setting; only uniquely mapping reads are reported |
| `--outFilterType BySJout` | Splice-junction-aware filtering | Used in RSEM-compatible mode |
| `--outFilterMismatchNmax 999` | Absolute mismatch limit | Effectively permissive absolute mismatch threshold |
| `--outFilterMismatchNoverLmax 0.04` | Mismatch fraction threshold | In RSEM-compatible mode, reads with mismatch fraction above 4% are excluded |
| `--outSAMtype BAM Unsorted` | BAM output | STAR writes an unsorted BAM directly |
| `--outSAMattributes NH HI AS NM MD` | SAM tags | Keeps multimapping and alignment score/mismatch tags useful for downstream QC |

### STAR `extra_star_align_args`

STAR multimapping behavior can be changed through:

```bash
--extra_star_align_args "--outFilterMultimapNmax 1"
```

This overrides the default multimapping threshold and runs STAR in strict unique-mapping mode.

Default behavior:

```bash
--outFilterMultimapNmax 20
```

Strict unique-only behavior:

```bash
--outFilterMultimapNmax 1
```

### STAR post-alignment BAM filtering

For genome-based quantification only, the STAR BAM is filtered after alignment.

Single-end:

```bash
samtools view -bS -F 4 -F 256 sample.Aligned.out.bam > sample.Aligned.filtered.out.bam
mv sample.Aligned.filtered.out.bam sample.Aligned.out.bam
```

Paired-end:

```bash
samtools view -bS -f 2 -F 4 -F 8 -F 256 sample.Aligned.out.bam > sample.Aligned.filtered.out.bam
mv sample.Aligned.filtered.out.bam sample.Aligned.out.bam
```

This post-alignment filter:

- removes unmapped reads;
- removes secondary alignments;
- removes mate-unmapped reads in paired-end mode;
- retains only proper pairs in paired-end mode.

This filter does **not** remove multimappers by itself. Multimappers are controlled by `--outFilterMultimapNmax`.

### RSEM-specific note

When `params.quantification != 'genome'`, for example in RSEM/transcriptome-based quantification, this STAR BAM filter is not applied.

This is intentional: RSEM uses multimapping and ambiguous reads probabilistically for abundance estimation. Therefore, unique-only or proper-pair filtering should not be applied before RSEM quantification.

---

## Important note on STAR `flagstat`

For STAR-aligned BAM files, `samtools flagstat` should not be used as the main estimate of the original FASTQ mapping rate.

If unmapped reads are not written into the BAM, `flagstat` may report nearly 100% mapped reads because it only sees alignments present in the BAM. The correct STAR mapping summary should be read from `Log.final.out`, especially:

- `Number of input reads`
- `Uniquely mapped reads number`
- `Uniquely mapped reads %`
- `Number of reads mapped to multiple loci`
- `% of reads mapped to multiple loci`
- `Number of reads mapped to too many loci`
- `% of reads mapped to too many loci`

For paired-end data, STAR reports input reads as read pairs/fragments, while `samtools flagstat` counts BAM records. Therefore, 52M input read pairs in STAR may correspond to about 104M BAM records in `flagstat`.

---

## Complete filtering workflow: HISAT2

```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃             HISAT2 ALIGNMENT AND BAM FILTERING                  ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

📁 INPUT: sample_R1.fq.gz (+ sample_R2.fq.gz if paired-end)
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 1: HISAT2 ALIGNMENT                                         │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  Single-end:                                                      │
│    hisat2 -x genome -U reads.fq.gz                                │
│                                                                   │
│  Paired-end:                                                      │
│    hisat2 -x genome -1 reads_1.fq.gz -2 reads_2.fq.gz              │
│      --no-mixed                                                   │
│      --no-discordant                                              │
│                                                                   │
│  Default args:                                                    │
│    --met-stderr --new-summary --dta                               │
│                                                                   │
│  Optional extra args:                                             │
│    params.extra_hisat2_align_args                                 │
│                                                                   │
│  Output: SAM stream                                               │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    SAM stream from HISAT2
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 2: SAMTOOLS VIEW FILTERING                                  │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  Single-end filtering:                                            │
│                                                                   │
│    samtools view -bS -F 4 -F 256 - > sample.bam                   │
│                                                                   │
│  Removes:                                                         │
│    -F 4     unmapped reads                                        │
│    -F 256   secondary alignments                                  │
│                                                                   │
│  Paired-end filtering:                                            │
│                                                                   │
│    samtools view -bS -f 2 -F 4 -F 8 -F 256 - > sample.bam          │
│                                                                   │
│  Keeps/removes:                                                   │
│    -f 2     keep only proper pairs                                │
│    -F 4     remove unmapped reads                                 │
│    -F 8     remove reads with unmapped mate                       │
│    -F 256   remove secondary alignments                           │
│                                                                   │
│  Additional HISAT2 paired-end constraints:                        │
│    --no-mixed       suppress unpaired/mixed alignments             │
│    --no-discordant  suppress discordant alignments                 │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    filtered HISAT2 BAM
    ├─ No unmapped reads
    ├─ No secondary alignments
    ├─ No mate-unmapped reads in paired-end mode
    ├─ Proper pairs only in paired-end mode
    ├─ No mixed alignments in paired-end mode
    └─ No discordant alignments in paired-end mode
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 3: BAM_SORT_STATS_SAMTOOLS                                  │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  samtools sort                                                    │
│  samtools index                                                   │
│  samtools stats                                                   │
│  samtools flagstat                                                │
│  samtools idxstats                                                │
│                                                                   │
│  Output: sorted/indexed genome BAM + QC metrics                    │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    sorted HISAT2 BAM
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 4: OPTIONAL DEDUPLICATION / DUPLICATE MARKING               │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  If params.with_umi = true:                                       │
│    BAM_DEDUP_UMI_HISAT2                                           │
│                                                                   │
│  Else if !params.skip_markduplicates:                             │
│    Picard MarkDuplicates                                          │
│                                                                   │
│  Note: Picard MarkDuplicates marks duplicate reads by default.     │
│  It does not remove them unless extra Picard arguments are used.   │
└───────────────────────────────────────────────────────────────────┘
```

---

## HISAT2 filtering interpretation

### Single-end HISAT2 command

```bash
hisat2 \
    -x $INDEX \
    -U reads.fq.gz \
    --summary-file sample.hisat2.summary.log \
    --threads $task.cpus \
    $seq_center \
    $unaligned \
    $args \
| samtools view -bS -F 4 -F 256 - > sample.bam
```

### Paired-end HISAT2 command

```bash
hisat2 \
    -x $INDEX \
    -1 reads_1.fq.gz \
    -2 reads_2.fq.gz \
    --summary-file sample.hisat2.summary.log \
    --threads $task.cpus \
    $seq_center \
    $unaligned \
    --no-mixed \
    --no-discordant \
    $args \
| samtools view -bS -f 2 -F 4 -F 8 -F 256 - > sample.bam
```

### HISAT2 default and extra args

The default HISAT2 arguments are:

```bash
--met-stderr --new-summary --dta
```

Meaning:

- `--met-stderr`: prints HISAT2 internal metrics to `stderr`;
- `--new-summary`: reports the alignment summary in the newer, machine-friendly format;
- `--dta`: reports alignments tailored for downstream transcriptome assembly, for example StringTie.

Additional HISAT2 arguments can be passed through:

```groovy
params.extra_hisat2_align_args
```

The final HISAT2 arguments are built as:

```groovy
ext.args = [
    '--met-stderr --new-summary --dta',
    params.extra_hisat2_align_args ?: ''
].join(' ').trim()
```

Example:

```bash
--extra_hisat2_align_args "-k 10"
```

This produces:

```bash
--met-stderr --new-summary --dta -k 10
```

---

## SAM flag reference for STAR and HISAT2 BAM filtering

| Flag | Hex | Meaning | Action in pipeline |
|------|-----|---------|--------------------|
| `2` | `0x0002` | Proper pair | Kept in paired-end mode with `-f 2` |
| `4` | `0x0004` | Read unmapped | Removed in single-end and paired-end |
| `8` | `0x0008` | Mate unmapped | Removed in paired-end only |
| `256` | `0x0100` | Secondary alignment | Removed in single-end and paired-end |

---

