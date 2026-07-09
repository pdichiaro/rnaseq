# RNA-seq Alignment Filtering Summary

## STAR and HISAT2 filtering workflow

This document summarizes how aligned reads are retained or excluded after genome alignment in the RNA-seq pipeline when using **STAR** or **HISAT2**.

The key point is that the two aligners are handled differently:

- **STAR**: filtering is mainly controlled inside STAR through `--outFilter*` alignment parameters. No additional `samtools view` MAPQ/flag filter is applied immediately after STAR alignment.
- **HISAT2**: HISAT2 output is piped directly into `samtools view`, where unmapped reads, mate-unmapped reads in paired-end mode, and secondary alignments are removed before BAM sorting.

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
│  Standard mode additionally uses:                                 │
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
    ├─ Multimapping is controlled by --outFilterMultimapNmax 20
    ├─ Mismatch filtering is stricter in RSEM mode
    └─ No explicit post-STAR samtools MAPQ filter is applied here
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 2: BAM_SORT_STATS_SAMTOOLS                                  │
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
    └─ Reads are not additionally filtered by MAPQ or SAM flags here
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 3: OPTIONAL DEDUPLICATION / DUPLICATE MARKING               │
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
| `--outFilterMultimapNmax 20` | Maximum number of loci a read can map to | Reads mapping to more than 20 loci are excluded by STAR |
| `--outFilterType BySJout` | Splice-junction-aware filtering | Used in RSEM-compatible mode |
| `--outFilterMismatchNmax 999` | Absolute mismatch limit | Effectively permissive absolute mismatch threshold |
| `--outFilterMismatchNoverLmax 0.04` | Mismatch fraction threshold | In RSEM-compatible mode, reads with mismatch fraction above 4% are excluded |
| `--outSAMtype BAM Unsorted` | BAM output | STAR writes an unsorted BAM directly |
| `--outSAMattributes NH HI AS NM MD` | SAM tags | Keeps multimapping and alignment score/mismatch tags useful for downstream QC |

### What STAR does **not** do in this pipeline

The STAR branch does **not** apply an additional post-alignment command such as:

```bash
samtools view -q <MAPQ> -F <FLAGS>
```

Therefore, in the default STAR workflow there is no explicit post-STAR:

- MAPQ threshold
- removal of secondary alignments by `samtools view -F 256`
- removal of supplementary alignments by `samtools view -F 2048`
- proper-pair filtering by `samtools view -f 2`
- insert-size / fragment-length filtering
- blacklist-region filtering

Filtering is therefore primarily determined by STAR alignment parameters and by any user-provided `params.extra_star_align_args`.

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
│    samtools view -bS -F 4 -F 8 -F 256 - > sample.bam              │
│                                                                   │
│  Removes:                                                         │
│    -F 4     unmapped reads                                        │
│    -F 8     reads with unmapped mate                              │
│    -F 256   secondary alignments                                  │
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
| samtools view -bS -F 4 -F 8 -F 256 - > sample.bam
```

---

## SAM flag reference for HISAT2 filtering

| Flag | Hex | Meaning | Action in HISAT2 branch |
|------|-----|---------|--------------------------|
| `4` | `0x0004` | Read unmapped | Removed in single-end and paired-end |
| `8` | `0x0008` | Mate unmapped | Removed in paired-end only |
| `256` | `0x0100` | Secondary alignment | Removed in single-end and paired-end |

---
