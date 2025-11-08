process TXIMPORT_MERGE_COUNTS {
    tag "tximport_merge_counts"
    label "process_medium"

    container 'docker://pdichiaro/env_r_ngs'

    input:
    tuple val(meta), path(gene_counts)
    tuple val(meta2), path(gene_tpm)
    tuple val(meta3), path(gene_lengths)
    tuple val(meta4), path(transcript_counts)
    tuple val(meta5), path(transcript_tpm)
    tuple val(meta6), path(transcript_lengths)
    path annotation_matrix

    output:
    path "tximport.merged.gene_counts.tsv"              , emit: counts_gene
    path "tximport.merged.gene_tpm.tsv"                 , emit: tpm_gene
    path "tximport.merged.gene_lengths.tsv"             , emit: lengths_gene
    path "tximport.merged.transcript_counts.tsv"        , emit: counts_transcript
    path "tximport.merged.transcript_tpm.tsv"           , emit: tpm_transcript
    path "tximport.merged.transcript_lengths.tsv"       , emit: lengths_transcript
    path "tximport_merge_summary.txt"                   , emit: summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "tximport"
    """
    tximport_merge_counts.r \\
        --gene_counts ${gene_counts} \\
        --gene_tpm ${gene_tpm} \\
        --gene_lengths ${gene_lengths} \\
        --transcript_counts ${transcript_counts} \\
        --transcript_tpm ${transcript_tpm} \\
        --transcript_lengths ${transcript_lengths} \\
        --annotation_matrix ${annotation_matrix} \\
        --output_prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "tximport"
    """
    touch ${prefix}.merged.gene_counts.tsv
    touch ${prefix}.merged.gene_tpm.tsv
    touch ${prefix}.merged.gene_lengths.tsv
    touch ${prefix}.merged.transcript_counts.tsv
    touch ${prefix}.merged.transcript_tpm.tsv
    touch ${prefix}.merged.transcript_lengths.tsv
    touch ${prefix}_merge_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}