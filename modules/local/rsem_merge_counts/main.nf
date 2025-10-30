process RSEM_MERGE_COUNTS {
    label "process_medium"

    container 'docker://pdichiaro/env_r_ngs'

    input:
    path ('genes/*')
    path ('isoforms/*')
    path annotation_matrix

    output:
    path "rsem.merged.gene_counts.tsv"      , emit: counts_gene
    path "rsem.merged.gene_tpm.tsv"         , emit: tpm_gene
    path "rsem.merged.transcript_counts.tsv", emit: counts_transcript
    path "rsem.merged.transcript_tpm.tsv"   , emit: tpm_transcript
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "rsem"
    """
    rsem_merge_counts.r \\
        -g genes \\
        -t isoforms \\
        -o ${prefix} \\
        -a ${annotation_matrix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """

    stub:
    """
    touch rsem.merged.gene_counts.tsv
    touch rsem.merged.gene_tpm.tsv
    touch rsem.merged.transcript_counts.tsv
    touch rsem.merged.transcript_tpm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
