process MULTIQC_GENOME_COUNTS {
    tag "genome_counts_mqc"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path summary_files

    output:
    path "*_mqc.tsv", emit: mqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "genome_counts"
    """
    genome_counts_multiqc.py \\
        ${summary_files} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "genome_counts"
    """
    touch ${prefix}_generalstats_mqc.tsv
    touch ${prefix}_distribution_mqc.tsv
    touch ${prefix}_genes_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
