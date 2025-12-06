process GENOME_SECTION_HEADER {
    tag "$quantifier"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    val quantifier

    output:
    path "*_section_header_mqc.txt", emit: section_header
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def section_title = quantifier.replaceAll('_', ' ').split(' ').collect { it.capitalize() }.join(' ')
    """
    cat <<-END_SECTION > ${quantifier}_genome_section_header_mqc.txt
    # section_id: 'genome_${quantifier}'
    # section_name: '${section_title} - Genome'
    # section_anchor: '${quantifier}_genome'
    # description: 'Genome-level count statistics for ${section_title}'
    END_SECTION

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${quantifier}_genome_section_header_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
