process DESEQ2_SECTION_HEADER {
    label 'process_single'
    tag "$quantifier"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    val quantifier  // e.g., 'kallisto', 'salmon', 'star_rsem', etc.

    output:
    path "*_section_header_mqc.txt", emit: section_header
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Convert quantifier to section anchor format
    def section_anchor = "deseq2-${quantifier.replaceAll('_', '-')}-qc"
    def section_title = quantifier.replaceAll('_', ' ').split().collect { it.capitalize() }.join(' ')
    """
    cat > ${quantifier}_deseq2_section_header_mqc.txt <<HEADER_EOF
#id: '${section_anchor}'
#section_name: 'DESeq2 ${section_title} QC'
#section_anchor: '${section_anchor}'
#description: 'Quality control plots from DESeq2 normalization using ${section_title} quantification data.'
#plot_type: 'html'
#pconfig:
#    namespace: 'DESeq2 QC'
HEADER_EOF

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
	END_VERSIONS
    """

    stub:
    """
    touch ${quantifier}_deseq2_section_header_mqc.txt

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
	END_VERSIONS
    """
}
