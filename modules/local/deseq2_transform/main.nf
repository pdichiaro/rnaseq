process DESEQ2_TRANSFORM {
    label 'process_single'
    tag "$deseq2_file"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path deseq2_file
    path pca_header
    path clustering_header
    path read_dist_header

    output:
    path "*_mqc.tsv", optional: true, emit: multiqc_files
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_name = deseq2_file.getName()
    def base_name = file_name.replaceAll(/\.txt$/, '')
    """
    # Determine file type and apply appropriate header
    if [[ "${file_name}" == *".pca."* ]]; then
        # PCA file - use PCA header
        cat ${pca_header} ${deseq2_file} > "${base_name}_mqc.tsv"
        echo "Created ${base_name}_mqc.tsv with PCA header"
    elif [[ "${file_name}" == *".sample.dists."* ]]; then
        # Sample distance file - use clustering header
        cat ${clustering_header} ${deseq2_file} > "${base_name}_mqc.tsv"
        echo "Created ${base_name}_mqc.tsv with clustering header"
    elif [[ "${file_name}" == *".read.distribution.normalized."* ]]; then
        # Read distribution file - use read distribution header
        cat ${read_dist_header} ${deseq2_file} > "${base_name}_mqc.tsv"
        echo "Created ${base_name}_mqc.tsv with read distribution header"
    else
        echo "Warning: Unknown file type: ${file_name}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch stub_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
