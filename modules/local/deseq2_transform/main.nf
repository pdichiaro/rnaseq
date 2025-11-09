process DESEQ2_TRANSFORM {
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path deseq2_files
    path pca_header
    path clustering_header

    output:
    path "*_mqc.tsv", optional: true, emit: multiqc_files
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Transform DESeq2 PCA files
    for pca_file in *.pca.vals.txt *pca.vals.txt; do
        if [ -f "\$pca_file" ]; then
            base_name=\$(basename "\$pca_file" .txt)
            cat ${pca_header} "\$pca_file" > "\${base_name}_mqc.tsv"
            echo "Created \${base_name}_mqc.tsv"
        fi
    done

    # Transform DESeq2 sample distance files
    for dist_file in *.sample.dists.txt *sample.dists.txt; do
        if [ -f "\$dist_file" ]; then
            base_name=\$(basename "\$dist_file" .txt)
            cat ${clustering_header} "\$dist_file" > "\${base_name}_mqc.tsv"
            echo "Created \${base_name}_mqc.tsv"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch stub_pca.vals_mqc.tsv
    touch stub_sample.dists_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
