process DESEQ2_TRANSFORM {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(run_dir)
    path pca_header
    path clustering_header

    output:
    tuple val(meta), path("*_mqc.tsv"), optional: true, emit: multiqc_files
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Find and transform DESeq2 PCA files
    find ${run_dir} -name "*.pca.vals.txt" -o -name "*pca.vals.txt" | while read pca_file; do
        if [ -f "\$pca_file" ]; then
            base_name=\$(basename "\$pca_file" .txt)
            cat ${pca_header} "\$pca_file" > "\${base_name}_mqc.tsv"
            echo "Created \${base_name}_mqc.tsv"
        fi
    done

    # Find and transform DESeq2 sample distance files
    find ${run_dir} -name "*.sample.dists.txt" -o -name "*sample.dists.txt" | while read dist_file; do
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
