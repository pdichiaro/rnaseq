process ANNOTATION_MATRIX {
    tag "$txdb_sqlite"
    label 'process_medium'

    container 'docker://pdichiaro/env_r_ngs'

    input:
    path txdb_sqlite
    path metadata

    output:
    path "*_annotation_matrix.txt"    , emit: annotation_matrix
    path "*_genomic_features.rds"     , emit: genomic_features
    path "*_summary_stats.txt"        , emit: summary_stats
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${txdb_sqlite.baseName}"
    """
    annotation_matrix.r \\
        -q ${txdb_sqlite} \\
        -m ${metadata} \\
        -o ${prefix} \\
        -c $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-GenomicFeatures: \$(Rscript -e "library(GenomicFeatures); cat(as.character(packageVersion('GenomicFeatures')))")
        bioconductor-GenomicAlignments: \$(Rscript -e "library(GenomicAlignments); cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-GenomicRanges: \$(Rscript -e "library(GenomicRanges); cat(as.character(packageVersion('GenomicRanges')))")
    END_VERSIONS
    """
}
