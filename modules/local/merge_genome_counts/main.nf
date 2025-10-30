process MERGE_GENOME_COUNTS {
    tag "merge_genome_counts"
    label 'process_medium'

    container 'docker://pdichiaro/env_r_ngs'

    input:
    path count_files
    val feature_types
    path annotation_matrix

    output:
    path "*_counts_merged.txt"    , emit: merged_counts
    path "*_merge_summary.txt"    , emit: summary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "genome"
    def features = feature_types ?: "transcript,intron,exon,5utr,3utr"
    """
    # Create input directory and copy files
    mkdir -p input_counts
    cp ${count_files} input_counts/

    merge_genome_counts.r \\
        -i input_counts \\
        -p "_combined_counts.txt" \\
        -o ${prefix} \\
        -t ${features} \\
        -a ${annotation_matrix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
