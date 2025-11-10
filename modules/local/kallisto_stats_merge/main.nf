process KALLISTO_STATS_MERGE {
    label 'process_single'
    
    publishDir "${params.outdir}/kallisto", mode: params.publish_dir_mode

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path stats_cnt_files, stageAs: 'input_cnt/*'
    path stats_pct_files, stageAs: 'input_pct/*'

    output:
    path "kallisto_alignment_stats_cnt.txt", emit: cnt_file
    path "kallisto_alignment_stats_pct.txt", emit: pct_file
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/bin/bash
    
    # Create merged counts file - just data, no comment headers
    echo "Sample\tTotal Processed\tPseudoaligned\tUnique" > kallisto_alignment_stats_cnt.txt
    
    for file in input_cnt/*_kallisto_alignment_stats_cnt_mqc.tsv; do
        if [ -f "\$file" ]; then
            tail -n +2 "\$file" >> kallisto_alignment_stats_cnt.txt
        fi
    done
    
    # Create merged percentages file - just data, no comment headers
    echo "Sample\tPseudoaligned (%)\tUnique (%)" > kallisto_alignment_stats_pct.txt
    
    for file in input_pct/*_kallisto_alignment_stats_pct_mqc.txt; do
        if [ -f "\$file" ]; then
            tail -n +2 "\$file" >> kallisto_alignment_stats_pct.txt
        fi
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch kallisto_alignment_stats_cnt.txt
    touch kallisto_alignment_stats_pct.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """
}
