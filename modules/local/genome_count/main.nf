process GENOME_COUNT {
    tag "$meta.id"
    label 'process_high'

    container 'docker://pdichiaro/env_r_ngs'

    input:
    tuple val(meta), path(bam), path(bai)
    path genomic_features

    output:
    tuple val(meta), path("*_combined_counts.txt")    , emit: combined_counts
    tuple val(meta), path("*_transcript_counts.txt")  , emit: transcript_counts
    tuple val(meta), path("*_intron_counts.txt")      , emit: intron_counts
    tuple val(meta), path("*_exon_counts.txt")        , emit: exon_counts
    tuple val(meta), path("*_5utr_counts.txt")        , emit: utr5_counts
    tuple val(meta), path("*_3utr_counts.txt")        , emit: utr3_counts
    tuple val(meta), path("*_counts.rds")             , emit: counts_rds
    tuple val(meta), path("*_summary.txt")            , emit: summary
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = meta.strandedness ?: 'non_specific'
    def single_end = meta.single_end ? 'single' : 'paired'
    
    // Convert strandedness to script format
    def strand_param = 'non_specific'
    if (strandedness == 'forward') {
        strand_param = 'forward'
    } else if (strandedness == 'reverse') {
        strand_param = 'reverse'
    }
    
    """
    genome_count.r \\
        -f ${genomic_features} \\
        -b ${bam} \\
        -i ${meta.id} \\
        -s ${strand_param} \\
        -p ${single_end} \\
        -o ${prefix} \\
        -c $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-GenomicFeatures: \$(Rscript -e "library(GenomicFeatures); cat(as.character(packageVersion('GenomicFeatures')))")
        bioconductor-GenomicAlignments: \$(Rscript -e "library(GenomicAlignments); cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-GenomicRanges: \$(Rscript -e "library(GenomicRanges); cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-Rsamtools: \$(Rscript -e "library(Rsamtools); cat(as.character(packageVersion('Rsamtools')))")
    END_VERSIONS
    """
}
