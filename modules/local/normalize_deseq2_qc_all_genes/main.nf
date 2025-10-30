process NORMALIZE_DESEQ2_QC_ALL_GENES {
    tag "deseq2_qc_all_genes"
    label "process_medium"

    container 'docker://pdichiaro/env_r_ngs'

    input:
    path counts
    path pca_header_multiqc
    path clustering_header_multiqc

    output:
    path "scaling_dat.txt"      , emit: scaling_factors
    path "normalized_counts.txt", emit: normalized_counts
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()
    prefix = task.ext.prefix ?: "deseq2_all_genes"
    """
    deseq2_qc_all_genes.r \\
        --count_file $counts \\
        --outdir ./ \\
        --outprefix $prefix \\
        --cores $task.cpus \\
        $args

    if [ -f "R_sessionInfo.log" ]; then
        # Handle PCA files
        sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" <$pca_header_multiqc > pca_header.tmp
        sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" pca_header.tmp
        cat pca_header.tmp *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv
        rm pca_header.tmp

        # Handle clustering files
        sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" <$clustering_header_multiqc > clustering_header.tmp
        sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" clustering_header.tmp
        cat clustering_header.tmp *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
        rm clustering_header.tmp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    prefix = task.ext.prefix ?: "deseq2_all_genes"
    """
    touch ${label_lower}.pca.vals_mqc.tsv
    touch ${label_lower}.sample.dists_mqc.tsv
    touch ${prefix}.dds.RData
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.plots.pdf
    touch ${prefix}.sample.dists.txt
    touch normalized_counts.txt
    touch scaling_dat.txt
    touch R_sessionInfo.log

    mkdir size_factors
    touch size_factors/${prefix}.size_factors.RData

    cat <<-END_VERSIONS >versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}