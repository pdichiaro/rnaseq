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
    path "*_mqc.txt", optional: true, emit: multiqc_files
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_name = deseq2_file.getName()
    // Create output name with _mqc.txt suffix (replace .txt with _mqc.txt)
    def output_name = file_name.replaceAll(/\.txt$/, '_mqc.txt')
    """
    # Detect quantifier and level from filename for unique IDs and section anchors
    # Filename patterns: kallisto.deseq2.all_genes.*, kallisto.deseq2.invariant_genes.*, etc.
    QUANTIFIER=""
    QUANTIFIER_SHORT=""
    LEVEL=""
    SECTION_NAME=""
    
    # Extract quantifier
    if [[ "${file_name}" == kallisto.* ]]; then
        QUANTIFIER="deseq2-kallisto-qc"
        QUANTIFIER_SHORT="kallisto"
        PARENT_NAME="DESeq2 Kallisto QC"
    elif [[ "${file_name}" == salmon.deseq2.* ]]; then
        QUANTIFIER="deseq2-salmon-qc"
        QUANTIFIER_SHORT="salmon"
        PARENT_NAME="DESeq2 Salmon QC"
    elif [[ "${file_name}" == star.rsem.* ]]; then
        QUANTIFIER="deseq2-star-rsem-qc"
        QUANTIFIER_SHORT="star_rsem"
        PARENT_NAME="DESeq2 STAR RSEM QC"
    elif [[ "${file_name}" == star.salmon.* ]]; then
        QUANTIFIER="deseq2-star-salmon-qc"
        QUANTIFIER_SHORT="star_salmon"
        PARENT_NAME="DESeq2 STAR Salmon QC"
    elif [[ "${file_name}" == star.genome.* ]]; then
        QUANTIFIER="deseq2-star-genome-qc"
        QUANTIFIER_SHORT="star_genome"
        PARENT_NAME="DESeq2 STAR Genome QC"
    elif [[ "${file_name}" == hisat2.genome.* ]]; then
        QUANTIFIER="deseq2-hisat2-genome-qc"
        QUANTIFIER_SHORT="hisat2"
        PARENT_NAME="DESeq2 HISAT2 QC"
    else
        echo "Warning: Could not determine quantifier from filename: ${file_name}"
        QUANTIFIER="deseq2-qc"
        QUANTIFIER_SHORT="unknown"
        PARENT_NAME="DESeq2 QC"
    fi
    
    # Extract level (all_genes or invariant_genes)
    if [[ "${file_name}" == *.all_genes.* ]]; then
        LEVEL="all_genes"
        SECTION_NAME="All Genes"
    elif [[ "${file_name}" == *.invariant_genes.* ]]; then
        LEVEL="invariant_genes"
        SECTION_NAME="Invariant Genes"
    else
        LEVEL="unknown"
        SECTION_NAME="Unknown Level"
    fi
    
    echo "Detected quantifier: \${QUANTIFIER_SHORT}, level: \${LEVEL}, section: \${QUANTIFIER}"
    echo "Output file will be: ${output_name}"

    # Determine number prefix based on gene set (01-04 for all_genes, 05-08 for invariant_genes)
    OFFSET=0
    if [[ "\${LEVEL}" == "invariant_genes" ]]; then
        OFFSET=4
    fi
    
    # Add appropriate header to each file type for MultiQC custom content module
    # Each plot gets nested under parent section with unique ID
    # Numbers added to SECTION_NAME for visible ordering (MultiQC sorts by section_name in nested sections)
    # Number ranges: 01-04 for All Genes, 05-08 for Invariant Genes
    # Plot titles remain clean without numbers
    # IMPORTANT: Check .pca.top*.vals.txt BEFORE .pca.vals.txt to avoid false matches
    if [[ "${file_name}" == *".pca.top"*".vals.txt" ]]; then
        # PCA top variable genes (pattern: *.pca.top500.vals.txt) - ORDER: 4 or 8
        PLOT_NUM=\$((4 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_pca_top500_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="\$(printf '%02d' \$PLOT_NUM). PCA Top 500 (\${SECTION_NAME})"
        PLOT_TITLE="PCA Top 500 (\${SECTION_NAME})"
        numbered_output="\$(printf '%02d' \$PLOT_NUM)_${output_name}"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${pca_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${numbered_output}"
        echo "Created \${numbered_output} with PCA-500 header (ID: \${PLOT_ID}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".pca.vals.txt" ]]; then
        # PCA all genes (pattern: *.pca.vals.txt) - ORDER: 3 or 7
        PLOT_NUM=\$((3 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_pca_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="\$(printf '%02d' \$PLOT_NUM). PCA (\${SECTION_NAME})"
        PLOT_TITLE="PCA (\${SECTION_NAME})"
        numbered_output="\$(printf '%02d' \$PLOT_NUM)_${output_name}"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${pca_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${numbered_output}"
        echo "Created \${numbered_output} with PCA header (ID: \${PLOT_ID}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".sample.dists."* ]]; then
        # Sample distance - ORDER: 2 or 6
        PLOT_NUM=\$((2 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_sample_distance_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="\$(printf '%02d' \$PLOT_NUM). Sample Distances (\${SECTION_NAME})"
        PLOT_TITLE="Sample Distances (\${SECTION_NAME})"
        numbered_output="\$(printf '%02d' \$PLOT_NUM)_${output_name}"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${clustering_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${numbered_output}"
        echo "Created \${numbered_output} with sample distance header (ID: \${PLOT_ID}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".read.distribution.normalized."* ]]; then
        # Read distribution - ORDER: 1 or 5
        PLOT_NUM=\$((1 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_read_distribution_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="\$(printf '%02d' \$PLOT_NUM). Read Distribution (\${SECTION_NAME})"
        PLOT_TITLE="Read Distribution (\${SECTION_NAME})"
        numbered_output="\$(printf '%02d' \$PLOT_NUM)_${output_name}"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${read_dist_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${numbered_output}"
        echo "Created \${numbered_output} with read distribution header (ID: \${PLOT_ID}, parent: \${QUANTIFIER})"
    else
        # Unknown file type - copy as-is with _mqc.txt suffix
        cp ${deseq2_file} "${output_name}"
        echo "Copied ${output_name} without header (unknown type)"
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    bash: \$(bash --version | head -n1 | awk '{print \$4}')
END_VERSIONS
    """

    stub:
    """
    touch stub_file.txt

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    bash: \$(bash --version | head -n1 | awk '{print \$4}')
END_VERSIONS
    """
}
