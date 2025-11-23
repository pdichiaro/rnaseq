process DEEPTOOLS_BIGWIG_NORM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::deeptools=3.5.1 bioconda::samtools=1.16.1"
    container 'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0'

    input:
    tuple val(meta), path(bam), path(bai), val(scaling)

    output:
    tuple val(meta), path("*.unstranded.norm.bw"), emit: unstranded_bw
    tuple val(meta), path("*.fwd.norm.bw")       , optional:true, emit: fw_bw
    tuple val(meta), path("*.rev.norm.bw")       , optional:true, emit: rev_bw
    path "versions.yml"                          , emit: versions
 
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pe = meta.single_end ? 'single' : 'paired'
    
    def strandedness = 'unstranded'
    if (meta.strandedness == 'forward') {
        strandedness = 'forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'reverse'
    }

    if(strandedness == 'unstranded'){
        """
        echo "Sample: ${meta.id}"
        echo "Library strandedness: unstranded"
        echo "Data type: $pe"
        echo "Scaling factor: $scaling"
        echo "Prefix: $prefix"
        
        bamCoverage \\
                --numberOfProcessors $task.cpus \\
                --binSize 1 \\
                --scaleFactor $scaling \\
                --bam $bam \\
                -o ${prefix}.unstranded.norm.bw

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
        END_VERSIONS
        """
    } else {
        if( pe == 'single' ){
            if(strandedness == 'forward'){
                """
                echo "Sample: ${meta.id}"
                echo "Library strandedness: forward"
                echo "Data type: single-end"
                echo "Scaling factor: $scaling"
                echo "Prefix: $prefix"
                
                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    -o ${prefix}.unstranded.norm.bw

                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    --samFlagExclude 16 \\
                    -o ${prefix}.fwd.norm.bw

                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    --samFlagInclude 16 \\
                    -o ${prefix}.rev.norm.bw

                cat <<-END_VERSIONS > versions.yml
                "${task.process}":
                    deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
                END_VERSIONS
                """
            } else {
                """
                # Read scaling factor from individual per-sample file
                
                echo "Sample: ${meta.id}"
                echo "Library strandedness: reverse"
                echo "Data type: single-end"
                echo "Scaling factor: $scaling"
                echo "Prefix: $prefix"
                
                # 1. Unstranded BigWig (always generate)
                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    -o ${prefix}.unstranded.norm.bw

                # 2. Forward strand BigWig (reverse library: FLIPPED logic)
                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    --samFlagInclude 16 \\
                    -o ${prefix}.fwd.norm.bw

                # 3. Reverse strand BigWig (reverse library: FLIPPED logic)
                bamCoverage \\
                    --bam $bam \\
                    --numberOfProcessors $task.cpus \\
                    --binSize 1 \\
                    --scaleFactor $scaling \\
                    --samFlagExclude 16 \\
                    -o ${prefix}.rev.norm.bw

                cat <<-END_VERSIONS > versions.yml
                "${task.process}":
                    deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
                END_VERSIONS
                """
            }
        }else{
            if( pe == 'paired' ){
                if(strandedness == 'forward'){
                    """
                    # Read scaling factor from individual per-sample file
                    
                    echo "Sample: ${meta.id}"
                    echo "Library strandedness: forward"
                    echo "Data type: paired-end"
                    echo "Scaling factor: $scaling"
                    echo "Prefix: $prefix"

                    # 1. Unstranded BigWig (always generate)
                    bamCoverage \\
                        --bam $bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.unstranded.norm.bw

                    # 2. Forward strand BigWig (forward library: normal logic)
                    
                    # include reads that are 2nd in a pair (128); exclude reads mapped to reverse strand (16)
                    samtools view -b -f 128 -F 16 $bam > ${prefix}.fwd1.bam
                    # exclude reads mapped to reverse strand (16) and first in a pair (64): 64 + 16 = 80
                    samtools view -b -f 80 $bam > ${prefix}.fwd2.bam
                    # combine the temporary files
                    samtools merge -f ${prefix}.fwd.bam ${prefix}.fwd1.bam ${prefix}.fwd2.bam
                    rm -rf ${prefix}.fwd1.bam ${prefix}.fwd2.bam
                    # index the filtered BAM file
                    samtools index ${prefix}.fwd.bam
                    # run bamCoverage
                    bamCoverage \\
                        --bam ${prefix}.fwd.bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.fwd.norm.bw
                    rm -rf ${prefix}.fwd.bam ${prefix}.fwd.bam.bai

                    # 3. Reverse strand BigWig (forward library: normal logic)
                    
                    # include reads that map to the reverse strand (16) and are second in a pair (128): 128 + 16 = 144
                    samtools view -b -f 144 $bam > ${prefix}.rev1.bam
                    # include reads that are first in a pair (64), but exclude those mapped to reverse strand (16)
                    samtools view -b -f 64 -F 16 $bam > ${prefix}.rev2.bam
                    # merge the temporary files
                    samtools merge -f ${prefix}.rev.bam ${prefix}.rev1.bam ${prefix}.rev2.bam
                    rm -rf ${prefix}.rev1.bam ${prefix}.rev2.bam
                    # index the merged, filtered BAM file
                    samtools index ${prefix}.rev.bam
                    # run bamCoverage
                    bamCoverage \\
                        --bam ${prefix}.rev.bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.rev.norm.bw
                    rm -rf ${prefix}.rev.bam ${prefix}.rev.bam.bai

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
                    END_VERSIONS
                    """
                } else {
                    """
                    # Read scaling factor from individual per-sample file
                    
                    echo "Sample: ${meta.id}"
                    echo "Library strandedness: reverse"
                    echo "Data type: paired-end"
                    echo "Scaling factor: $scaling"
                    echo "Prefix: $prefix"

                    # 1. Unstranded BigWig (always generate)
                    bamCoverage \\
                        --bam $bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.unstranded.norm.bw

                    # 2. Forward strand BigWig (reverse library: FLIPPED logic)
                    
                    # For reverse library: reads mapped as 'reverse' actually represent forward transcripts
                    # include reads that map to the reverse strand (16) and are second in a pair (128): 128 + 16 = 144
                    samtools view -b -f 144 $bam > ${prefix}.fwd1.bam
                    # include reads that are first in a pair (64), but exclude those mapped to reverse strand (16)
                    samtools view -b -f 64 -F 16 $bam > ${prefix}.fwd2.bam
                    # merge the temporary files
                    samtools merge -f ${prefix}.fwd.bam ${prefix}.fwd1.bam ${prefix}.fwd2.bam
                    rm -rf ${prefix}.fwd1.bam ${prefix}.fwd2.bam
                    # index the merged, filtered BAM file
                    samtools index ${prefix}.fwd.bam
                    # run bamCoverage
                    bamCoverage \\
                        --bam ${prefix}.fwd.bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.fwd.norm.bw
                    rm -rf ${prefix}.fwd.bam ${prefix}.fwd.bam.bai

                    # 3. Reverse strand BigWig (reverse library: FLIPPED logic)
                    
                    # For reverse library: reads mapped as 'forward' actually represent reverse transcripts
                    # include reads that are 2nd in a pair (128); exclude reads mapped to reverse strand (16)
                    samtools view -b -f 128 -F 16 $bam > ${prefix}.rev1.bam
                    # exclude reads mapped to reverse strand (16) and first in a pair (64): 64 + 16 = 80
                    samtools view -b -f 80 $bam > ${prefix}.rev2.bam
                    # combine the temporary files
                    samtools merge -f ${prefix}.rev.bam ${prefix}.rev1.bam ${prefix}.rev2.bam
                    rm -rf ${prefix}.rev1.bam ${prefix}.rev2.bam
                    # index the filtered BAM file
                    samtools index ${prefix}.rev.bam
                    # run bamCoverage
                    bamCoverage \\
                        --bam ${prefix}.rev.bam \\
                        --numberOfProcessors $task.cpus \\
                        --binSize 1 \\
                        --scaleFactor $scaling \\
                        -o ${prefix}.rev.norm.bw
                    rm -rf ${prefix}.rev.bam ${prefix}.rev.bam.bai

                    cat <<-END_VERSIONS > versions.yml
                    "${task.process}":
                        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
                    END_VERSIONS
                    """
                }
            }
        }
    }
}

