#!/usr/bin/env nextflow

// Test to understand channel behavior with path-only outputs

process MAKE_FILES_1 {
    output:
    path "*.txt", emit: files
    
    script:
    """
    touch file1_process1.txt
    touch file2_process1.txt
    """
}

process MAKE_FILES_2 {
    output:
    path "*.txt", emit: files
    
    script:
    """
    touch file1_process2.txt
    touch file2_process2.txt
    """
}

process MAKE_FILES_3 {
    output:
    path "*.txt", emit: files
    
    script:
    """
    touch file1_process3.txt
    """
}

process COLLECT_ALL {
    publishDir "results", mode: 'copy'
    
    input:
    path files
    
    output:
    path "collected_files.txt"
    
    script:
    """
    echo "Files received:" > collected_files.txt
    for f in ${files}; do
        echo "  - \$f" >> collected_files.txt
    done
    """
}

workflow {
    // Test 1: Mix without flatten
    ch_files = Channel.empty()
    
    MAKE_FILES_1()
    MAKE_FILES_2()
    MAKE_FILES_3()
    
    ch_files = ch_files.mix(MAKE_FILES_1.out.files)
    ch_files = ch_files.mix(MAKE_FILES_2.out.files)
    ch_files = ch_files.mix(MAKE_FILES_3.out.files)
    
    // Try without flatten
    ch_files.view { "Channel item (no flatten): $it (class: ${it.class.name})" }
    
    COLLECT_ALL(ch_files.collect())
    
    // Try with flatten
    ch_files.flatten().view { "Channel item (with flatten): $it (class: ${it.class.name})" }
}
