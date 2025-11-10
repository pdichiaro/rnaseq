process KALLISTO_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(run_info_json)

    output:
    path "*_mqc.tsv"      , emit: multiqc_files
    path "*_mqc.txt"      , emit: multiqc_txt_files
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import json
    import sys
    
    # Read the run_info.json file
    with open('${run_info_json}', 'r') as f:
        data = json.load(f)
    
    # Extract the values we want to plot
    sample_id = '${meta.id}'
    n_processed = data.get('n_processed', 0)
    n_pseudoaligned = data.get('n_pseudoaligned', 0)
    n_unique = data.get('n_unique', 0)
    
    # Calculate percentages
    pct_pseudoaligned = (n_pseudoaligned / n_processed * 100) if n_processed > 0 else 0
    pct_unique = (n_unique / n_processed * 100) if n_processed > 0 else 0
    
    # Create counts file for bargraph (TSV format)
    counts_file = '${prefix}_kallisto_alignment_stats_cnt_mqc.tsv'
    with open(counts_file, 'w') as out:
        # Write column headers only (for merge to detect)
        out.write("Sample\\tTotal Processed\\tPseudoaligned\\tUnique\\n")
        # Write the data with counts
        out.write(f"{sample_id}\\t{n_processed}\\t{n_pseudoaligned}\\t{n_unique}\\n")
    
    # Create percentages file for bargraph (TXT format)
    pct_file = '${prefix}_kallisto_alignment_stats_pct_mqc.txt'
    with open(pct_file, 'w') as out:
        # Write column headers only (for merge to detect)
        out.write("Sample\\tPseudoaligned (%)\\tUnique (%)\\n")
        # Write the data with percentages
        out.write(f"{sample_id}\\t{pct_pseudoaligned:.2f}\\t{pct_unique:.2f}\\n")
    
    # Debug: print what was created
    print(f"Created {counts_file} and {pct_file} for sample {sample_id}")
    print(f"  Total Processed: {n_processed}")
    print(f"  Pseudoaligned: {n_pseudoaligned} ({pct_pseudoaligned:.2f}%)")
    print(f"  Unique: {n_unique} ({pct_unique:.2f}%)")
    
    # Create versions file
    with open('versions.yml', 'w') as v:
        v.write('"${task.process}":\\n')
        v.write(f'    python: "{sys.version.split()[0]}"\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_kallisto_alignment_stats_cnt_mqc.tsv
    touch ${prefix}_kallisto_alignment_stats_pct_mqc.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
