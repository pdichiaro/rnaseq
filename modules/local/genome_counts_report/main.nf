process GENOME_COUNTS_REPORT {
    tag "genome_counts"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path summary_files

    output:
    path "genome_read_counts.txt", emit: report
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat <<'EOF' > create_report.py
#!/usr/bin/env python3
import os
import re

# Read all summary files
summary_files = [f for f in os.listdir('.') if f.endswith('_summary.txt')]

# Create report
with open('genome_read_counts.txt', 'w') as out:
    out.write("=" * 80 + "\\n")
    out.write("GENOME READ COUNTS REPORT\\n")
    out.write("=" * 80 + "\\n\\n")
    out.write("This report shows the distribution of reads across genomic features\\n")
    out.write("for each sample in the analysis.\\n\\n")
    
    for summary_file in sorted(summary_files):
        sample_id = summary_file.replace('_summary.txt', '')
        out.write("-" * 80 + "\\n")
        out.write(f"Sample: {sample_id}\\n")
        out.write("-" * 80 + "\\n")
        
        with open(summary_file, 'r') as f:
            content = f.read()
            out.write(content)
        
        out.write("\\n")
    
    out.write("=" * 80 + "\\n")
    out.write("END OF REPORT\\n")
    out.write("=" * 80 + "\\n")

print(f"Processed {len(summary_files)} samples")
EOF

    python3 create_report.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
