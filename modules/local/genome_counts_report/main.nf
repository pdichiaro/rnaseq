process GENOME_COUNTS_REPORT {
    tag "genome_counts"
    label 'process_single'
    publishDir "${params.outdir}/${params.aligner}/genome", mode: params.publish_dir_mode

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path count_files

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

# Find all count files
all_files = os.listdir('.')

# Group files by sample
samples = {}
for f in all_files:
    if '_transcript_counts.txt' in f:
        sample = f.replace('_transcript_counts.txt', '')
        if sample not in samples:
            samples[sample] = {}
        samples[sample]['transcript'] = f
    elif '_intron_counts.txt' in f:
        sample = f.replace('_intron_counts.txt', '')
        if sample not in samples:
            samples[sample] = {}
        samples[sample]['intron'] = f
    elif '_exon_counts.txt' in f:
        sample = f.replace('_exon_counts.txt', '')
        if sample not in samples:
            samples[sample] = {}
        samples[sample]['exon'] = f
    elif '_5utr_counts.txt' in f:
        sample = f.replace('_5utr_counts.txt', '')
        if sample not in samples:
            samples[sample] = {}
        samples[sample]['5utr'] = f
    elif '_3utr_counts.txt' in f:
        sample = f.replace('_3utr_counts.txt', '')
        if sample not in samples:
            samples[sample] = {}
        samples[sample]['3utr'] = f

def count_reads(filename):
    # Count total reads in a count file (sum of second column)
    total = 0
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('gene_id'):
                    continue
                parts = line.split('\\t')
                if len(parts) >= 2:
                    try:
                        total += int(float(parts[1]))
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Error reading {filename}: {e}")
    return total

# Create report
with open('genome_read_counts.txt', 'w') as out:
    out.write("=" * 80 + "\\n")
    out.write("GENOME READ COUNTS REPORT\\n")
    out.write("=" * 80 + "\\n\\n")
    out.write("This report shows the total number of reads mapped to different\\n")
    out.write("genomic features (exons, introns, UTRs, etc.) for each sample.\\n\\n")
    
    # Table header
    out.write("-" * 80 + "\\n")
    utr5_label = "5'UTR"
    utr3_label = "3'UTR"
    out.write(f"{'Sample':<25} {'Exon':<12} {'Intron':<12} {utr5_label:<12} {utr3_label:<12} {'Transcript':<12}\\n")
    out.write("-" * 80 + "\\n")
    
    # Process each sample
    for sample in sorted(samples.keys()):
        files = samples[sample]
        
        exon_count = count_reads(files.get('exon', '')) if 'exon' in files else 0
        intron_count = count_reads(files.get('intron', '')) if 'intron' in files else 0
        utr5_count = count_reads(files.get('5utr', '')) if '5utr' in files else 0
        utr3_count = count_reads(files.get('3utr', '')) if '3utr' in files else 0
        transcript_count = count_reads(files.get('transcript', '')) if 'transcript' in files else 0
        
        out.write(f"{sample:<25} {exon_count:<12,} {intron_count:<12,} {utr5_count:<12,} {utr3_count:<12,} {transcript_count:<12,}\\n")
    
    out.write("-" * 80 + "\\n")
    out.write("\\n")
    out.write("Note: Counts represent the total number of reads assigned to each feature type.\\n")
    out.write("      - Exon: Reads mapped to exonic regions\\n")
    out.write("      - Intron: Reads mapped to intronic regions\\n")
    out.write("      - 5'UTR: Reads mapped to 5' untranslated regions\\n")
    out.write("      - 3'UTR: Reads mapped to 3' untranslated regions\\n")
    out.write("      - Transcript: Reads mapped to full transcript regions\\n")
    out.write("\\n")
    out.write("=" * 80 + "\\n")
    out.write("END OF REPORT\\n")
    out.write("=" * 80 + "\\n")

print(f"Processed {len(samples)} samples")
EOF

    python3 create_report.py

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //g')
END_VERSIONS
    """
}
