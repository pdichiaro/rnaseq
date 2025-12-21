#!/usr/bin/env python3

"""
Format genome counts summary files for MultiQC custom content.
Converts genome_count summary into MultiQC-compatible format showing:
- Read distribution across genomic features (exon, intron, UTRs)
- Genes detected per feature type
- General statistics for main report table
"""

import argparse
import sys
import os

def parse_genome_count_summary(summary_files, sample_ids=None):
    """
    Parse one or more genome count summary files.
    
    Args:
        summary_files: List of summary file paths
        sample_ids: Optional list of sample IDs (if not provided, extract from files)
    
    Returns:
        Dictionary with sample data
    """
    all_data = {}
    
    for i, summary_file in enumerate(summary_files):
        # Get sample ID
        if sample_ids and i < len(sample_ids):
            sample_id = sample_ids[i]
        else:
            # Extract from filename: SAMPLE_summary.txt
            sample_id = os.path.basename(summary_file).replace('_summary.txt', '')
        
        # Parse the summary file
        feature_data = {}
        try:
            with open(summary_file, 'r') as f:
                # Skip header
                header = f.readline().strip().split('\t')
                
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            sample = parts[0]
                            feature_type = parts[1]
                            total_reads = int(parts[2])
                            genes_with_reads = int(parts[3])
                            
                            feature_data[feature_type] = {
                                'total_reads': total_reads,
                                'genes_detected': genes_with_reads
                            }
            
            all_data[sample_id] = feature_data
            
        except Exception as e:
            print(f"Error parsing {summary_file}: {e}", file=sys.stderr)
            continue
    
    return all_data


def write_multiqc_general_stats(all_data, output_file):
    """
    Write general statistics table for MultiQC main table.
    Shows % exonic, % intronic, and total genes detected.
    """
    
    mqc_header = """# id: 'genome_counts_general_stats'
# section_name: 'Genome Feature Counts'
# description: 'Read distribution across genomic features from genome-based quantification'
# plot_type: 'generalstats'
# pconfig:
#     percent_exonic:
#         title: '% Exonic'
#         description: 'Percentage of reads mapping to exons'
#         namespace: 'Genome Counts'
#         max: 100
#         min: 0
#         scale: 'RdYlGn'
#         format: '{:.1f}%'
#     percent_intronic:
#         title: '% Intronic'
#         description: 'Percentage of reads mapping to introns'
#         namespace: 'Genome Counts'
#         max: 100
#         min: 0
#         scale: 'RdYlGn-rev'
#         format: '{:.1f}%'
#     exon_intron_ratio:
#         title: 'Exon/Intron Ratio'
#         description: 'Ratio of exonic to intronic reads (higher is better)'
#         namespace: 'Genome Counts'
#         min: 0
#         scale: 'Blues'
#         format: '{:.2f}'
#     genes_detected_exon:
#         title: 'Genes (Exon)'
#         description: 'Number of genes with exonic reads'
#         namespace: 'Genome Counts'
#         format: '{:,.0f}'
"""
    
    with open(output_file, 'w') as out:
        out.write(mqc_header + '\n')
        
        # Write header
        out.write('Sample\tpercent_exonic\tpercent_intronic\texon_intron_ratio\tgenes_detected_exon\n')
        
        # Write data for each sample
        for sample_id, features in all_data.items():
            # Calculate totals
            total_reads = sum(f['total_reads'] for f in features.values())
            
            if total_reads > 0:
                exon_reads = features.get('exon', {}).get('total_reads', 0)
                intron_reads = features.get('intron', {}).get('total_reads', 0)
                
                percent_exonic = (exon_reads / total_reads) * 100
                percent_intronic = (intron_reads / total_reads) * 100
                
                # Calculate exon/intron ratio (avoid division by zero)
                if intron_reads > 0:
                    exon_intron_ratio = exon_reads / intron_reads
                else:
                    exon_intron_ratio = exon_reads  # If no intron reads, use exon count
                
                genes_exon = features.get('exon', {}).get('genes_detected', 0)
                
                out.write(f"{sample_id}\t{percent_exonic:.2f}\t{percent_intronic:.2f}\t"
                         f"{exon_intron_ratio:.2f}\t{genes_exon}\n")


def write_multiqc_bargraph(all_data, output_file):
    """
    Write bargraph data showing read distribution across all features.
    """
    
    mqc_header = """# id: 'genome_feature_distribution'
# section_name: 'Genomic Feature Distribution'
# description: 'Distribution of reads across different genomic features (transcript, exon, intron, UTRs)'
# plot_type: 'bargraph'
# pconfig:
#     id: 'genome_feature_reads_plot'
#     title: 'Reads per Genomic Feature'
#     ylab: 'Number of Reads'
#     cpswitch_counts_label: 'Number of Reads'
"""
    
    with open(output_file, 'w') as out:
        out.write(mqc_header + '\n')
        
        # Get all feature types (should be consistent across samples)
        feature_types = ['transcript', 'exon', 'intron', '5utr', '3utr']
        
        # Write header
        out.write('Sample\t' + '\t'.join(feature_types) + '\n')
        
        # Write data for each sample
        for sample_id, features in all_data.items():
            values = []
            for ft in feature_types:
                reads = features.get(ft, {}).get('total_reads', 0)
                values.append(str(reads))
            
            out.write(f"{sample_id}\t" + '\t'.join(values) + '\n')


def write_multiqc_table(all_data, output_file):
    """
    Write detailed table showing genes detected per feature.
    """
    
    mqc_header = """# id: 'genome_genes_detected'
# section_name: 'Genes Detected per Feature'
# description: 'Number of genes with at least one read in each genomic feature type'
# plot_type: 'table'
# pconfig:
#     id: 'genome_genes_table'
#     title: 'Genes Detected by Feature Type'
"""
    
    with open(output_file, 'w') as out:
        out.write(mqc_header + '\n')
        
        feature_types = ['transcript', 'exon', 'intron', '5utr', '3utr']
        
        # Write header
        out.write('Sample\t' + '\t'.join([f'{ft}_genes' for ft in feature_types]) + '\n')
        
        # Write data for each sample
        for sample_id, features in all_data.items():
            values = []
            for ft in feature_types:
                genes = features.get(ft, {}).get('genes_detected', 0)
                values.append(str(genes))
            
            out.write(f"{sample_id}\t" + '\t'.join(values) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Format genome counts summary files for MultiQC custom content'
    )
    parser.add_argument(
        'summary_files',
        nargs='+',
        help='One or more genome count summary files (*_summary.txt)'
    )
    parser.add_argument(
        '-o', '--output-prefix',
        default='genome_counts',
        help='Output file prefix (default: genome_counts)'
    )
    parser.add_argument(
        '-s', '--sample-ids',
        nargs='+',
        help='Optional sample IDs (if not provided, extracted from filenames)'
    )
    
    args = parser.parse_args()
    
    # Parse all summary files
    all_data = parse_genome_count_summary(args.summary_files, args.sample_ids)
    
    if not all_data:
        print("Error: No data parsed from summary files", file=sys.stderr)
        sys.exit(1)
    
    # Write MultiQC outputs
    write_multiqc_general_stats(all_data, f"{args.output_prefix}_generalstats_mqc.tsv")
    write_multiqc_bargraph(all_data, f"{args.output_prefix}_distribution_mqc.tsv")
    write_multiqc_table(all_data, f"{args.output_prefix}_genes_mqc.tsv")
    
    print(f"Generated MultiQC files for {len(all_data)} samples")


if __name__ == '__main__':
    main()
