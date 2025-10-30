#!/bin/bash

echo "=== Detailed Parameter Usage Analysis ==="
echo ""

# List of potentially unused parameters to investigate
params=(
    "config_profile_contact"
    "config_profile_description" 
    "config_profile_name"
    "config_profile_url"
    "custom_config_base"
    "custom_config_version"
    "deseq2_vst"
    "extra_fastp_args"
    "extra_fqlint_args"
    "extra_kallisto_quant_args"
    "extra_salmon_quant_args"
    "extra_trimgalore_args"
    "featurecounts_feature_type"
    "help"
    "help_full"
    "igenomes_ignore"
    "multiqc_title"
    "pipelines_testdata_base_path"
    "pseudo_aligner_kmer_size"
    "publish_dir_mode"
    "save_bbsplit_reads"
    "save_merged_fastq"
    "save_non_ribo_reads"
    "save_reference"
    "save_umi_intermeds"
    "show_hidden"
    "skip_gtf_transcript_filter"
    "stringtie_ignore_gtf"
    "trace_report_suffix"
    "umitools_extract_method"
    "umitools_umi_separator"
)

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

truly_unused=()
actually_used=()

for param in "${params[@]}"; do
    echo "Checking: $param"
    
    # More comprehensive search including different patterns
    usage=$(grep -r -l "${param}" \
        --include="*.nf" \
        --include="*.config" \
        --include="*.groovy" \
        --include="*.yml" \
        --include="*.yaml" \
        . 2>/dev/null | grep -v "check_param_usage.sh" | grep -v "detailed_param_check.sh")
    
    # Also check for alternative patterns like params['param'] or variable assignments
    usage2=$(grep -r -E "(params\[.${param}.\]|${param}\s*=|\\\$\{${param}\})" \
        --include="*.nf" \
        --include="*.config" \
        . 2>/dev/null)
    
    if [[ -z "$usage" && -z "$usage2" ]]; then
        truly_unused+=("$param")
        echo -e "  ${RED}✗ Not found${NC}"
    else
        actually_used+=("$param")
        echo -e "  ${GREEN}✓ Found in:${NC}"
        if [[ -n "$usage" ]]; then
            echo "$usage" | sed 's/^/    /'
        fi
        if [[ -n "$usage2" ]]; then
            echo "  Usage patterns:"
            echo "$usage2" | head -3 | sed 's/^/    /'
        fi
    fi
    echo ""
done

echo "=== Final Analysis ==="
echo -e "Actually used parameters: ${GREEN}${#actually_used[@]}${NC}"
for param in "${actually_used[@]}"; do
    echo "  ✓ $param"
done

echo ""
echo -e "Truly unused parameters: ${RED}${#truly_unused[@]}${NC}"
for param in "${truly_unused[@]}"; do
    echo "  ✗ $param"
done