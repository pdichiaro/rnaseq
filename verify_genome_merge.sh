#!/bin/bash

# Genome Quantification Merge Verification Script
# Usage: bash verify_genome_merge.sh [results_directory]

RESULTS_DIR="${1:-results}"
OUTDIR="$RESULTS_DIR/star/genome"

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║   Genome Quantification Merge Verification                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

pass_count=0
fail_count=0
warn_count=0

# Function to print status
print_pass() {
    echo -e "${GREEN}✓${NC} $1"
    ((pass_count++))
}

print_fail() {
    echo -e "${RED}✗${NC} $1"
    ((fail_count++))
}

print_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
    ((warn_count++))
}

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1. DIRECTORY STRUCTURE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -d "$OUTDIR" ]; then
    print_pass "Genome output directory exists: $OUTDIR"
    echo
    echo "Contents:"
    ls -lh "$OUTDIR/" | head -25
else
    print_fail "Genome output directory not found: $OUTDIR"
    echo "Please check your results directory path."
    exit 1
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "2. PER-SAMPLE FOLDER ORGANIZATION"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Find sample folders (exclude deseq2 and deeptools)
sample_folders=$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -type d ! -name "deseq2" ! -name "deeptools" 2>/dev/null)

if [ -n "$sample_folders" ]; then
    sample_count=$(echo "$sample_folders" | wc -l)
    print_pass "Found $sample_count sample folder(s)"
    echo
    echo "Sample folders:"
    echo "$sample_folders" | while read folder; do
        sample_name=$(basename "$folder")
        file_count=$(ls -1 "$folder" 2>/dev/null | wc -l)
        echo "  • $sample_name/ ($file_count files)"
    done
else
    print_fail "No sample folders found"
    print_warn "Expected: genome/{sample_id}/ folders for per-sample results"
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "3. SAMPLE FOLDER CONTENTS (FIRST SAMPLE)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

first_sample=$(echo "$sample_folders" | head -1)
if [ -n "$first_sample" ]; then
    sample_name=$(basename "$first_sample")
    echo "Checking: $sample_name/"
    echo
    
    expected_files=(
        "${sample_name}_exon_counts.txt"
        "${sample_name}_transcript_counts.txt"
        "${sample_name}_intron_counts.txt"
        "${sample_name}_5utr_counts.txt"
        "${sample_name}_3utr_counts.txt"
        "${sample_name}_combined_counts.txt"
        "${sample_name}_summary.txt"
    )
    
    for file in "${expected_files[@]}"; do
        if [ -f "$first_sample/$file" ]; then
            size=$(stat -f%z "$first_sample/$file" 2>/dev/null || stat -c%s "$first_sample/$file" 2>/dev/null)
            if [ "$size" -gt 0 ]; then
                print_pass "$file (${size} bytes)"
            else
                print_fail "$file (EMPTY FILE)"
            fi
        else
            print_fail "$file (MISSING)"
        fi
    done
else
    print_warn "No sample folder to check"
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "4. MERGED COUNT FILES AT ROOT"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

expected_merged=(
    "genome_exon_counts_merged.txt"
    "genome_transcript_counts_merged.txt"
    "genome_intron_counts_merged.txt"
    "genome_5utr_counts_merged.txt"
    "genome_3utr_counts_merged.txt"
)

for file in "${expected_merged[@]}"; do
    if [ -f "$OUTDIR/$file" ]; then
        size=$(stat -f%z "$OUTDIR/$file" 2>/dev/null || stat -c%s "$OUTDIR/$file" 2>/dev/null)
        if [ "$size" -gt 0 ]; then
            print_pass "$file (${size} bytes)"
        else
            print_fail "$file (EMPTY FILE)"
        fi
    else
        print_fail "$file (MISSING)"
    fi
done

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "5. MERGE SUMMARY FILES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

expected_summaries=(
    "genome_exon_merge_summary.txt"
    "genome_transcript_merge_summary.txt"
    "genome_intron_merge_summary.txt"
    "genome_5utr_merge_summary.txt"
    "genome_3utr_merge_summary.txt"
)

for file in "${expected_summaries[@]}"; do
    if [ -f "$OUTDIR/$file" ]; then
        size=$(stat -f%z "$OUTDIR/$file" 2>/dev/null || stat -c%s "$OUTDIR/$file" 2>/dev/null)
        if [ "$size" -gt 0 ]; then
            print_pass "$file (${size} bytes)"
        else
            print_warn "$file (EMPTY FILE)"
        fi
    else
        print_warn "$file (MISSING)"
    fi
done

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "6. MERGED FILE DIMENSIONS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

for file in "${expected_merged[@]}"; do
    if [ -f "$OUTDIR/$file" ]; then
        lines=$(wc -l < "$OUTDIR/$file")
        cols=$(head -n 1 "$OUTDIR/$file" | awk '{print NF}')
        genes=$((lines - 1))  # Subtract header
        samples=$((cols - 6))  # Subtract 6 fixed columns
        
        echo "• $(basename $file):"
        echo "    Genes: $genes"
        echo "    Samples: $samples"
        echo "    Total columns: $cols"
        
        if [ "$genes" -gt 0 ] && [ "$samples" -gt 0 ]; then
            print_pass "Valid matrix dimensions"
        else
            print_fail "Invalid matrix dimensions"
        fi
        echo
    fi
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "7. SAMPLE NAMES IN MERGED FILES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f "$OUTDIR/genome_exon_counts_merged.txt" ]; then
    echo "Header from exon counts merged file:"
    head -n 1 "$OUTDIR/genome_exon_counts_merged.txt"
    echo
    
    # Extract sample names (columns 7+)
    header=$(head -n 1 "$OUTDIR/genome_exon_counts_merged.txt")
    sample_names=$(echo "$header" | awk '{for(i=7;i<=NF;i++) print $i}')
    
    echo "Samples in merged file:"
    echo "$sample_names" | while read name; do
        echo "  • $name"
    done
    
    # Check if sample folders match
    echo
    echo "Checking if sample folders match merged file samples..."
    echo "$sample_names" | while read name; do
        if [ -d "$OUTDIR/$name" ]; then
            print_pass "Folder exists for sample: $name"
        else
            print_warn "No folder found for sample: $name"
        fi
    done
else
    print_fail "Cannot check sample names - merged file not found"
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "8. DESEQ2 QC OUTPUTS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

deseq2_qc_dir="$OUTDIR/deseq2/all_genes/Quality_Control"
if [ -d "$deseq2_qc_dir" ]; then
    print_pass "DESeq2 QC directory exists"
    
    pdf_count=$(ls -1 "$deseq2_qc_dir"/*.pdf 2>/dev/null | wc -l)
    if [ "$pdf_count" -gt 0 ]; then
        print_pass "Found $pdf_count QC PDF file(s)"
        ls -1 "$deseq2_qc_dir"/*.pdf 2>/dev/null | while read pdf; do
            echo "  • $(basename $pdf)"
        done
    else
        print_warn "No QC PDF files found"
    fi
else
    print_warn "DESeq2 QC directory not found (may be disabled)"
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "9. DEEPTOOLS BIGWIG FILES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

deeptools_dir="$OUTDIR/deeptools/all_genes"
if [ -d "$deeptools_dir" ]; then
    print_pass "DeepTools directory exists"
    
    bw_count=$(ls -1 "$deeptools_dir"/*.bw 2>/dev/null | wc -l)
    if [ "$bw_count" -gt 0 ]; then
        print_pass "Found $bw_count BigWig file(s)"
        ls -1 "$deeptools_dir"/*.bw 2>/dev/null | while read bw; do
            echo "  • $(basename $bw)"
        done
    else
        print_warn "No BigWig files found"
    fi
else
    print_warn "DeepTools directory not found (may be disabled)"
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "10. STRUCTURE CONSISTENCY CHECK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Count items at root level (should be clean)
root_items=$(ls -1 "$OUTDIR" | wc -l)
sample_folders_count=$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -type d ! -name "deseq2" ! -name "deeptools" 2>/dev/null | wc -l)
merged_files_count=$(ls -1 "$OUTDIR"/*.txt 2>/dev/null | wc -l)

echo "Root directory organization:"
echo "  • Sample folders: $sample_folders_count"
echo "  • Merged files: $merged_files_count"
echo "  • Total items: $root_items"
echo

if [ "$root_items" -lt 30 ]; then
    print_pass "Root directory is well organized (< 30 items)"
else
    print_warn "Root directory may be cluttered ($root_items items)"
fi

# Check if RSEM exists for comparison
rsem_dir="$RESULTS_DIR/star/rsem"
if [ -d "$rsem_dir" ]; then
    echo
    echo "Comparing with RSEM structure:"
    rsem_sample_folders=$(find "$rsem_dir" -mindepth 1 -maxdepth 1 -type d ! -name "deseq2" ! -name "deeptools" 2>/dev/null | wc -l)
    echo "  • RSEM sample folders: $rsem_sample_folders"
    echo "  • Genome sample folders: $sample_folders_count"
    
    if [ "$rsem_sample_folders" -eq "$sample_folders_count" ]; then
        print_pass "Consistent sample count with RSEM"
    else
        print_warn "Sample count mismatch with RSEM"
    fi
fi

echo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "VERIFICATION SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo
echo -e "${GREEN}✓ Passed:${NC} $pass_count"
echo -e "${YELLOW}⚠ Warnings:${NC} $warn_count"
echo -e "${RED}✗ Failed:${NC} $fail_count"
echo

if [ "$fail_count" -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  ✓ ALL CRITICAL CHECKS PASSED                                 ║${NC}"
    echo -e "${GREEN}║  Genome quantification merge completed successfully!          ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║  ✗ SOME CHECKS FAILED                                          ║${NC}"
    echo -e "${RED}║  Please review the errors above                                ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
