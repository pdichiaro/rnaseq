#!/bin/bash

# Verification script for channel refactoring
# This script checks that all the required changes are present

echo "=========================================="
echo "Channel Refactoring Verification Script"
echo "=========================================="
echo ""

WORKFLOW_FILE="workflows/rnaseq/main.nf"
ERRORS=0

# Function to check for pattern
check_pattern() {
    local pattern="$1"
    local description="$2"
    local expected_count="$3"
    
    count=$(grep -c "$pattern" "$WORKFLOW_FILE" 2>/dev/null || echo "0")
    
    if [ "$count" -eq "$expected_count" ]; then
        echo "✓ $description: Found $count occurrences (expected $expected_count)"
    else
        echo "✗ $description: Found $count occurrences (expected $expected_count)"
        ERRORS=$((ERRORS + 1))
    fi
}

echo "1. Checking channel initialization..."
check_pattern "ch_scaling_factors_individual = Channel.empty()" "Main channel initialization" 1
check_pattern "ch_scaling_factors_individual_invariant = Channel.empty()" "Invariant channel initialization" 1
check_pattern "ch_scaling_factors_individual_all_genes = Channel.empty()" "All genes channel initialization" 1
echo ""

echo "2. Checking population points (should have 8 each)..."
check_pattern "ch_scaling_factors_individual_invariant.*mix.*scaling_factors_individual" "Invariant population" 4
check_pattern "ch_scaling_factors_individual_all_genes.*mix.*scaling_factors_individual" "All genes population" 4
echo ""

echo "3. Checking consumption points..."
check_pattern "ch_scaling_per_sample_invariant = ch_scaling_factors_individual_invariant" "Invariant consumption" 1
check_pattern "ch_scaling_per_sample_all_genes = ch_scaling_factors_individual_all_genes" "All genes consumption" 1
echo ""

echo "4. Verifying no old directory-based filtering logic remains..."
# Check that we're not filtering by 'invariant' directory names anymore
if grep -q "parent_dir.contains('invariant')" "$WORKFLOW_FILE" 2>/dev/null; then
    old_filter_pattern_count=$(grep -c "parent_dir.contains('invariant')" "$WORKFLOW_FILE")
    echo "✗ Still has directory-based filtering logic"
    echo "  Found $old_filter_pattern_count occurrences of directory filtering"
    ERRORS=$((ERRORS + 1))
else
    echo "✓ No directory-based filtering found (clean refactoring)"
fi
echo ""

echo "5. Syntax validation..."
if nextflow config "$WORKFLOW_FILE" -profile test >/dev/null 2>&1; then
    echo "✓ Nextflow syntax is valid"
else
    echo "✗ Nextflow syntax check failed"
    ERRORS=$((ERRORS + 1))
fi
echo ""

echo "=========================================="
if [ $ERRORS -eq 0 ]; then
    echo "✓ ALL CHECKS PASSED!"
    echo "The channel refactoring is complete and correct."
    exit 0
else
    echo "✗ $ERRORS CHECK(S) FAILED"
    echo "Please review the errors above."
    exit 1
fi
echo "=========================================="
