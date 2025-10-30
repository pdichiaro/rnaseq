#!/bin/bash

echo "=== Checking Parameter Usage in nf-core/rnaseq Pipeline ==="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Arrays to track results
unused_params=()
used_params=()

echo "Analyzing parameter usage..."
echo ""

while IFS= read -r param; do
    param_trimmed=$(echo "$param" | xargs)  # Remove whitespace
    
    if [[ -z "$param_trimmed" ]]; then
        continue
    fi
    
    # Search for parameter usage in all .nf files and .config files
    usage_count=$(grep -r "params\.${param_trimmed}\b" \
        --include="*.nf" \
        --include="*.config" \
        --exclude="nextflow.config" \
        . 2>/dev/null | wc -l)
    
    # Also check for usage in test configs
    test_usage=$(grep -r "params\.${param_trimmed}\b" \
        --include="*.config" \
        conf/ 2>/dev/null | wc -l)
    
    total_usage=$((usage_count + test_usage))
    
    if [[ $total_usage -eq 0 ]]; then
        unused_params+=("$param_trimmed")
    else
        used_params+=("$param_trimmed")
    fi
    
done < all_params.txt

echo -e "${GREEN}✓ Used Parameters (${#used_params[@]}):${NC}"
for param in "${used_params[@]}"; do
    echo "  ✓ $param"
done

echo ""
echo -e "${RED}✗ Potentially Unused Parameters (${#unused_params[@]}):${NC}"
for param in "${unused_params[@]}"; do
    echo "  ✗ $param"
done

echo ""
echo "=== Summary ==="
echo "Total parameters defined: $((${#used_params[@]} + ${#unused_params[@]}))"
echo -e "Used parameters: ${GREEN}${#used_params[@]}${NC}"
echo -e "Potentially unused parameters: ${RED}${#unused_params[@]}${NC}"

# Additional check for hardcoded values that could be parameterized
echo ""
echo "=== Checking for potential hardcoded values ==="

# Look for common patterns that might benefit from parameterization
echo "Searching for potential hardcoded values..."
grep -r -n --include="*.nf" -E "(cpus.*[0-9]+|memory.*[0-9]+|time.*[0-9]+)" . | head -10