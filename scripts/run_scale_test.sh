#!/bin/bash

# Scale Performance Test Runner
# Runs the scale performance test and optionally exports to CSV

set -euo pipefail

tmp_output="$(mktemp)"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}Scale Performance Test: O(√t) Space Complexity Verification${NC}"
echo -e "${CYAN}================================================================${NC}\n"

cleanup() {
    rm -f "$tmp_output"
}
trap cleanup EXIT

# Check if CSV export is requested
if [ "${1:-}" = "--csv" ]; then
    echo -e "${BLUE}Running with CSV export enabled...${NC}\n"
    SCALE_TEST_CSV=1 cargo run --example scale_performance_test 2>&1 \
        | tee scale_test_results.txt \
        | tee "$tmp_output"
    echo -e "\n${GREEN}✓ Results saved to scale_test_results.txt${NC}"
    echo -e "${GREEN}✓ CSV data included in output${NC}"
else
    echo -e "${BLUE}Running scale performance test...${NC}\n"
    cargo run --example scale_performance_test 2>&1 | tee "$tmp_output"
    echo -e "\n${BLUE}Tip: Use --csv flag to export results to CSV format${NC}"
    echo -e "${BLUE}Example: $0 --csv${NC}"
fi

if ! grep -q "✓ All tests satisfy space bound" "$tmp_output"; then
    echo "Space bound verification missing from output" >&2
    exit 1
fi

if ! grep -q "✓ Space scaling:" "$tmp_output" || ! grep -q "Scaling is sublinear" "$tmp_output"; then
    echo "Sublinear scaling verification missing from output" >&2
    exit 1
fi

