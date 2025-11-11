#!/bin/bash

# Comprehensive test runner for sqrt-space-sim
# Provides clear output and summary

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}╔════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║   sqrt-space-sim Test Suite           ║${NC}"
echo -e "${CYAN}╚════════════════════════════════════════╝${NC}"
echo ""

# Check cargo
if ! command -v cargo &> /dev/null; then
    echo -e "${RED}✗ Error: cargo not found${NC}"
    exit 1
fi

# Build first
echo -e "${BLUE}Building project...${NC}"
if cargo build --lib --quiet 2>&1; then
    echo -e "${GREEN}✓ Build successful${NC}"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi
echo ""

# Test results
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Function to run a test and capture results
run_test() {
    local test_name=$1
    local test_filter=$2
    local test_args=$3  # Additional cargo test arguments
    
    TESTS_RUN=$((TESTS_RUN + 1))
    echo -e "${BLUE}Running: ${test_name}${NC}"
    
    # Build command - handle empty test filter for --lib
    local cmd="cargo test $test_args"
    if [ -n "$test_filter" ]; then
        cmd="$cmd --test $test_filter"
    fi
    cmd="$cmd -- --nocapture"
    
    if eval "$cmd" 2>&1 | tee /tmp/test_$$.log; then
        if grep -q "test result: ok" /tmp/test_$$.log; then
            local passed=$(grep "test result: ok" /tmp/test_$$.log | grep -oE "[0-9]+ passed" | grep -oE "[0-9]+" | head -1)
            local failed=$(grep "test result: ok" /tmp/test_$$.log | grep -oE "[0-9]+ failed" | grep -oE "[0-9]+" | head -1 || echo "0")
            
            if [ "${failed:-0}" = "0" ]; then
                echo -e "${GREEN}✓ ${test_name}: PASSED (${passed} tests)${NC}"
                TESTS_PASSED=$((TESTS_PASSED + 1))
            else
                echo -e "${YELLOW}⚠ ${test_name}: ${passed} passed, ${failed} failed${NC}"
                TESTS_FAILED=$((TESTS_FAILED + 1))
            fi
        else
            echo -e "${RED}✗ ${test_name}: FAILED${NC}"
            TESTS_FAILED=$((TESTS_FAILED + 1))
        fi
    else
        echo -e "${RED}✗ ${test_name}: ERROR${NC}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    echo ""
}

# Run all test suites
echo -e "${CYAN}═══════════════════════════════════════${NC}"
echo -e "${CYAN}Running Test Suites${NC}"
echo -e "${CYAN}═══════════════════════════════════════${NC}"
echo ""

run_test "Unit Tests" "" "--lib"
run_test "Correctness Tests" "correctness" ""
run_test "Space Bounds Tests" "space_bounds" ""
run_test "Integration Tests" "integration" ""
run_test "Algebra Tests" "algebra" ""
run_test "Ledger Tests" "ledger" ""

# Summary
echo -e "${CYAN}═══════════════════════════════════════${NC}"
echo -e "${CYAN}Test Summary${NC}"
echo -e "${CYAN}═══════════════════════════════════════${NC}"
echo ""
echo -e "Test Suites Run: ${TESTS_RUN}"
echo -e "${GREEN}Passed: ${TESTS_PASSED}${NC}"
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "${RED}Failed: ${TESTS_FAILED}${NC}"
else
    echo -e "${GREEN}Failed: ${TESTS_FAILED}${NC}"
fi
echo ""

# Cleanup
rm -f /tmp/test_$$.log

# Exit code
if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "${YELLOW}Note: Some tests may fail due to known limitations:${NC}"
    echo -e "${YELLOW}  - Block boundary reconstruction not fully implemented${NC}"
    echo -e "${YELLOW}  - Interface checking may fail for adjacent blocks${NC}"
    echo ""
    exit 1
else
    echo -e "${GREEN}All tests passed! ✓${NC}"
    exit 0
fi

