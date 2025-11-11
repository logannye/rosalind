#!/bin/bash

# Test runner script for sqrt-space-sim
# Runs all tests and provides clear output

set -e  # Exit on error

echo "=========================================="
echo "sqrt-space-sim Test Runner"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test counters
PASSED=0
FAILED=0
TOTAL=0

# Function to run a test suite
run_test_suite() {
    local suite_name=$1
    local test_pattern=$2
    
    echo -e "${BLUE}Running $suite_name tests...${NC}"
    echo "----------------------------------------"
    
    if cargo test --test "$test_pattern" -- --nocapture 2>&1 | tee /tmp/test_output_$$.log; then
        local passed=$(grep -c "test result: ok" /tmp/test_output_$$.log || echo "0")
        local failed=$(grep -c "test result: FAILED" /tmp/test_output_$$.log || echo "0")
        
        PASSED=$((PASSED + passed))
        FAILED=$((FAILED + failed))
        TOTAL=$((TOTAL + passed + failed))
        
        if [ "$failed" -eq 0 ]; then
            echo -e "${GREEN}✓ $suite_name: All tests passed${NC}"
        else
            echo -e "${RED}✗ $suite_name: $failed test(s) failed${NC}"
        fi
    else
        echo -e "${RED}✗ $suite_name: Test suite failed${NC}"
        FAILED=$((FAILED + 1))
        TOTAL=$((TOTAL + 1))
    fi
    
    echo ""
}

# Check if cargo is available
if ! command -v cargo &> /dev/null; then
    echo -e "${RED}Error: cargo not found. Please install Rust.${NC}"
    exit 1
fi

# Build the project first
echo -e "${BLUE}Building project...${NC}"
if ! cargo build --lib 2>&1 | grep -q "Finished"; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi
echo -e "${GREEN}Build successful${NC}"
echo ""

# Run unit tests (from src/)
echo -e "${BLUE}Running unit tests...${NC}"
echo "----------------------------------------"
if cargo test --lib -- --nocapture 2>&1 | tee /tmp/unit_test_output_$$.log; then
    unit_passed=$(grep -c "test result: ok" /tmp/unit_test_output_$$.log || echo "0")
    unit_failed=$(grep -c "test result: FAILED" /tmp/unit_test_output_$$.log || echo "0")
    PASSED=$((PASSED + unit_passed))
    FAILED=$((FAILED + unit_failed))
    TOTAL=$((TOTAL + unit_passed + unit_failed))
    
    if [ "$unit_failed" -eq 0 ]; then
        echo -e "${GREEN}✓ Unit tests: All passed${NC}"
    else
        echo -e "${RED}✗ Unit tests: $unit_failed test(s) failed${NC}"
    fi
else
    echo -e "${RED}✗ Unit tests: Failed${NC}"
    FAILED=$((FAILED + 1))
    TOTAL=$((TOTAL + 1))
fi
echo ""

# Run integration tests
run_test_suite "Correctness" "correctness"
run_test_suite "Space Bounds" "space_bounds"
run_test_suite "Integration" "integration"

# Summary
echo "=========================================="
echo -e "${BLUE}Test Summary${NC}"
echo "=========================================="
echo -e "Total tests: $TOTAL"
echo -e "${GREEN}Passed: $PASSED${NC}"
if [ "$FAILED" -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED${NC}"
else
    echo -e "${GREEN}Failed: $FAILED${NC}"
fi
echo ""

# Clean up temp files
rm -f /tmp/test_output_$$.log /tmp/unit_test_output_$$.log

# Exit with appropriate code
if [ "$FAILED" -gt 0 ]; then
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
else
    echo -e "${GREEN}All tests passed! ✓${NC}"
    exit 0
fi

