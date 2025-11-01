#!/bin/bash
# Quick tests that are safe to run on login node
# These are small unit tests with minimal data

cd /proj/raulab/users/brian/packrat

echo "=================================="
echo "PackRat Quick Tests (Login Safe)"
echo "=================================="
echo "Date: $(date)"
echo "Directory: $(pwd)"
echo ""

# Check if R is available
if ! command -v R &> /dev/null; then
    echo "Error: R is not available. Please load R module:"
    echo "  module load r/4.4.0"
    exit 1
fi

echo "Running unit tests with testthat..."
Rscript -e "
  suppressPackageStartupMessages({
    library(testthat)
  })

  # Check if package can be loaded
  tryCatch({
    library(PackRat)
    cat('✓ Package loaded successfully\n\n')
  }, error = function(e) {
    cat('✗ Failed to load package. Installing...\n')
    devtools::install('.', quiet = TRUE)
  })

  # Run tests
  results <- test_dir('tests/testthat', reporter = 'summary')

  # Summary
  cat('\n================================\n')
  cat('Test Summary:\n')
  cat(sprintf('Passed: %d\n', sum(results\$passed)))
  cat(sprintf('Failed: %d\n', sum(results\$failed)))
  cat(sprintf('Skipped: %d\n', sum(results\$skipped)))

  # Exit with error if tests failed
  if (sum(results\$failed) > 0) {
    quit(status = 1)
  }
"

echo ""
echo "Quick tests completed at $(date)"