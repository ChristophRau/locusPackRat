#!/bin/bash
#SBATCH --job-name=packrat_integration
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=/proj/raulab/users/brian/packrat/tests/slurm/test_results/integration_%j.out
#SBATCH --error=/proj/raulab/users/brian/packrat/tests/slurm/test_results/integration_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bgural@unc.edu

# Integration tests for PackRat
# These tests use larger datasets and full workflows

echo "=================================="
echo "PackRat Integration Tests"
echo "=================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Load required modules
module load r/4.4.0

# Change to package directory
cd /proj/raulab/users/brian/packrat

# Create output directory if it doesn't exist
mkdir -p tests/integration/output

# Run integration tests
echo "Running integration tests..."
Rscript -e "
  # Load package
  library(PackRat)

  # Source integration test file
  source('tests/integration/test_full_pipeline.R')

  # Run tests
  tryCatch({
    run_integration_tests()
    cat('\n✓ All integration tests passed\n')
  }, error = function(e) {
    cat('\n✗ Integration tests failed:\n')
    print(e)
    quit(status = 1)
  })
"

echo ""
echo "Integration tests completed at $(date)"
echo "Total runtime: $SECONDS seconds"