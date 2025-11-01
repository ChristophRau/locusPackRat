# Makefile for PackRat package development and testing
# Run with: make [target]

.PHONY: help test test-quick test-integration test-deps test-all \
        install clean check docs build setup-minimal

# Default target - show help
help:
	@echo "PackRat Development Makefile"
	@echo "============================"
	@echo ""
	@echo "Available targets:"
	@echo "  make test-quick     - Run quick unit tests (login node safe)"
	@echo "  make test-deps      - Test dependency installation from scratch"
	@echo "  make test-integration - Submit integration tests to SLURM"
	@echo "  make test-all       - Run all tests"
	@echo "  make setup-minimal  - Create minimal R environment for testing"
	@echo "  make install        - Install package locally"
	@echo "  make check          - Run R CMD check"
	@echo "  make docs           - Generate documentation with roxygen2"
	@echo "  make build          - Build package tarball"
	@echo "  make clean          - Clean test outputs and build artifacts"
	@echo ""
	@echo "Testing workflow:"
	@echo "  1. make test-quick  (for quick validation)"
	@echo "  2. make test-deps   (test dependency handling)"
	@echo "  3. make test-integration (for thorough testing)"

# Quick tests - safe for login node
test-quick:
	@echo "Running quick tests (login node safe)..."
	@bash tests/run_quick_tests.sh

# Test dependency installation from minimal environment
test-deps:
	@echo "Testing dependency installation from minimal R environment..."
	@echo "Note: This requires the packrat-minimal conda environment."
	@echo "If not created, run: make setup-minimal"
	@Rscript tests/test_dependency_install.R

# Setup minimal R environment
setup-minimal:
	@echo "Setting up minimal R environment for dependency testing..."
	@bash tests/setup_minimal.sh

# Integration tests - requires compute node
test-integration:
	@echo "Submitting integration tests to SLURM..."
	@sbatch tests/slurm/run_integration_tests.sh
	@echo "Job submitted. Check status with: squeue -u $$USER"
	@echo "View results in: tests/slurm/test_results/"

# Run all tests
test-all: test-quick
	@echo ""
	@echo "Quick tests complete. Submitting integration tests..."
	@$(MAKE) test-integration

# Alias for test
test: test-quick

# Install package locally
install:
	@echo "Installing PackRat package..."
	@R CMD INSTALL . --no-multiarch --with-keep.source

# R CMD check
check:
	@echo "Running R CMD check..."
	@R CMD build . --no-build-vignettes
	@R CMD check PackRat_*.tar.gz --no-manual --no-vignettes
	@rm -f PackRat_*.tar.gz

# Generate documentation with roxygen2
docs:
	@echo "Generating documentation..."
	@Rscript -e "devtools::document()"

# Build package
build: docs
	@echo "Building package..."
	@R CMD build .

# Clean test outputs and build artifacts
clean:
	@echo "Cleaning test outputs and build artifacts..."
	@rm -f tests/slurm/test_results/*.out
	@rm -f tests/slurm/test_results/*.err
	@rm -rf tests/integration/output/*
	@rm -rf PackRat.Rcheck/
	@rm -f PackRat_*.tar.gz
	@rm -f .Rhistory .RData
	@echo "Clean complete"

# Interactive development session
dev:
	@echo "Starting interactive R session for development..."
	@echo "Run devtools::load_all() to load the package"
	@R

# Check SLURM job status
status:
	@echo "Current SLURM jobs:"
	@squeue -u $$USER --name="packrat*"
	@echo ""
	@echo "Recent test results:"
	@ls -lht tests/slurm/test_results/*.out 2>/dev/null | head -5 || echo "No test results found"

# View latest test output
view-results:
	@if [ -f "$$(ls -t tests/slurm/test_results/*.out 2>/dev/null | head -1)" ]; then \
		tail -50 $$(ls -t tests/slurm/test_results/*.out | head -1); \
	else \
		echo "No test results found"; \
	fi

# Setup renv for dependency management
setup-renv:
	@echo "Setting up renv for dependency management..."
	@Rscript -e "if (!requireNamespace('renv')) install.packages('renv'); renv::init()"

# Update dependencies
update-deps:
	@echo "Updating package dependencies..."
	@Rscript -e "renv::snapshot()"