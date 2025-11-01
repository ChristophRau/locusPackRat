# PackRat Quick Start Guide

## Setup Development Environment

### Option 1: Conda Environment (Recommended)
```bash
# Create and activate conda environment
./setup_environment.sh
conda activate packrat-dev

# Test installation
Rscript tests/test_installation.R
```

### Option 2: System R with Modules
```bash
# Load R module
module load r/4.4.0

# Install in R
R
> install.packages(c("data.table", "dplyr", "openxlsx", "testthat"))
> devtools::install(".")
```

## Testing Your Changes

### Quick Tests (Login Node Safe)
```bash
make test-quick  # Run unit tests
```

### Full Tests (Requires Compute Node)
```bash
make test-integration  # Submit to SLURM
make status           # Check job status
make view-results     # View results
```

## Basic Usage

```r
library(PackRat)

# Load your data
genes <- read.csv("genes.csv")
gwasHits <- read.csv("gwas_hits.csv")

# Filter genes to keep only GWAS hits
candidates <- filterGenes(
  inputTable = genes,
  referenceTable = gwasHits,
  by = "gene_id",
  joinType = "semi"
)

# Add filters
filtered <- filterGenes(
  inputTable = candidates,
  filters = list(
    makeFilter("biotype", "==", "protein_coding"),
    makeFilter("expression", ">", 100)
  )
)

# Create Excel report
createGeneWorkbook(filtered, "results.xlsx")
```

## Development Workflow

1. **Make changes** to R files
2. **Test quickly**: `make test-quick`
3. **Test thoroughly**: `make test-integration`
4. **Build package**: `make build`
5. **Check package**: `make check`

## Package Structure

- `R/` - Source code (use camelCase naming)
- `tests/` - Test scripts and data
- `inst/testdata/` - Small example datasets
- `environment.yml` - Conda specification
- `Makefile` - Development commands

## Key Functions

- `filterGenes()` - Relational filtering with flexible joins
- `makeFilter()` - Create filter specifications
- `buildGeneTable()` - Combine annotations
- `createGeneWorkbook()` - Export to Excel
- `plotLocus()` - Create visualizations

## Help & Documentation

- Run `?filterGenes` in R for function help
- See `examples/` for usage examples
- Check `.claude_notes/` for development details

## Common Issues

### Dependency Conflicts
If you get library errors, ensure you're using EITHER conda OR modules, not both:
```bash
conda deactivate  # If using modules
# OR
module unload r   # If using conda
```

### Missing Packages
Install missing packages:
```r
install.packages(c("package1", "package2"))
```

### SLURM Jobs Not Running
Check queue and errors:
```bash
squeue -u $USER
cat tests/slurm/test_results/*.err
```