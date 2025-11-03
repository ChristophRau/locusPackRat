# Locus Zoom Plot Test Run Summary
**Date**: January 3, 2025

## Test Status: BLOCKED

### Issue
The tests cannot run because the required Bioconductor annotation package is not installed in the system R library.

**Missing package**: `TxDb.Mmusculus.UCSC.mm39.knownGene`

### What We Fixed
1. ✅ Created test fixtures with correct data structure:
   - `test_genes_mouse.csv` - 10 genes on chr1
   - `test_scan_data.rds` - QTL scan results with proper structure (LOD as named vector, pos as list with Mb element, allele.effects matrix)
   - `test_threshold_data.rds` - Thresholds with correct naming convention (`HR_Ctrl_threshold`)

2. ✅ Fixed data type issues:
   - Converted pos column to integer explicitly
   - Verified all data structures match expected format

3. ✅ Resolved library path conflicts:
   - System R library: `/nas/longleaf/rhel9/apps/r/4.4.0/lib64/R/library`
   - Successfully loaded plotgardener, data.table, dplyr, RColorBrewer

### Current Blocker

The `generateLocusZoomPlot()` function requires:
```r
assembly <- plotgardener::assembly(
  Genome = "mm39_GRCm39",
  TxDb   = "TxDb.Mmusculus.UCSC.mm39.knownGene",
  OrgDb  = "org.Mm.eg.db"
)
```

**Error**: `'assembly' not available as a default. Please make a assembly object with assembly() or pick an assembly from the defaults listed with genomes().`

The TxDb package is not installed in the system library.

### Next Steps

**Option 1**: Install TxDb in system library (requires admin/module maintainer)
```r
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene")
```

**Option 2**: Use mm10 assembly if TxDb.Mmusculus.UCSC.mm10.knownGene is available

**Option 3**: Run tests in original working environment where packages are installed

**Option 4**: Skip plot generation tests and focus on Excel/data tests first (test-packet-excel.R)

## Test Files Created

- ✅ `tests/testthat/fixtures/test_genes_mouse.csv`
- ✅ `tests/testthat/fixtures/test_scan_data.rds`
- ✅ `tests/testthat/fixtures/test_threshold_data.rds`
- ✅ `tests/testthat/fixtures/create_test_fixtures.R`
- ✅ `tests/testthat/test-packet-zoom-plot.R` (3 tests defined)
- ✅ `tests/run_zoom_test_system.R` (simplified test runner)

## Debugging Output

Successfully created test data frame for plotManhattan:
```
DEBUG: miqtl_df_for_plot dimensions: 101 x 6
DEBUG: Column classes: character, character, integer, numeric, character, numeric
           marker chr      pos      lod chrom            p
marker_1 marker_1   1 95000000 3.411288  chr1 0.0003878935
marker_2 marker_2   1 95100000 2.909279  chr1 0.0012323127
marker_3 marker_3   1 95200000 3.266316  chr1 0.0005416072
```

All data structures are correct. Issue is purely the missing TxDb package.
