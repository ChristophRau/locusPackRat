# Next Steps for Testing - January 3, 2025

## Current Status

### Completed ✅
1. Test fixtures created with correct data structures
2. Test file `test-packet-zoom-plot.R` written (3 tests)
3. Library paths configured to access both system (plotgardener) and miqtl-env (TxDb packages)
4. All test data verified as correct format

### Blocked 🚫

**plotManhattan Issue**: plotgardener v1.12.0 throws error "'pos' column must be an integer or numeric" even though the column IS integer.

This appears to be an environment/version issue. The column is verified as:
- `class(pos)` = "integer"
- `is.numeric(pos)` = TRUE
- `is.integer(pos)` = TRUE

But plotManhattan still rejects it.

## Recommended Next Actions

### Option 1: Test in Original Working Environment ⭐ RECOMMENDED
Run the tests in the same environment where the original `21_makeLociPackets.R` script successfully runs. This likely has all the correct package versions and dependencies.

**Command**:
```bash
# From the original scripts directory
cd /proj/raulab/users/brian/cc_gwas/scripts/packets
# Check what R environment is used
which R
R --version
```

Then replicate that environment for our tests.

### Option 2: Skip Plot Tests, Focus on Excel Tests
The Excel generation tests don't require plotgardener and should run immediately.

**Action**: Implement `tests/testthat/test-packet-excel.R` (5 tests defined in test_design.md)

Required fixtures (need to create):
- `test_orthology.csv`
- `test_associations.csv`
- `test_mouse_pheno.csv`
- `test_merged_gene_info.csv`
- `test_rna_expression.csv`

### Option 3: Debug plotgardener Further
Investigate plotgardener source code or try different data.frame construction methods. This is lower priority since the original script works.

## Files Ready for Testing

### Locus Zoom Plot Tests (blocked)
- ✅ `tests/testthat/fixtures/test_genes_mouse.csv`
- ✅ `tests/testthat/fixtures/test_scan_data.rds`
- ✅ `tests/testthat/fixtures/test_threshold_data.rds`
- ✅ `tests/testthat/test-packet-zoom-plot.R`
- ✅ `tests/run_zoom_test_system.R`

### Excel Tests (ready to implement)
- ⏳ Need to create additional fixtures
- ⏳ Need to implement `tests/testthat/test-packet-excel.R`

### Integration Tests (ready to implement)
- ⏳ Need `tests/testthat/test-packet-core.R`
- ⏳ Need `tests/testthat/helper-setup.R`

## Priority Queue

1. **PRIORITY 1**: Identify and use the working R environment from original scripts
2. **PRIORITY 2**: Create Excel test fixtures and implement Excel tests
3. **PRIORITY 3**: Package rename (genePackRat → locusPackRat)
4. **PRIORITY 4**: Integration tests
5. **PRIORITY 5**: R CMD check

## plotgardener Debug Log

**Tested configurations**:
- ✅ System library only: TxDb missing
- ✅ miqtl-env first: OpenSSL conflicts
- ✅ System first + miqtl-env: TxDb found, but plotManhattan errors
- ✅ Explicit integer conversion: Still errors
- ✅ Without xfield/yfield parameters: Still errors

**Verified**:
- Data frame structure matches documentation requirements
- Column types are correct (integer, numeric, character as specified)
- Assembly object creates successfully
- TxDb packages load correctly

**Conclusion**: Issue is likely plotgardener version mismatch or internal bug, not our data/code.
