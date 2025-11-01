# Unit tests for filterGenes function
# These are small, fast tests safe to run on login node

library(testthat)
library(PackRat)

test_that("filterGenes handles semi-joins correctly", {
  # Create small test data
  genes <- data.frame(
    gene_id = c("G1", "G2", "G3", "G4", "G5"),
    symbol = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
    stringsAsFactors = FALSE
  )

  reference <- data.frame(
    gene_id = c("G1", "G3", "G5"),
    value = c(100, 200, 300),
    stringsAsFactors = FALSE
  )

  # Test semi-join (keep only matches)
  result <- filterGenes(
    inputTable = genes,
    referenceTable = reference,
    by = "gene_id",
    joinType = "semi"
  )

  expect_equal(nrow(result), 3)
  expect_equal(sort(result$gene_id), c("G1", "G3", "G5"))
  expect_equal(ncol(result), 2)  # No columns added from reference
})

test_that("filterGenes handles anti-joins correctly", {
  genes <- data.frame(
    gene_id = c("G1", "G2", "G3", "G4", "G5"),
    symbol = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
    stringsAsFactors = FALSE
  )

  exclude <- data.frame(
    gene_id = c("G2", "G4"),
    stringsAsFactors = FALSE
  )

  # Test anti-join (exclude matches)
  result <- filterGenes(
    inputTable = genes,
    referenceTable = exclude,
    by = "gene_id",
    joinType = "anti"
  )

  expect_equal(nrow(result), 3)
  expect_equal(sort(result$gene_id), c("G1", "G3", "G5"))
})

test_that("filterGenes handles different column names", {
  genes <- data.frame(
    mouse_id = c("M1", "M2", "M3"),
    mouse_symbol = c("Abc1", "Def2", "Ghi3"),
    stringsAsFactors = FALSE
  )

  human <- data.frame(
    human_gene = c("M1", "M3"),
    human_symbol = c("ABC1", "GHI3"),
    stringsAsFactors = FALSE
  )

  # Test with column name mapping
  result <- filterGenes(
    inputTable = genes,
    referenceTable = human,
    by = c("mouse_id" = "human_gene"),
    joinType = "inner"
  )

  expect_equal(nrow(result), 2)
  expect_true("human_symbol" %in% names(result))
  expect_equal(result$mouse_symbol, c("Abc1", "Ghi3"))
})

test_that("filterGenes applies custom filters correctly", {
  genes <- data.frame(
    gene_id = paste0("G", 1:10),
    expression = c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95),
    biotype = rep(c("protein_coding", "lincRNA"), 5),
    stringsAsFactors = FALSE
  )

  # Test with custom filters
  result <- filterGenes(
    inputTable = genes,
    filters = list(
      makeFilter("expression", ">", 30),
      makeFilter("biotype", "==", "protein_coding")
    )
  )

  expect_equal(nrow(result), 3)  # G4, G6, G8, G10 have >30, but only G4, G6, G8 are protein_coding
  expect_true(all(result$expression > 30))
  expect_true(all(result$biotype == "protein_coding"))
})

test_that("filterGenes handles function filters", {
  genes <- data.frame(
    gene_id = paste0("G", 1:10),
    expr_heart = runif(10, 0, 100),
    expr_liver = runif(10, 0, 100),
    stringsAsFactors = FALSE
  )

  # Test with custom function
  result <- filterGenes(
    inputTable = genes,
    filters = list(
      function(df) df$expr_heart > df$expr_liver
    )
  )

  expect_true(all(result$expr_heart > result$expr_liver))
})

test_that("makeFilter creates correct specifications", {
  # Test basic filter
  f1 <- makeFilter("column1", "==", "value1")
  expect_equal(f1$column, "column1")
  expect_equal(f1$condition, "==")
  expect_equal(f1$value, "value1")
  expect_equal(f1$na.rm, FALSE)

  # Test with na.rm
  f2 <- makeFilter("column2", ">", 100, na.rm = TRUE)
  expect_equal(f2$na.rm, TRUE)

  # Test between condition
  expect_error(makeFilter("col", "between", 5))  # Should fail - needs 2 values
  f3 <- makeFilter("col", "between", c(5, 10))
  expect_equal(length(f3$value), 2)

  # Test is_na doesn't need value
  f4 <- makeFilter("col", "is_na")
  expect_null(f4$value)
})

test_that("filterGenes handles empty results gracefully", {
  genes <- data.frame(
    gene_id = c("G1", "G2", "G3"),
    value = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  # Filter that matches nothing
  result <- filterGenes(
    inputTable = genes,
    filters = list(makeFilter("value", ">", 100))
  )

  expect_equal(nrow(result), 0)
  expect_equal(names(result), names(genes))
})

test_that("filterGenes preserves data.table class", {
  skip_if_not_installed("data.table")

  dt <- data.table::data.table(
    gene_id = c("G1", "G2"),
    value = c(10, 20)
  )

  result <- filterGenes(
    inputTable = dt,
    filters = list(makeFilter("value", ">", 5))
  )

  expect_true(inherits(result, "data.table"))
})