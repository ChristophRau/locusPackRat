# Integration tests for PackRat
# These tests use larger datasets and test complete workflows
# Should be run on compute nodes via SLURM

library(PackRat)

run_integration_tests <- function() {

  cat("\n=== PackRat Integration Tests ===\n\n")

  # Test 1: Full filtering pipeline with multiple joins
  test_filtering_pipeline <- function() {
    cat("Test 1: Full filtering pipeline...\n")

    # Generate test data
    set.seed(123)
    n_genes <- 1000

    genes <- data.frame(
      gene_id = paste0("ENSMUSG", sprintf("%08d", 1:n_genes)),
      symbol = paste0("Gene", 1:n_genes),
      chr = sample(1:19, n_genes, replace = TRUE),
      biotype = sample(c("protein_coding", "lincRNA", "pseudogene"),
                      n_genes, replace = TRUE, prob = c(0.7, 0.2, 0.1)),
      stringsAsFactors = FALSE
    )

    # Expression data (subset of genes)
    expressed_genes <- sample(genes$gene_id, 600)
    expression <- data.frame(
      gene_id = expressed_genes,
      TPM_heart = rlnorm(600, meanlog = 2, sdlog = 2),
      TPM_liver = rlnorm(600, meanlog = 2, sdlog = 2),
      stringsAsFactors = FALSE
    )

    # GWAS hits (smaller subset)
    gwas_genes <- sample(expressed_genes, 100)
    gwas <- data.frame(
      gene_id = gwas_genes,
      pvalue = 10^-runif(100, 4, 12),
      trait = sample(c("cardiac", "metabolic"), 100, replace = TRUE),
      stringsAsFactors = FALSE
    )

    # Test the pipeline
    result <- genes %>%
      filterGenes(referenceTable = gwas, by = "gene_id", joinType = "semi") %>%
      filterGenes(referenceTable = expression, by = "gene_id", joinType = "left") %>%
      filterGenes(filters = list(
        makeFilter("biotype", "==", "protein_coding"),
        makeFilter("TPM_heart", ">", 10)
      ))

    # Validate results
    stopifnot(nrow(result) > 0)
    stopifnot(nrow(result) < nrow(gwas))  # Should be filtered
    stopifnot(all(result$biotype == "protein_coding"))
    stopifnot(all(result$TPM_heart > 10, na.rm = TRUE))

    cat("  ✓ Passed: ", nrow(genes), "->", nrow(result), "genes\n")
    return(TRUE)
  }

  # Test 2: Complex filtering with custom functions
  test_complex_filtering <- function() {
    cat("Test 2: Complex filtering with custom functions...\n")

    # Generate data
    set.seed(456)
    genes <- data.frame(
      gene_id = paste0("G", 1:500),
      expr_control = runif(500, 0, 100),
      expr_treatment = runif(500, 0, 100),
      pvalue = 10^-runif(500, 0, 10),
      stringsAsFactors = FALSE
    )

    # Calculate fold change
    genes$log2fc <- log2(genes$expr_treatment / (genes$expr_control + 1))

    # Complex filter
    result <- filterGenes(
      genes,
      filters = list(
        # Significant p-value
        makeFilter("pvalue", "<", 0.05),
        # Custom function for fold change
        function(df) abs(df$log2fc) > 1,
        # Expression threshold in at least one condition
        function(df) df$expr_control > 20 | df$expr_treatment > 20
      )
    )

    stopifnot(nrow(result) > 0)
    stopifnot(all(result$pvalue < 0.05))
    stopifnot(all(abs(result$log2fc) > 1))

    cat("  ✓ Passed: ", nrow(genes), "->", nrow(result), "genes\n")
    return(TRUE)
  }

  # Test 3: Performance with larger datasets
  test_performance <- function() {
    cat("Test 3: Performance test with 10,000 genes...\n")

    n_genes <- 10000
    start_time <- Sys.time()

    # Create large dataset
    genes <- data.frame(
      gene_id = paste0("ENSMUSG", sprintf("%08d", 1:n_genes)),
      symbol = paste0("Gene", 1:n_genes),
      value = runif(n_genes, 0, 1000),
      category = sample(LETTERS[1:5], n_genes, replace = TRUE),
      stringsAsFactors = FALSE
    )

    # Create reference table
    ref <- data.frame(
      gene_id = sample(genes$gene_id, 5000),
      annotation = sample(c("important", "candidate"), 5000, replace = TRUE),
      stringsAsFactors = FALSE
    )

    # Filter
    result <- filterGenes(
      genes,
      referenceTable = ref,
      by = "gene_id",
      joinType = "inner",
      filters = list(
        makeFilter("value", ">", 500),
        makeFilter("category", "%in%", c("A", "B", "C"))
      )
    )

    end_time <- Sys.time()
    runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

    stopifnot(nrow(result) > 0)
    stopifnot(runtime < 10)  # Should complete in under 10 seconds

    cat("  ✓ Passed: Processed", n_genes, "genes in", round(runtime, 2), "seconds\n")
    return(TRUE)
  }

  # Test 4: Anti-join functionality
  test_anti_join <- function() {
    cat("Test 4: Anti-join (exclusion) functionality...\n")

    # Create gene list
    all_genes <- data.frame(
      gene_id = paste0("G", 1:100),
      essential = sample(c(TRUE, FALSE), 100, replace = TRUE),
      stringsAsFactors = FALSE
    )

    # Genes to exclude
    exclude_list <- data.frame(
      gene_id = paste0("G", seq(5, 95, by = 10)),  # G5, G15, G25, etc.
      reason = "low_quality",
      stringsAsFactors = FALSE
    )

    # Apply anti-join
    result <- filterGenes(
      all_genes,
      referenceTable = exclude_list,
      by = "gene_id",
      joinType = "anti"
    )

    # Verify exclusion worked
    stopifnot(nrow(result) == nrow(all_genes) - nrow(exclude_list))
    stopifnot(!any(result$gene_id %in% exclude_list$gene_id))

    cat("  ✓ Passed: Excluded", nrow(exclude_list), "genes from", nrow(all_genes), "\n")
    return(TRUE)
  }

  # Test 5: Different column name mapping
  test_column_mapping <- function() {
    cat("Test 5: Column mapping with different names...\n")

    # Mouse genes
    mouse <- data.frame(
      mouse_ensembl = paste0("ENSMUSG", 1:50),
      mouse_symbol = paste0("Mouse", 1:50),
      expression = runif(50, 0, 100),
      stringsAsFactors = FALSE
    )

    # Human orthologs (different column names)
    human <- data.frame(
      human_ensembl = paste0("ENSG", 1:30),
      mouse_id = sample(mouse$mouse_ensembl, 30),
      conservation = runif(30, 70, 100),
      stringsAsFactors = FALSE
    )

    # Map columns
    result <- filterGenes(
      mouse,
      referenceTable = human,
      by = c("mouse_ensembl" = "mouse_id"),
      joinType = "inner",
      filters = list(
        makeFilter("conservation", ">", 85),
        makeFilter("expression", ">", 50)
      )
    )

    stopifnot("conservation" %in% names(result))
    stopifnot(all(result$conservation > 85))
    stopifnot(all(result$expression > 50))

    cat("  ✓ Passed: Mapped", nrow(result), "genes between species\n")
    return(TRUE)
  }

  # Run all tests
  tests <- list(
    test_filtering_pipeline,
    test_complex_filtering,
    test_performance,
    test_anti_join,
    test_column_mapping
  )

  passed <- 0
  failed <- 0

  for (i in seq_along(tests)) {
    tryCatch({
      tests[[i]]()
      passed <- passed + 1
    }, error = function(e) {
      cat("  ✗ FAILED:", e$message, "\n")
      failed <- failed + 1
    })
  }

  cat("\n=== Integration Test Summary ===\n")
  cat("Passed:", passed, "\n")
  cat("Failed:", failed, "\n")

  if (failed > 0) {
    stop("Some integration tests failed")
  }

  return(TRUE)
}

# Run if executed directly
if (!interactive()) {
  run_integration_tests()
}