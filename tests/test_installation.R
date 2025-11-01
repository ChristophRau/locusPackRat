#!/usr/bin/env Rscript

# Test script to verify PackRat can be installed with all dependencies
# This ensures that dependency management is working correctly

cat("=====================================\n")
cat("PackRat Installation Test\n")
cat("=====================================\n\n")

# Function to test package installation
test_package_installation <- function() {

  # Step 1: Check R version
  r_version <- getRversion()
  cat("R Version:", as.character(r_version), "\n")

  if (r_version < "4.3.0") {
    stop("PackRat requires R >= 4.3.0. Current version: ", r_version)
  }

  # Step 2: Test installation (doesn't actually install, just checks)
  cat("\nChecking if PackRat can be installed...\n")

  # Get DESCRIPTION file
  desc <- read.dcf("DESCRIPTION")

  # Check Imports
  imports <- strsplit(desc[, "Imports"], ",\\s*")[[1]]
  imports <- gsub("\\s*\\([^)]+\\)", "", imports)  # Remove version specs
  imports <- trimws(imports)

  cat("\nRequired packages (Imports):\n")
  missing_imports <- c()
  for (pkg in imports) {
    if (pkg == "") next
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat("  ✓", pkg, "\n")
    } else {
      cat("  ✗", pkg, "(not installed)\n")
      missing_imports <- c(missing_imports, pkg)
    }
  }

  # Check Suggests
  if ("Suggests" %in% colnames(desc)) {
    suggests <- strsplit(desc[, "Suggests"], ",\\s*")[[1]]
    suggests <- gsub("\\s*\\([^)]+\\)", "", suggests)
    suggests <- trimws(suggests)

    cat("\nOptional packages (Suggests):\n")
    missing_suggests <- c()
    for (pkg in suggests) {
      if (pkg == "") next
      if (requireNamespace(pkg, quietly = TRUE)) {
        cat("  ✓", pkg, "\n")
      } else {
        cat("  ○", pkg, "(not installed - optional)\n")
        missing_suggests <- c(missing_suggests, pkg)
      }
    }
  }

  # Step 3: Try to install PackRat (local install)
  if (length(missing_imports) > 0) {
    cat("\n⚠ Missing required packages:\n")
    cat("  ", paste(missing_imports, collapse = ", "), "\n")
    cat("\nTo install missing packages, run:\n")
    cat("  install.packages(c(",
        paste0('"', missing_imports, '"', collapse = ", "),
        "))\n")
    return(FALSE)
  }

  cat("\n✓ All required dependencies are available\n")

  # Step 4: Actually try to load/install the package
  cat("\nAttempting to install PackRat locally...\n")

  tryCatch({
    # Try to install
    install.packages(".", repos = NULL, type = "source", quiet = TRUE)

    # Try to load
    library(PackRat)

    cat("✓ PackRat installed successfully!\n")

    # Test basic functionality
    cat("\nTesting basic functionality...\n")

    # Test 1: Create simple data
    test_data <- data.frame(
      gene_id = c("G1", "G2", "G3"),
      value = c(10, 20, 30)
    )

    # Test 2: Use filterGenes
    result <- filterGenes(
      inputTable = test_data,
      filters = list(makeFilter("value", ">", 15))
    )

    if (nrow(result) == 2) {
      cat("  ✓ filterGenes works correctly\n")
    } else {
      cat("  ✗ filterGenes test failed\n")
      return(FALSE)
    }

    # Clean up
    detach("package:PackRat", unload = TRUE)

    return(TRUE)

  }, error = function(e) {
    cat("✗ Installation failed:\n")
    cat("  ", e$message, "\n")
    return(FALSE)
  })
}

# Run the test
success <- test_package_installation()

if (success) {
  cat("\n=====================================\n")
  cat("✓ All tests passed!\n")
  cat("=====================================\n")
  quit(status = 0)
} else {
  cat("\n=====================================\n")
  cat("✗ Tests failed. See errors above.\n")
  cat("=====================================\n")
  quit(status = 1)
}