#!/usr/bin/env Rscript

# Test script to check PackRat dependency installation from a minimal R environment
# This simulates what happens when a user installs PackRat without pre-installed dependencies

cat("================================================\n")
cat("PackRat Dependency Installation Test\n")
cat("================================================\n\n")

# Function to check if package is installed
check_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Function to try installing a package
try_install <- function(pkg) {
  tryCatch({
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

# Step 1: Check current R environment
cat("R Environment:\n")
cat("  R Version:", as.character(getRversion()), "\n")
cat("  Library paths:\n")
for (path in .libPaths()) {
  cat("    ", path, "\n")
}
cat("\n")

# Step 2: Read PackRat's DESCRIPTION to get dependencies
desc_file <- "DESCRIPTION"
if (!file.exists(desc_file)) {
  stop("DESCRIPTION file not found. Run this from the PackRat package directory.")
}

desc <- read.dcf(desc_file)

# Parse Imports
imports <- character()
if ("Imports" %in% colnames(desc)) {
  imports_raw <- strsplit(desc[, "Imports"], ",\\s*")[[1]]
  imports <- gsub("\\s*\\([^)]+\\)", "", imports_raw)
  imports <- trimws(imports[imports != ""])
}

# Parse Suggests
suggests <- character()
if ("Suggests" %in% colnames(desc)) {
  suggests_raw <- strsplit(desc[, "Suggests"], ",\\s*")[[1]]
  suggests <- gsub("\\s*\\([^)]+\\)", "", suggests_raw)
  suggests <- trimws(suggests[suggests != ""])
}

cat("PackRat Dependencies:\n")
cat("  Imports (required):", paste(imports, collapse = ", "), "\n")
cat("  Suggests (optional):", paste(suggests, collapse = ", "), "\n\n")

# Step 3: Check which packages are already installed
cat("Checking installed packages...\n")
cat("\nRequired packages (Imports):\n")
missing_imports <- character()
for (pkg in imports) {
  if (check_package(pkg)) {
    cat("  ✓", pkg, "- already installed\n")
  } else {
    cat("  ✗", pkg, "- NOT installed\n")
    missing_imports <- c(missing_imports, pkg)
  }
}

cat("\nOptional packages (Suggests):\n")
missing_suggests <- character()
for (pkg in suggests) {
  if (check_package(pkg)) {
    cat("  ✓", pkg, "- already installed\n")
  } else {
    cat("  ○", pkg, "- not installed (optional)\n")
    missing_suggests <- c(missing_suggests, pkg)
  }
}

# Step 4: Try to install missing required packages
if (length(missing_imports) > 0) {
  cat("\n================================================\n")
  cat("Installing missing required packages...\n")
  cat("================================================\n\n")

  for (pkg in missing_imports) {
    cat("Installing", pkg, "... ")
    if (try_install(pkg)) {
      cat("✓ SUCCESS\n")
    } else {
      cat("✗ FAILED\n")
    }
  }
}

# Step 5: Try to install PackRat itself
cat("\n================================================\n")
cat("Installing PackRat package...\n")
cat("================================================\n\n")

install_result <- tryCatch({
  # Try to install PackRat
  install.packages(".", repos = NULL, type = "source", quiet = FALSE)

  # Try to load it
  library(PackRat)

  TRUE
}, error = function(e) {
  cat("Error during installation:\n")
  cat("  ", e$message, "\n")
  FALSE
})

# Step 6: Test basic functionality if installation succeeded
if (install_result) {
  cat("\n✓ PackRat installed successfully!\n\n")

  cat("Testing basic functionality...\n")

  # Create simple test data
  test_data <- data.frame(
    gene_id = c("G1", "G2", "G3"),
    value = c(10, 20, 30)
  )

  # Test filterGenes
  test_passed <- tryCatch({
    result <- filterGenes(
      inputTable = test_data,
      filters = list(makeFilter("value", ">", 15))
    )
    nrow(result) == 2
  }, error = function(e) {
    FALSE
  })

  if (test_passed) {
    cat("  ✓ filterGenes() works\n")
  } else {
    cat("  ✗ filterGenes() failed\n")
  }

} else {
  cat("\n✗ PackRat installation failed\n")
}

# Step 7: Final summary
cat("\n================================================\n")
cat("Summary\n")
cat("================================================\n\n")

# Re-check all packages
currently_missing_imports <- character()
for (pkg in imports) {
  if (!check_package(pkg)) {
    currently_missing_imports <- c(currently_missing_imports, pkg)
  }
}

if (length(currently_missing_imports) == 0 && install_result) {
  cat("✓ All required dependencies installed\n")
  cat("✓ PackRat installed successfully\n")
  cat("\nPackRat is ready to use!\n")
} else {
  if (length(currently_missing_imports) > 0) {
    cat("✗ Missing required packages:", paste(currently_missing_imports, collapse = ", "), "\n")
  }
  if (!install_result) {
    cat("✗ PackRat installation failed\n")
  }
  cat("\nTo manually install missing packages:\n")
  cat("  install.packages(c(", paste0('"', currently_missing_imports, '"', collapse = ", "), "))\n")
}