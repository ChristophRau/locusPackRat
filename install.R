#!/usr/bin/env Rscript

# Simple one-line installation script for PackRat with automatic dependencies
# Users can source this file to install everything automatically

# Check and install dependencies
deps <- c("data.table", "dplyr", "openxlsx")
new_deps <- deps[!sapply(deps, requireNamespace, quietly = TRUE)]
if (length(new_deps) > 0) {
  cat("Installing dependencies:", paste(new_deps, collapse = ", "), "\n")
  install.packages(new_deps, repos = "https://cloud.r-project.org")
}

# Install PackRat
cat("Installing PackRat...\n")
install.packages(".", repos = NULL, type = "source")

# Verify
if (requireNamespace("PackRat", quietly = TRUE)) {
  cat("\n✓ PackRat installed successfully!\n")
  cat("Load with: library(PackRat)\n")
} else {
  stop("Installation failed")
}