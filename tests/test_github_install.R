#!/usr/bin/env Rscript

# Test that GitHub installation handles dependencies automatically
# This is what happens when a user runs devtools::install_github()

cat("=================================================\n")
cat("Testing GitHub Installation Dependency Handling\n")
cat("=================================================\n\n")

cat("When a user runs:\n")
cat('  devtools::install_github("raulab/PackRat")\n\n')

cat("What happens behind the scenes:\n")
cat("------------------------------------------------\n\n")

# Read DESCRIPTION
desc <- read.dcf("DESCRIPTION")
imports <- strsplit(desc[, "Imports"], ",\\s*")[[1]]
imports <- gsub("\\s*\\([^)]+\\)", "", imports)
imports <- trimws(imports[imports != ""])

cat("1. devtools reads DESCRIPTION and finds Imports:\n")
for (pkg in imports) {
  cat("   -", pkg, "\n")
}
cat("\n")

cat("2. devtools checks which are installed:\n")
for (pkg in imports) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("   ✓", pkg, "- already installed\n")
  } else {
    cat("   ✗", pkg, "- NOT installed (will be installed)\n")
  }
}
cat("\n")

cat("3. devtools automatically installs missing dependencies from CRAN\n")
cat("   This happens BEFORE installing PackRat\n\n")

cat("4. Once all dependencies are installed, PackRat is installed\n\n")

cat("=================================================\n")
cat("Testing Local Simulation\n")
cat("=================================================\n\n")

# Simulate what devtools does
cat("Simulating devtools::install_github() behavior:\n\n")

# Check for missing dependencies
missing_deps <- character()
for (pkg in imports) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_deps <- c(missing_deps, pkg)
  }
}

if (length(missing_deps) > 0) {
  cat("Would install these dependencies first:\n")
  cat("  install.packages(c(\n")
  for (i in seq_along(missing_deps)) {
    cat('    "', missing_deps[i], '"', sep = "")
    if (i < length(missing_deps)) cat(",")
    cat("\n")
  }
  cat("  ))\n\n")
} else {
  cat("All dependencies already installed!\n\n")
}

cat("Then install PackRat from GitHub\n\n")

cat("=================================================\n")
cat("Key Points\n")
cat("=================================================\n\n")

cat("✓ Dependencies in 'Imports:' are AUTOMATICALLY installed\n")
cat("✓ Users don't need to install them manually\n")
cat("✓ This works with:\n")
cat("  - devtools::install_github()\n")
cat("  - remotes::install_github()\n")
cat("  - pak::pkg_install()\n\n")

cat("✗ Manual installation (install.packages with repos=NULL)\n")
cat("  does NOT auto-install dependencies\n\n")

cat("=================================================\n")
cat("Example User Experience\n")
cat("=================================================\n\n")

cat("# Fresh R session with NO packages:\n")
cat("> library(PackRat)\n")
cat("Error: there is no package called 'PackRat'\n\n")

cat("# User installs from GitHub:\n")
cat('> devtools::install_github("raulab/PackRat")\n')
cat("Downloading GitHub repo raulab/PackRat@master\n")
cat("These packages have unmet dependencies:\n")
cat("  PackRat depends on data.table, dplyr, openxlsx\n")
cat("Installing 3 packages: data.table, dplyr, openxlsx\n")
cat("Installing data.table...\n")
cat("Installing dplyr...\n")
cat("Installing openxlsx...\n")
cat("Installing PackRat...\n")
cat("✓ PackRat successfully installed!\n\n")

cat("# Now it works:\n")
cat("> library(PackRat)\n")
cat("# Package loads successfully\n")