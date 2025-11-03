#!/usr/bin/env Rscript
# Simple test runner using BOTH system R libraries AND miqtl-env libraries

# Use both library paths: system first (has plotgardener), then miqtl-env (has TxDb packages)
.libPaths(c(
  "/nas/longleaf/rhel9/apps/r/4.4.0/lib64/R/library",
  "/nas/longleaf/home/bgural/mambaforge/envs/miqtl-env/lib/R/library"
))

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(plotgardener)
  library(RColorBrewer)
})

# Source the package functions
source("/proj/raulab/users/brian/packrat/R/packet_helpers.R")
source("/proj/raulab/users/brian/packrat/R/packet_core.R")

cat("Testing generateLocusZoomPlot - Basic plot generation\n\n")

# Load test data
locus_info <- data.frame(
  chr = 1,
  peak_pos = 101.2,
  start_pos = 100,  # Clearer than "upper_pos_lod_drop"
  end_pos = 102,    # Clearer than "lower_pos_lod_drop"
  trait = "HR",
  drug = "Ctrl",
  max_lod = 7.8
)

scan_data <- readRDS("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_scan_data.rds")
threshold_data <- readRDS("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_threshold_data.rds")
genes <- fread("/proj/raulab/users/brian/packrat/tests/testthat/fixtures/test_genes_mouse.csv")

cat("Loaded test data:\n")
cat("  - Locus info: chr", locus_info$chr, "peak at", locus_info$peak_pos, "Mb\n")
cat("  - Scan data: ", length(scan_data$HR_Ctrl$LOD), "positions\n")
cat("  - Genes: ", nrow(genes), "genes\n\n")

top_genes <- data.frame(
  gene = c("Gata4", "Nkx2-5"),
  color = "#e34a33"
)

overlaps <- data.frame(
  chrom = "chr1",
  start = 100e6,
  end = 102e6,
  strand = "-",
  trait = "HR",
  drug = "Ctrl",
  traitXdrug = "HR: Ctrl"
)

output_file <- "/proj/raulab/users/brian/packrat/tests/test_output_2025-01-03/test_locus_zoom.pdf"

cat("Scan data structure:\n")
cat("  - HR_Ctrl LOD length:", length(scan_data$HR_Ctrl$LOD), "\n")
cat("  - HR_Ctrl pos$Mb length:", length(scan_data$HR_Ctrl$pos$Mb), "\n")
cat("  - HR_Ctrl chr length:", length(scan_data$HR_Ctrl$chr), "\n")
cat("  - Sample pos$Mb values:", paste(head(scan_data$HR_Ctrl$pos$Mb, 3), collapse=", "), "\n\n")

cat("Generating plot...\n")
result <- tryCatch({
  generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlaps,
    output_file = output_file,
    assembly = NULL  # Uses default mm39
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  traceback()
  return(NULL)
})

if (!is.null(result)) {
  cat("\nTest Results:\n")
  cat("  - File exists:", file.exists(output_file), "\n")
  if (file.exists(output_file)) {
    cat("  - File size:", file.size(output_file), "bytes\n")
    cat("  - Return value matches output_file:", result == output_file, "\n")
  }
  cat("\nTEST PASSED!\n")
} else {
  cat("\nTEST FAILED!\n")
}
