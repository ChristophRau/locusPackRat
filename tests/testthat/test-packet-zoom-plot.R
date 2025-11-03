test_that("generateLocusZoomPlot creates PDF file", {
  # Load test data
  locus_info <- data.frame(
    chr = 1,
    peak_pos = 101.2,
    start_pos = 100,  # Start of locus region in Mb
    end_pos = 102,    # End of locus region in Mb
    trait = "HR",
    drug = "Ctrl",
    max_lod = 7.8
  )

  scan_data <- readRDS(test_path("fixtures/test_scan_data.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data.rds"))
  genes <- data.table::fread(test_path("fixtures/test_genes_mouse.csv"))

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

  output_file <- file.path(
    "/proj/raulab/users/brian/packrat/tests/test_output_2025-01-03",
    "test_locus_zoom.pdf"
  )

  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlaps,
    output_file = output_file,
    assembly = NULL  # Uses default mm39
  )

  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 1000)  # At least 1KB
  expect_equal(result, output_file)
})

test_that("generateLocusZoomPlot handles empty gene list", {
  # Load test data
  locus_info <- data.frame(
    chr = 1,
    peak_pos = 101.2,
    start_pos = 100,
    end_pos = 102,
    trait = "HR",
    drug = "Ctrl",
    max_lod = 7.8
  )

  scan_data <- readRDS(test_path("fixtures/test_scan_data.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data.rds"))
  genes <- data.table::data.table()  # Empty gene list

  top_genes <- data.frame(
    gene = character(0),
    color = character(0)
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

  output_file <- file.path(
    "/proj/raulab/users/brian/packrat/tests/test_output_2025-01-03",
    "test_locus_zoom_empty.pdf"
  )

  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlaps,
    output_file = output_file,
    assembly = NULL
  )

  # Should still create plot without gene track
  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 1000)
})

test_that("generateLocusZoomPlot accepts custom assembly", {
  skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")

  # Load test data
  locus_info <- data.frame(
    chr = 1,
    peak_pos = 101.2,
    start_pos = 100,
    end_pos = 102,
    trait = "HR",
    drug = "Ctrl",
    max_lod = 7.8
  )

  scan_data <- readRDS(test_path("fixtures/test_scan_data.rds"))
  threshold_data <- readRDS(test_path("fixtures/test_threshold_data.rds"))
  genes <- data.table::fread(test_path("fixtures/test_genes_mouse.csv"))

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

  # Create custom assembly object (mm10 instead of default mm39)
  custom_assembly <- plotgardener::assembly(
    Genome = "mm10",
    TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    OrgDb = "org.Mm.eg.db"
  )

  output_file <- file.path(
    "/proj/raulab/users/brian/packrat/tests/test_output_2025-01-03",
    "test_locus_zoom_custom_assembly.pdf"
  )

  result <- generateLocusZoomPlot(
    locus_info = locus_info,
    scan_data = scan_data,
    threshold_data = threshold_data,
    genes_in_locus = genes,
    top_genes_in_locus = top_genes,
    overlapping_loci = overlaps,
    output_file = output_file,
    assembly = custom_assembly
  )

  expect_true(file.exists(output_file))
  expect_gt(file.size(output_file), 1000)
})
