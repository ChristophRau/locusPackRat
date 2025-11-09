#' Generate LocusZoom Plot from locusPackRat Project
#'
#' Creates a LocusZoom-style visualization plot from project data
#'
#' @param region_id Character: ID of the region to plot (NULL for first region or gene-based region)
#' @param project_dir Character: Path to locusPackRat project directory
#' @param scan_table Character: Name of supplementary table containing scan/QTL data
#' @param highlight_genes Character vector: Gene symbols to highlight in the plot
#' @param priority_table Character: Name of supplementary table containing gene priorities/scores
#' @param output_file Character: Path for output PDF (auto-generated if NULL)
#' @param assembly Assembly object for plotgardener (auto-detected from project if NULL)
#' @param plot_params List of plot parameters (uses defaults if NULL)
#'
#' @return Invisible TRUE on success
#'
#' @importFrom data.table fread fwrite setDT
#' @importFrom jsonlite read_json
#' @importFrom plotgardener pageCreate plotManhattan plotGenes plotRanges plotText annoYaxis annoHighlight
#'
#' @export
generateLocusZoomPlot_v2 <- function(
  region_id = NULL,
  project_dir = ".",
  scan_table = NULL,
  highlight_genes = NULL,
  priority_table = NULL,
  output_file = NULL,
  assembly = NULL,
  plot_params = NULL
) {

  # Check required packages
  if (!requireNamespace("plotgardener", quietly = TRUE)) {
    stop("Package 'plotgardener' is required for plotting. Please install it with:\n",
         "  BiocManager::install('plotgardener')")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required. Please install it.")
  }

  # Define color palettes (same as original)
  STRAIN_COLORS <- c(
    "#1B9E77",  # A/J - teal
    "#D95F02",  # C57BL/6J - orange
    "#7570B3",  # 129S1/SvImJ - purple
    "#E7298A",  # NOD/ShiLtJ - magenta
    "#66A61E",  # NZO/HILtJ - green
    "#E6AB02",  # CAST/EiJ - gold
    "#A6761D",  # PWK/PhJ - brown
    "#666666"   # WSB/EiJ - gray
  )

  LOCI_COLORS <- c(
    "#8DD3C7",  # light teal
    "#FFFFB3",  # light yellow
    "#BEBADA",  # light purple
    "#FB8072",  # salmon
    "#80B1D3",  # light blue
    "#FDB462",  # peach
    "#B3DE69",  # light green
    "#FCCDE5"   # light pink
  )

  GENE_HIGHLIGHT_COLOR <- "#e34a33"  # red for highlighted genes
  GENE_BACKGROUND_COLOR <- "#fdbb84"  # light orange for background genes

  # Set default plot parameters
  if (is.null(plot_params)) {
    plot_params <- list(
      page_width = 10.5,
      page_height = 5.5,
      x = 4.25,
      plot_width = 8,
      plot_height = 1,
      plot_y = 0.5
    )
  }

  # --- Load project configuration ---
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  if (!dir.exists(packrat_dir)) {
    stop("No .locusPackRat directory found. Run initPackRat() first.")
  }

  config_file <- file.path(packrat_dir, "config.json")
  if (!file.exists(config_file)) {
    stop("Config file not found. Project may be corrupted.")
  }
  config <- jsonlite::read_json(config_file)

  message(sprintf("Generating LocusZoom plot from %s %s %s project...",
                  config$species, config$genome, config$mode))

  # --- Determine region to plot ---
  if (config$mode == "region") {
    regions_file <- file.path(packrat_dir, "input/regions.csv")
    if (!file.exists(regions_file)) {
      stop("Regions file not found")
    }
    regions <- fread(regions_file)

    if (!is.null(region_id)) {
      locus_info <- regions[region_id == region_id]
      if (nrow(locus_info) == 0) {
        stop("Region ID '", region_id, "' not found")
      }
    } else {
      # Default to first region
      locus_info <- regions[1]
      region_id <- locus_info$region_id
    }
  } else {
    # Gene mode - create synthetic region
    genes_file <- file.path(packrat_dir, "input/genes.csv")
    if (!file.exists(genes_file)) {
      stop("Genes file not found")
    }
    all_genes <- fread(genes_file)

    # Remove genes without coordinates
    all_genes <- all_genes[!is.na(chr) & !is.na(start) & !is.na(end)]

    if (nrow(all_genes) == 0) {
      stop("No genes with coordinates found")
    }

    # Create region encompassing all genes on the most common chromosome
    chr_counts <- table(all_genes$chr)
    main_chr <- names(chr_counts)[which.max(chr_counts)]
    chr_genes <- all_genes[chr == main_chr]

    locus_info <- data.table(
      region_id = paste0("chr", main_chr, "_all_genes"),
      chr = main_chr,
      start = min(chr_genes$start, na.rm = TRUE),
      end = max(chr_genes$end, na.rm = TRUE)
    )
  }

  # --- Load scan data ---
  scan_data <- NULL
  threshold_val <- NULL

  if (!is.null(scan_table)) {
    scan_file <- file.path(packrat_dir, "supplementary", paste0(scan_table, ".csv"))
    if (!file.exists(scan_file)) {
      warning("Scan table '", scan_table, "' not found")
    } else {
      scan_data <- fread(scan_file)
      message("  Loaded scan data from table: ", scan_table)

      # Look for threshold column or use default
      if ("threshold" %in% names(scan_data)) {
        threshold_val <- unique(scan_data$threshold)[1]
      } else if ("significance_threshold" %in% names(scan_data)) {
        threshold_val <- unique(scan_data$significance_threshold)[1]
      } else {
        threshold_val <- 5  # Default LOD threshold
      }
    }
  } else {
    # Auto-detect scan tables
    available_tables <- listPackRatTables(project_dir)
    scan_tables <- available_tables[grepl("scan|qtl|lod", table_name, ignore.case = TRUE)]

    if (nrow(scan_tables) > 0) {
      scan_table <- scan_tables$table_name[1]
      message("  Auto-detected scan table: ", scan_table)
      scan_file <- file.path(packrat_dir, "supplementary", paste0(scan_table, ".csv"))
      scan_data <- fread(scan_file)

      # Look for threshold
      if ("threshold" %in% names(scan_data)) {
        threshold_val <- unique(scan_data$threshold)[1]
      } else {
        threshold_val <- 5
      }
    }
  }

  # --- Load genes in the region ---
  genes_file <- file.path(packrat_dir, "input/genes.csv")
  if (file.exists(genes_file)) {
    all_genes <- fread(genes_file)

    # Filter for genes in the locus region
    genes_in_locus <- all_genes[
      chr == locus_info$chr &
      end >= locus_info$start &
      start <= locus_info$end
    ]

    message("  Found ", nrow(genes_in_locus), " genes in region")
  } else {
    genes_in_locus <- data.table()
  }

  # --- Determine genes to highlight ---
  top_genes_in_locus <- data.table()

  if (!is.null(priority_table)) {
    priority_file <- file.path(packrat_dir, "supplementary", paste0(priority_table, ".csv"))
    if (file.exists(priority_file)) {
      priority_data <- fread(priority_file)

      # Merge with genes in locus
      if ("gene_symbol" %in% names(priority_data) && nrow(genes_in_locus) > 0) {
        genes_in_locus <- merge(genes_in_locus, priority_data,
                               by = "gene_symbol", all.x = TRUE)

        # Find score column
        score_cols <- names(priority_data)[grepl("score|priority|rank", names(priority_data), ignore.case = TRUE)]
        if (length(score_cols) > 0) {
          score_col <- score_cols[1]
          setorderv(genes_in_locus, score_col, order = -1)
          top_genes_in_locus <- genes_in_locus[1:min(10, .N)]
        }
      }
    }
  } else if (!is.null(highlight_genes)) {
    # Use explicitly provided genes
    top_genes_in_locus <- genes_in_locus[gene_symbol %in% highlight_genes]
  } else {
    # Default to top 10 genes by position
    if (nrow(genes_in_locus) > 0) {
      top_genes_in_locus <- genes_in_locus[1:min(10, .N)]
    }
  }

  # Create gene highlight table for plotgardener
  if (nrow(top_genes_in_locus) > 0) {
    gene_highlights <- data.table(
      gene = top_genes_in_locus$gene_symbol,
      color = GENE_HIGHLIGHT_COLOR
    )
  } else {
    gene_highlights <- NULL
  }

  # --- Find overlapping regions (if in region mode) ---
  overlapping_loci <- data.table()
  if (config$mode == "region" && exists("regions")) {
    overlapping_loci <- regions[
      chr == locus_info$chr &
      start <= locus_info$end &
      end >= locus_info$start &
      region_id != locus_info$region_id
    ]

    if (nrow(overlapping_loci) > 0) {
      message("  Found ", nrow(overlapping_loci), " overlapping regions")
    }
  }

  # --- Set up genome assembly ---
  if (is.null(assembly)) {
    if (config$genome == "mm39") {
      assembly <- plotgardener::assembly(
        Genome = "mm39_GRCm39",
        TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",
        OrgDb = "org.Mm.eg.db"
      )
    } else if (config$genome == "hg38") {
      assembly <- plotgardener::assembly(
        Genome = "hg38",
        TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        OrgDb = "org.Hs.eg.db"
      )
    } else if (config$genome == "mm10") {
      assembly <- plotgardener::assembly(
        Genome = "mm10",
        TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
        OrgDb = "org.Mm.eg.db"
      )
    } else if (config$genome == "hg19") {
      assembly <- plotgardener::assembly(
        Genome = "hg19",
        TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
        OrgDb = "org.Hs.eg.db"
      )
    }
  }

  # --- Determine output file ---
  if (is.null(output_file)) {
    output_dir <- file.path(packrat_dir, "output/plots")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    if (!is.null(region_id)) {
      plot_name <- paste0("locuszoom_", region_id, ".pdf")
    } else {
      plot_name <- paste0("locuszoom_chr", locus_info$chr, "_",
                         round(locus_info$start/1e6), "-",
                         round(locus_info$end/1e6), "Mb.pdf")
    }
    output_file <- file.path(output_dir, plot_name)
  }

  # --- Create the plot ---
  message("Creating plot: ", output_file)

  # Calculate plot region (add 0.5 Mb padding)
  plot_start_bp <- max(0, locus_info$start - 5e5)
  plot_end_bp <- locus_info$end + 5e5
  bounds_bp <- c(locus_info$start, locus_info$end)

  # Start plotting
  pdf(output_file, width = plot_params$page_width, height = plot_params$page_height)
  plotgardener::pageCreate(
    width = plot_params$page_width,
    height = plot_params$page_height,
    default.units = "inches",
    showGuides = FALSE
  )

  # Set up genome parameters
  params_genome <- plotgardener::pgParams(
    assembly = assembly,
    chrom = paste0("chr", locus_info$chr),
    chromstart = plot_start_bp,
    chromend = plot_end_bp
  )

  # --- Plot scan data if available ---
  if (!is.null(scan_data) && nrow(scan_data) > 0) {
    # Prepare data for Manhattan plot
    plot_data <- scan_data[chr == locus_info$chr |
                          chromosome == locus_info$chr |
                          chrom == paste0("chr", locus_info$chr)]

    # Standardize column names
    if ("chromosome" %in% names(plot_data)) {
      setnames(plot_data, "chromosome", "chr")
    }
    if ("position" %in% names(plot_data)) {
      setnames(plot_data, "position", "pos")
    }
    if ("lod" %in% names(plot_data)) {
      plot_data[, p := 10^(-lod)]
    } else if ("p_value" %in% names(plot_data)) {
      setnames(plot_data, "p_value", "p")
    }

    # Ensure chromosome format
    if (!grepl("^chr", plot_data$chr[1])) {
      plot_data[, chrom := paste0("chr", chr)]
    } else {
      plot_data[, chrom := chr]
    }

    # Select required columns
    plot_data <- plot_data[, .(chrom, pos, p)]

    # Determine y-axis limits
    ylim <- c(0, max(-log10(plot_data$p), -log10(10^(-threshold_val)), 5, na.rm = TRUE) + 1)

    # Plot Manhattan
    miqtl_plot <- plotgardener::plotManhattan(
      data = plot_data,
      params = params_genome,
      range = ylim,
      trans = "-log10",
      sigVal = 10^(-threshold_val),
      x = plot_params$x,
      y = plot_params$plot_y,
      width = plot_params$plot_width,
      height = plot_params$plot_height,
      just = c("center", "top"),
      fill = "#a6cee3",
      sigCol = "#1f78b4",
      sigLine = TRUE,
      baseline = TRUE,
      default.units = "inches"
    )

    # Add Y axis
    plotgardener::annoYaxis(
      plot = miqtl_plot,
      at = pretty(ylim),
      axisLine = TRUE,
      fontsize = 8
    )

    # Add Y axis label
    plotgardener::plotText(
      label = "LOD Score",
      x = plot_params$x - plot_params$plot_width/2 - 0.3,
      y = plot_params$plot_y + plot_params$plot_height / 2,
      rot = 90,
      fontsize = 8,
      just = "center",
      default.units = "inches"
    )

    # Highlight significant region
    plotgardener::annoHighlight(
      plot = miqtl_plot,
      chrom = paste0("chr", locus_info$chr),
      chromstart = floor(min(bounds_bp)),
      chromend = ceiling(max(bounds_bp)),
      fill = "#fb9a99",
      y = plot_params$plot_y,
      height = plot_params$plot_height,
      just = c("left", "top"),
      default.units = "inches",
      alpha = 0.2,
      params = params_genome
    )

    current_y <- plot_params$plot_y + plot_params$plot_height + 0.2
  } else {
    current_y <- plot_params$plot_y
  }

  # --- Plot overlapping regions if any ---
  if (nrow(overlapping_loci) > 1) {
    # Create color palette function
    loci_palette <- function(n) {
      if (n <= length(LOCI_COLORS)) {
        return(LOCI_COLORS[1:n])
      } else {
        return(grDevices::colorRampPalette(LOCI_COLORS)(n))
      }
    }

    plotgardener::plotRanges(
      data = overlapping_loci,
      params = params_genome,
      fill = plotgardener::colorby("region_id", palette = loci_palette),
      x = plot_params$x,
      y = current_y,
      width = plot_params$plot_width,
      height = 0.5,
      just = c("center", "top"),
      default.units = "inches"
    )

    plotgardener::plotText(
      label = "Overlapping Regions",
      x = plot_params$x - plot_params$plot_width/2 - 0.3,
      y = current_y + 0.25,
      rot = 90,
      fontsize = 8,
      just = "center",
      default.units = "inches"
    )

    current_y <- current_y + 0.7
  }

  # --- Plot genes ---
  # Extract gene names for geneOrder
  gene_order <- if(nrow(top_genes_in_locus) > 0) {
    top_genes_in_locus$gene_symbol
  } else {
    NULL
  }

  # Try to plot genes
  gene_plot <- tryCatch({
    plotgardener::plotGenes(
      params = params_genome,
      x = plot_params$x,
      y = current_y,
      width = plot_params$plot_width,
      height = 1,
      just = c("center", "top"),
      default.units = "inches",
      geneOrder = gene_order,
      fontsize = 6,
      geneHighlights = gene_highlights,
      geneBackground = GENE_BACKGROUND_COLOR
    )
  }, error = function(e) {
    # If highlighting fails, try without it
    message("Note: Gene highlighting failed, plotting without highlights")
    plotgardener::plotGenes(
      params = params_genome,
      x = plot_params$x,
      y = current_y,
      width = plot_params$plot_width,
      height = 1,
      just = c("center", "top"),
      default.units = "inches",
      fontsize = 6
    )
  })

  # Add genome label
  plotgardener::plotGenomeLabel(
    params = params_genome,
    x = plot_params$x,
    y = current_y + 1.1,
    length = plot_params$plot_width,
    just = c("center", "top"),
    default.units = "inches"
  )

  # Add title
  title_text <- if (!is.null(region_id)) {
    paste0("LocusZoom Plot: ", region_id)
  } else {
    paste0("LocusZoom Plot: Chr", locus_info$chr, ":",
           round(locus_info$start/1e6, 1), "-",
           round(locus_info$end/1e6, 1), " Mb")
  }

  plotgardener::plotText(
    label = title_text,
    x = plot_params$page_width / 2,
    y = 0.3,
    fontsize = 12,
    fontface = "bold",
    just = "center",
    default.units = "inches"
  )

  # Add subtitle with project info
  plotgardener::plotText(
    label = paste0(config$species, " ", config$genome, " | ", config$mode, " mode"),
    x = plot_params$page_width / 2,
    y = 0.5,
    fontsize = 10,
    just = "center",
    default.units = "inches"
  )

  dev.off()
  message("Plot saved to: ", output_file)

  invisible(TRUE)
}

#' Generate LocusZoom plots for all regions in a project
#'
#' @param project_dir Character: Path to locusPackRat project directory
#' @param scan_table Character: Name of supplementary table containing scan/QTL data
#' @param ... Additional arguments passed to generateLocusZoomPlot_v2
#'
#' @return Invisible list of output files created
#' @export
generateAllLocusZoomPlots <- function(project_dir = ".", scan_table = NULL, ...) {
  packrat_dir <- file.path(project_dir, ".locusPackRat")
  config <- jsonlite::read_json(file.path(packrat_dir, "config.json"))

  output_files <- character()

  if (config$mode == "region") {
    regions <- fread(file.path(packrat_dir, "input/regions.csv"))

    message("Generating plots for ", nrow(regions), " regions...")
    for (i in seq_len(nrow(regions))) {
      rid <- regions$region_id[i]
      message("\n[", i, "/", nrow(regions), "] Processing ", rid)

      output_file <- generateLocusZoomPlot_v2(
        region_id = rid,
        project_dir = project_dir,
        scan_table = scan_table,
        ...
      )
      output_files <- c(output_files, output_file)
    }
  } else {
    # Single plot for gene mode
    output_file <- generateLocusZoomPlot_v2(
      project_dir = project_dir,
      scan_table = scan_table,
      ...
    )
    output_files <- output_file
  }

  message("\nGenerated ", length(output_files), " plots")
  invisible(output_files)
}