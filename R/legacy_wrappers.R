#' Legacy Wrapper Functions for Backward Compatibility
#'
#' These functions provide backward compatibility for the old packet_core workflow.
#' They are DEPRECATED and will be removed in a future version.
#'

#' Generate LocusZoom Plot (Legacy)
#'
#' @description
#' \lifecycle{deprecated}
#' This function is deprecated. Please use generateLocusZoomPlot_v2() with the new
#' project-based workflow instead.
#'
#' @inheritParams generateLocusZoomPlot
#' @return Plot file path
#'
#' @export
generateLocusZoomPlot_legacy <- function(
  locus_info,
  scan_data,
  threshold_data,
  genes_in_locus,
  top_genes_in_locus,
  overlapping_loci,
  output_file,
  assembly = NULL,
  plot_params = NULL
) {

  .Deprecated("generateLocusZoomPlot_v2",
              msg = paste0(
                "generateLocusZoomPlot is deprecated.\n",
                "Please use the new project-based workflow:\n",
                "  1. initPackRat() to create a project\n",
                "  2. addRatTable() to add scan data\n",
                "  3. generateLocusZoomPlot_v2() to create plots\n"
              ))

  # Create temporary project
  temp_dir <- tempfile("locuszoom_legacy_")
  dir.create(temp_dir)

  tryCatch({
    # Initialize project with the locus as a region
    region_data <- data.table::data.table(
      chr = locus_info$chr,
      start = if("start_pos" %in% names(locus_info)) locus_info$start_pos * 1e6 else
             min(locus_info$upper_pos_lod_drop, locus_info$lower_pos_lod_drop) * 1e6,
      end = if("end_pos" %in% names(locus_info)) locus_info$end_pos * 1e6 else
           max(locus_info$upper_pos_lod_drop, locus_info$lower_pos_lod_drop) * 1e6,
      region_id = paste0(locus_info$trait, "_", locus_info$drug)
    )

    # Add other columns from locus_info
    for (col in setdiff(names(locus_info), c("chr", "start_pos", "end_pos",
                                              "upper_pos_lod_drop", "lower_pos_lod_drop"))) {
      region_data[[col]] <- locus_info[[col]]
    }

    # Initialize project
    initPackRat(
      data = region_data,
      mode = "region",
      species = "mouse",
      genome = "mm39",
      project_dir = temp_dir,
      force = TRUE
    )

    # Add genes if provided
    if (!is.null(genes_in_locus) && nrow(genes_in_locus) > 0) {
      addRatTable(
        data = genes_in_locus,
        table_name = "genes_in_locus",
        link_type = "region",
        project_dir = temp_dir
      )
    }

    # Add scan data if provided
    if (!is.null(scan_data)) {
      # Convert scan data to table format
      scan_table <- data.table::data.table()

      # Handle the complex scan_data structure
      scan_key <- paste0(locus_info$trait, "_", locus_info$drug)
      if (scan_key %in% names(scan_data)) {
        current_scan <- scan_data[[scan_key]]
        scan_table <- data.table::data.table(
          marker = names(current_scan$LOD),
          chr = as.character(current_scan$chr),
          pos = current_scan$pos$Mb * 1e6,
          lod = current_scan$LOD
        )
      }

      if (nrow(scan_table) > 0) {
        # Add threshold info
        if (!is.null(threshold_data)) {
          threshold_key <- paste0(locus_info$trait, "_", locus_info$drug, "_threshold")
          if (threshold_key %in% names(threshold_data)) {
            scan_table[, threshold := threshold_data[[threshold_key]]]
          }
        }

        addRatTable(
          data = scan_table,
          table_name = "scan_data",
          link_type = "region",
          project_dir = temp_dir
        )
      }
    }

    # Generate plot using new function
    generateLocusZoomPlot_v2(
      project_dir = temp_dir,
      scan_table = "scan_data",
      highlight_genes = if(!is.null(top_genes_in_locus)) top_genes_in_locus$gene else NULL,
      output_file = output_file,
      assembly = assembly,
      plot_params = plot_params
    )

  }, finally = {
    # Clean up temp directory
    unlink(temp_dir, recursive = TRUE)
  })

  invisible(output_file)
}

#' Generate Locus Packet (Legacy)
#'
#' @description
#' \lifecycle{deprecated}
#' This function is deprecated. Please use the new project-based workflow instead.
#'
#' @inheritParams generateLocusPacket
#' @return List of output files created
#'
#' @export
generateLocusPacket_legacy <- function(
  locus_cluster,
  input_path = "data/",
  output_path = "results/qtl_packets/",
  scan_data,
  threshold_data,
  ...
) {

  .Deprecated(msg = paste0(
    "generateLocusPacket is deprecated.\n",
    "Please use the new project-based workflow:\n",
    "  1. initPackRat() to create a project with your regions\n",
    "  2. addRatTable() to add supplementary data\n",
    "  3. generateGeneSheet() for Excel outputs\n",
    "  4. generateLocusZoomPlot_v2() for plots\n"
  ))

  stop("generateLocusPacket has been deprecated. Please use the new workflow.")
}

#' Generate Gene Info Excel (Legacy)
#'
#' @description
#' \lifecycle{deprecated}
#' This function is deprecated. Please use generateGeneSheet() instead.
#'
#' @inheritParams generateGeneInfoExcel
#' @return Excel file path
#'
#' @export
generateGeneInfoExcel_legacy <- function(
  genes_in_locus,
  loci_info,
  merged_gene_info,
  output_file,
  ...
) {

  .Deprecated("generateGeneSheet",
              msg = paste0(
                "generateGeneInfoExcel is deprecated.\n",
                "Please use generateGeneSheet() which provides:\n",
                "  - Dynamic supplementary table inclusion\n",
                "  - Flexible filtering and highlighting\n",
                "  - Multiple output formats (CSV and Excel)\n",
                "  - Better integration with project structure\n"
              ))

  # Create temporary project
  temp_dir <- tempfile("geneinfo_legacy_")
  dir.create(temp_dir)

  tryCatch({
    # Initialize with genes
    if (nrow(genes_in_locus) > 0) {
      initPackRat(
        data = genes_in_locus,
        mode = "gene",
        species = "mouse",
        genome = "mm39",
        project_dir = temp_dir,
        force = TRUE
      )

      # Add merged gene info as supplementary table
      if (!is.null(merged_gene_info)) {
        addRatTable(
          data = merged_gene_info,
          table_name = "gene_annotations",
          link_type = "gene",
          project_dir = temp_dir
        )
      }

      # Generate Excel using new function
      generateGeneSheet(
        format = "excel",
        output_file = output_file,
        include_supplementary = TRUE,
        project_dir = temp_dir
      )
    } else {
      message("No genes found in locus for Excel generation.")
    }

  }, finally = {
    # Clean up
    unlink(temp_dir, recursive = TRUE)
  })

  invisible(output_file)
}