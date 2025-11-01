#' Create LocusZoom-style plot for genomic region
#' 
#' @description
#' Generate publication-quality plots of genomic loci showing genes,
#' variants, and association statistics using plotgardener.
#' 
#' @param locus_chr Chromosome of the locus
#' @param locus_start Start position of the locus
#' @param locus_end End position of the locus
#' @param genes_df Data frame with gene annotations
#' @param variants_df Optional data frame with variant positions
#' @param association_df Optional data frame with association statistics
#' @param expression_df Optional data frame with expression values
#' @param highlight_genes Character vector of genes to highlight
#' @param output_file Path to save plot (PNG or PDF)
#' @param plot_height Height of plot in inches
#' @param plot_width Width of plot in inches
#' @param genome_build Genome build (e.g., "mm10", "hg38")
#' 
#' @return Invisible NULL (creates plot as side effect)
#' 
#' @examples
#' \dontrun{
#' genes <- data.frame(
#'   gene = c("Abcb10", "Acsl5"),
#'   chr = c(8, 8),
#'   start = c(28500000, 28700000),
#'   end = c(28550000, 28750000),
#'   strand = c("+", "-")
#' )
#' plot_locus(chr = 8, start = 28000000, end = 29000000,
#'           genes_df = genes, output_file = "locus.png")
#' }
#' 
#' @export
plot_locus <- function(locus_chr,
                      locus_start,
                      locus_end,
                      genes_df = NULL,
                      variants_df = NULL,
                      association_df = NULL,
                      expression_df = NULL,
                      highlight_genes = NULL,
                      output_file = NULL,
                      plot_height = 8,
                      plot_width = 12,
                      genome_build = "mm10") {
  
  # Check for plotgardener
  if (!requireNamespace("plotgardener", quietly = TRUE)) {
    message("plotgardener not available. Falling back to base R plotting.")
    return(plot_locus_base(locus_chr, locus_start, locus_end, genes_df,
                          variants_df, association_df, expression_df,
                          highlight_genes, output_file, plot_height, plot_width))
  }
  
  # Set up page
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = plot_width, height = plot_height)
    } else {
      png(output_file, width = plot_width, height = plot_height, 
          units = "in", res = 300)
    }
  }
  
  plotgardener::pageCreate(width = plot_width, height = plot_height, 
                          default.units = "inches")
  
  # Define layout
  y_pos <- 0.5
  panel_height <- 1.5
  
  # Panel 1: Association plot (if provided)
  if (!is.null(association_df)) {
    # Manhattan-style plot
    assoc_plot <- plotgardener::plotManhattan(
      data = association_df,
      chrom = locus_chr,
      chromstart = locus_start,
      chromend = locus_end,
      assembly = genome_build,
      x = 0.5, y = y_pos,
      width = plot_width - 1, height = panel_height,
      just = c("left", "top")
    )
    
    # Add significance line
    abline(h = -log10(5e-8), lty = 2, col = "red")
    
    y_pos <- y_pos + panel_height + 0.25
  }
  
  # Panel 2: Gene track
  if (!is.null(genes_df)) {
    # Prepare gene data for plotgardener
    gene_colors <- rep("darkblue", nrow(genes_df))
    if (!is.null(highlight_genes)) {
      gene_colors[genes_df$gene %in% highlight_genes] <- "red"
    }
    
    gene_plot <- plotgardener::plotGenes(
      chrom = paste0("chr", locus_chr),
      chromstart = locus_start,
      chromend = locus_end,
      assembly = genome_build,
      x = 0.5, y = y_pos,
      width = plot_width - 1, height = panel_height,
      just = c("left", "top"),
      geneOrder = "strand",
      fill = gene_colors
    )
    
    # Add gene labels for highlighted genes
    if (!is.null(highlight_genes)) {
      highlighted <- genes_df[genes_df$gene %in% highlight_genes, ]
      for (i in 1:nrow(highlighted)) {
        plotgardener::plotText(
          label = highlighted$gene[i],
          x = (highlighted$start[i] + highlighted$end[i]) / 2,
          y = y_pos + panel_height/2,
          just = "center",
          fontsize = 8,
          fontface = "bold"
        )
      }
    }
    
    y_pos <- y_pos + panel_height + 0.25
  }
  
  # Panel 3: Expression heatmap (if provided)
  if (!is.null(expression_df)) {
    # Create expression matrix
    expr_matrix <- as.matrix(expression_df[, -1])  # Assuming first column is gene names
    rownames(expr_matrix) <- expression_df[, 1]
    
    # Filter to genes in region
    genes_in_region <- genes_df$gene[genes_df$chr == locus_chr & 
                                     genes_df$start >= locus_start & 
                                     genes_df$end <= locus_end]
    expr_matrix <- expr_matrix[rownames(expr_matrix) %in% genes_in_region, , drop = FALSE]
    
    if (nrow(expr_matrix) > 0) {
      # Create heatmap
      plotgardener::plotHeatmap(
        data = expr_matrix,
        x = 0.5, y = y_pos,
        width = plot_width - 1, height = panel_height,
        just = c("left", "top"),
        palette = colorRampPalette(c("white", "red"))(100)
      )
      
      y_pos <- y_pos + panel_height + 0.25
    }
  }
  
  # Panel 4: Variant track (if provided)
  if (!is.null(variants_df)) {
    # Plot variants as vertical lines
    variants_in_region <- variants_df[variants_df$chr == locus_chr &
                                     variants_df$pos >= locus_start &
                                     variants_df$pos <= locus_end, ]
    
    if (nrow(variants_in_region) > 0) {
      plot(0, type = "n", 
           xlim = c(locus_start, locus_end),
           ylim = c(0, 1),
           xlab = sprintf("Chromosome %s position (Mb)", locus_chr),
           ylab = "Variants",
           main = "Coding Variants")
      
      # Color by variant type
      var_colors <- ifelse(grepl("missense", variants_in_region$type), "orange",
                          ifelse(grepl("nonsense|stop", variants_in_region$type), "red",
                                "gray"))
      
      segments(x0 = variants_in_region$pos, y0 = 0,
              x1 = variants_in_region$pos, y1 = 0.8,
              col = var_colors, lwd = 2)
      
      legend("topright", 
             legend = c("Missense", "Nonsense", "Other"),
             col = c("orange", "red", "gray"),
             lwd = 2, bty = "n")
    }
  }
  
  # Add title
  plotgardener::plotText(
    label = sprintf("Locus: Chr%s:%s-%s Mb",
                   locus_chr,
                   format(locus_start/1e6, digits = 3),
                   format(locus_end/1e6, digits = 3)),
    x = plot_width/2, y = 0.25,
    just = "center",
    fontsize = 14,
    fontface = "bold"
  )
  
  # Close device if saving to file
  if (!is.null(output_file)) {
    dev.off()
    message(sprintf("Plot saved to %s", output_file))
  }
  
  invisible(NULL)
}

#' Create locus plot using base R (fallback)
#' 
#' @description
#' Fallback plotting function using base R graphics when plotgardener
#' is not available.
#' 
#' @inheritParams plot_locus
#' @keywords internal
plot_locus_base <- function(locus_chr,
                           locus_start,
                           locus_end,
                           genes_df = NULL,
                           variants_df = NULL,
                           association_df = NULL,
                           expression_df = NULL,
                           highlight_genes = NULL,
                           output_file = NULL,
                           plot_height = 8,
                           plot_width = 12) {
  
  # Set up device
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = plot_width, height = plot_height)
    } else {
      png(output_file, width = plot_width * 100, height = plot_height * 100)
    }
  }
  
  # Determine number of panels
  n_panels <- sum(!is.null(association_df), !is.null(genes_df), 
                 !is.null(variants_df), !is.null(expression_df))
  
  if (n_panels == 0) {
    stop("No data provided to plot")
  }
  
  # Set up multi-panel layout
  par(mfrow = c(n_panels, 1), mar = c(4, 4, 2, 2))
  
  # Panel 1: Association plot
  if (!is.null(association_df)) {
    plot(association_df$pos, -log10(association_df$p_value),
         xlim = c(locus_start, locus_end),
         xlab = sprintf("Chromosome %s position", locus_chr),
         ylab = "-log10(P)",
         main = "Association Statistics",
         pch = 19, col = "darkblue")
    abline(h = -log10(5e-8), lty = 2, col = "red")
  }
  
  # Panel 2: Gene track
  if (!is.null(genes_df)) {
    genes_in_region <- genes_df[genes_df$chr == locus_chr &
                               genes_df$start <= locus_end &
                               genes_df$end >= locus_start, ]
    
    plot(0, type = "n",
         xlim = c(locus_start, locus_end),
         ylim = c(0, nrow(genes_in_region) + 1),
         xlab = sprintf("Chromosome %s position", locus_chr),
         ylab = "",
         yaxt = "n",
         main = "Gene Annotations")
    
    # Plot genes as rectangles
    for (i in 1:nrow(genes_in_region)) {
      gene_color <- ifelse(genes_in_region$gene[i] %in% highlight_genes, 
                          "red", "darkblue")
      rect(xleft = genes_in_region$start[i],
           xright = genes_in_region$end[i],
           ybottom = i - 0.4,
           ytop = i + 0.4,
           col = gene_color,
           border = "black")
      
      # Add gene name
      text(x = (genes_in_region$start[i] + genes_in_region$end[i])/2,
           y = i,
           labels = genes_in_region$gene[i],
           cex = 0.7,
           col = "white",
           font = 2)
    }
  }
  
  # Panel 3: Variant track
  if (!is.null(variants_df)) {
    variants_in_region <- variants_df[variants_df$chr == locus_chr &
                                     variants_df$pos >= locus_start &
                                     variants_df$pos <= locus_end, ]
    
    if (nrow(variants_in_region) > 0) {
      # Create variant density plot
      hist(variants_in_region$pos,
           breaks = seq(locus_start, locus_end, length.out = 50),
           xlim = c(locus_start, locus_end),
           xlab = sprintf("Chromosome %s position", locus_chr),
           ylab = "Variant count",
           main = "Variant Distribution",
           col = "gray80",
           border = "gray60")
    }
  }
  
  # Panel 4: Expression heatmap
  if (!is.null(expression_df)) {
    # Simple representation - bar plot of mean expression
    genes_in_region <- genes_df$gene[genes_df$chr == locus_chr &
                                     genes_df$start <= locus_end &
                                     genes_df$end >= locus_start]
    
    expr_subset <- expression_df[expression_df$gene %in% genes_in_region, ]
    
    if (nrow(expr_subset) > 0) {
      # Calculate mean expression
      expr_means <- rowMeans(expr_subset[, -1], na.rm = TRUE)
      names(expr_means) <- expr_subset$gene
      
      barplot(expr_means,
              main = "Gene Expression",
              ylab = "Mean Expression",
              las = 2,
              col = ifelse(names(expr_means) %in% highlight_genes, "red", "darkblue"))
    }
  }
  
  # Close device
  if (!is.null(output_file)) {
    dev.off()
    message(sprintf("Plot saved to %s", output_file))
  }
  
  invisible(NULL)
}

#' Create summary visualization dashboard
#' 
#' @description
#' Generate a multi-panel dashboard summarizing locus characteristics.
#' 
#' @param locus_data Locus information
#' @param gene_data Gene annotations and data
#' @param output_file Path to save plot
#' 
#' @return Invisible NULL
#' 
#' @export
plot_locus_dashboard <- function(locus_data,
                                gene_data,
                                output_file = NULL) {
  
  # Set up device
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 12, height = 8)
    } else {
      png(output_file, width = 1200, height = 800)
    }
  }
  
  # Create 2x2 layout
  par(mfrow = c(2, 2), mar = c(5, 4, 3, 2))
  
  # Panel 1: Gene biotype distribution
  biotype_table <- table(gene_data$biotype)
  pie(biotype_table,
      main = "Gene Biotypes",
      col = rainbow(length(biotype_table)),
      labels = sprintf("%s\n(%d)", names(biotype_table), biotype_table))
  
  # Panel 2: Expression distribution
  if (any(grepl("express|CPM|TPM", names(gene_data)))) {
    expr_col <- names(gene_data)[grepl("express|CPM|TPM", names(gene_data))][1]
    hist(log10(gene_data[[expr_col]] + 1),
         main = "Expression Distribution",
         xlab = "log10(Expression + 1)",
         ylab = "Number of genes",
         col = "skyblue",
         border = "darkblue")
  } else {
    plot.new()
    text(0.5, 0.5, "No expression data", cex = 1.5)
  }
  
  # Panel 3: Evidence summary
  evidence_summary <- data.frame(
    Evidence = c("Human ortholog", "Expression > 1", "Has eQTL", "Has variants"),
    Count = c(
      sum(!is.na(gene_data$human_ortholog), na.rm = TRUE),
      sum(gene_data$expression > 1, na.rm = TRUE),
      sum(gene_data$has_eqtl == TRUE, na.rm = TRUE),
      sum(gene_data$n_variants > 0, na.rm = TRUE)
    )
  )
  
  barplot(evidence_summary$Count,
          names.arg = evidence_summary$Evidence,
          main = "Evidence Summary",
          ylab = "Number of genes",
          col = "coral",
          las = 2)
  
  # Panel 4: Priority score distribution
  if ("priority_score" %in% names(gene_data)) {
    hist(gene_data$priority_score,
         main = "Priority Score Distribution",
         xlab = "Priority Score",
         ylab = "Number of genes",
         col = "lightgreen",
         border = "darkgreen")
  } else {
    # Calculate simple priority score
    gene_data$priority_score <- 0
    if ("biotype" %in% names(gene_data)) {
      gene_data$priority_score <- gene_data$priority_score + 
        ifelse(gene_data$biotype == "protein_coding", 1, 0)
    }
    if ("human_ortholog" %in% names(gene_data)) {
      gene_data$priority_score <- gene_data$priority_score + 
        ifelse(!is.na(gene_data$human_ortholog), 1, 0)
    }
    
    hist(gene_data$priority_score,
         main = "Priority Score Distribution",
         xlab = "Priority Score",
         ylab = "Number of genes",
         col = "lightgreen",
         border = "darkgreen")
  }
  
  # Add overall title
  mtext(sprintf("Locus Dashboard: %s", locus_data$locus_id[1]),
        outer = TRUE, cex = 1.5, line = -2)
  
  # Close device
  if (!is.null(output_file)) {
    dev.off()
    message(sprintf("Dashboard saved to %s", output_file))
  }
  
  invisible(NULL)
}