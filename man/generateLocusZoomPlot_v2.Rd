\name{generateLocusZoomPlot_v2}
\alias{generateLocusZoomPlot_v2}
\title{Generate LocusZoom-style Plot from locusPackRat Project}
\description{
  Creates a LocusZoom-style visualization for a genomic region, displaying
  association statistics (e.g., QTL LOD scores) along with gene annotations.
  Works with the new project-based workflow and automatically integrates
  supplementary data.
}
\usage{
generateLocusZoomPlot_v2(region_id = NULL, project_dir = ".",
                         scan_table = NULL, highlight_genes = NULL,
                         priority_table = NULL, output_file = NULL,
                         assembly = NULL, plot_params = NULL)
}
\arguments{
  \item{region_id}{ID of the region to plot (for region-based projects). If NULL,
    uses the first region.}
  \item{project_dir}{Path to the project directory. Defaults to current directory.}
  \item{scan_table}{Name of the supplementary table containing scan/association
    data. If NULL, attempts to auto-detect.}
  \item{highlight_genes}{Vector of gene symbols to highlight on the plot.}
  \item{priority_table}{Name of supplementary table with gene priorities/scores
    for enhanced visualization.}
  \item{output_file}{Path for the output plot file. If NULL, auto-generates in
    output/plots/ directory.}
  \item{assembly}{Genome assembly for gene annotations. If NULL, uses project
    settings.}
  \item{plot_params}{List of plotting parameters to customize the visualization.}
}
\value{
  Invisibly returns the path to the created plot file. Creates a multi-panel
  plot with:
  \itemize{
    \item Top panel: Association statistics (e.g., LOD scores)
    \item Bottom panel: Gene annotations with transcripts
    \item Optional highlighting of specific genes
    \item Threshold lines if available
  }
}
\note{
  This function requires the plotgardener package from Bioconductor. Install it
  using: BiocManager::install("plotgardener")
}
\examples{
\dontrun{
# Generate plot for first region with auto-detected scan data
generateLocusZoomPlot_v2(project_dir = "qtl_project")

# Generate plot with specific scan table and highlighted genes
generateLocusZoomPlot_v2(
  region_id = "chr5_qtl",
  scan_table = "qtl_scans",
  highlight_genes = c("Gene1", "Gene2"),
  output_file = "results/chr5_locus.pdf",
  project_dir = "qtl_project"
)

# Generate plot with custom parameters
generateLocusZoomPlot_v2(
  project_dir = "qtl_project",
  plot_params = list(
    height = 8,
    width = 10,
    point_size = 0.5,
    gene_label_size = 8
  )
)
}
}
\seealso{
  \code{\link{initPackRat}}, \code{\link{addRatTable}}
}
\author{
  Brian Gural, Todd Kimball, Anh Luu, Christoph D. Rau
}