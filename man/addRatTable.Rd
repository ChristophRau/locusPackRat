\name{addRatTable}
\alias{addRatTable}
\title{Add Supplementary Data Table to locusPackRat Project}
\description{
  Adds a supplementary data table to an existing locusPackRat project. The table
  can be linked to either genes or regions and will be automatically integrated
  when generating outputs.
}
\usage{
addRatTable(data, table_name, link_type = c("gene", "region"),
            link_by = NULL, project_dir = ".", verbose = TRUE)
}
\arguments{
  \item{data}{A data frame containing supplementary data to add to the project.}
  \item{table_name}{Name for this supplementary table (will be used as filename
    and reference in outputs).}
  \item{link_type}{Either "gene" to link by gene symbols or "region" to link by
    genomic coordinates.}
  \item{link_by}{Column name(s) to use for linking. If NULL, automatically detects
    based on column names. For genes, looks for gene_symbol, gene_name, or symbol.
    For regions, requires chr, start, and end columns.}
  \item{project_dir}{Path to the project directory. Defaults to current directory.}
  \item{verbose}{If TRUE, prints progress messages. Default is TRUE.}
}
\value{
  Invisibly returns TRUE on success. The table is saved in the
  .locusPackRat/supplementary/ directory and the config.json is updated.
}
\examples{
# Add expression data to a gene project
expression_data <- data.frame(
  gene_symbol = c("BRCA1", "TP53", "EGFR"),
  tissue = c("brain", "liver", "kidney"),
  expression_level = c(100, 250, 50)
)
addRatTable(expression_data, table_name = "expression",
            link_type = "gene", project_dir = "my_project")

# Add QTL scan results to a region project
scan_results <- data.frame(
  chr = rep(5, 100),
  pos = seq(10000000, 15000000, length.out = 100),
  lod = rnorm(100, mean = 3, sd = 1)
)
addRatTable(scan_results, table_name = "qtl_scans",
            link_type = "region", project_dir = "qtl_project")
}
\seealso{
  \code{\link{initPackRat}}, \code{\link{listPackRatTables}},
  \code{\link{generateGeneSheet}}
}
\author{
  Brian Gural, Todd Kimball, Anh Luu, Christoph D. Rau
}