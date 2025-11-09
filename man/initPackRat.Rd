\name{initPackRat}
\alias{initPackRat}
\title{Initialize a locusPackRat Project}
\description{
  Creates a new project directory structure for managing genomic analyses with
  persistent data storage. The project can be initialized with either genes or
  genomic regions and automatically loads coordinate and orthology information.
}
\usage{
initPackRat(data, mode = c("gene", "region"), species = c("human", "mouse"),
            genome = c("hg38", "hg19", "mm39", "mm10"),
            project_dir = ".", force = FALSE, verbose = TRUE)
}
\arguments{
  \item{data}{A data frame containing either gene symbols or genomic regions.
    For gene mode, must have a column with gene symbols. For region mode, must
    have chr, start, and end columns.}
  \item{mode}{Either "gene" for gene-based projects or "region" for genomic
    region-based projects (e.g., QTL intervals).}
  \item{species}{Either "human" or "mouse" to specify the primary species.}
  \item{genome}{The reference genome version: "hg38" or "hg19" for human,
    "mm39" or "mm10" for mouse.}
  \item{project_dir}{Path where the project should be created. Defaults to
    current directory.}
  \item{force}{If TRUE, overwrites existing project. Default is FALSE.}
  \item{verbose}{If TRUE, prints progress messages. Default is TRUE.}
}
\value{
  Invisibly returns the path to the created project directory. Creates a
  .locusPackRat subdirectory with the following structure:
  \itemize{
    \item input/ - Contains genes.csv or regions.csv and orthology.csv
    \item supplementary/ - Directory for additional data tables
    \item output/ - Directory for generated outputs
    \item config.json - Project configuration file
  }
}
\examples{
# Initialize with gene list
genes <- data.frame(
  gene_symbol = c("BRCA1", "TP53", "EGFR"),
  expression = c(100, 250, 50)
)
initPackRat(genes, mode = "gene", species = "human",
            genome = "hg38", project_dir = "my_project")

# Initialize with QTL regions
regions <- data.frame(
  chr = c(5, 10, 12),
  start = c(10000000, 50000000, 75000000),
  end = c(15000000, 55000000, 80000000),
  trait = c("height", "weight", "BMI")
)
initPackRat(regions, mode = "region", species = "mouse",
            genome = "mm39", project_dir = "qtl_project")
}
\seealso{
  \code{\link{addRatTable}}, \code{\link{generateGeneSheet}},
  \code{\link{listPackRatTables}}
}
\author{
  Brian Gural, Todd Kimball, Anh Luu, Christoph D. Rau
}