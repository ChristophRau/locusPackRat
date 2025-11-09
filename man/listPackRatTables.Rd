\name{listPackRatTables}
\alias{listPackRatTables}
\title{List Supplementary Tables in locusPackRat Project}
\description{
  Lists all supplementary data tables that have been added to a locusPackRat
  project using addRatTable().
}
\usage{
listPackRatTables(project_dir = ".")
}
\arguments{
  \item{project_dir}{Path to the project directory. Defaults to current directory.}
}
\value{
  A character vector of supplementary table names, or NULL if no tables exist.
  Returns an empty vector with a message if the project doesn't exist.
}
\examples{
# List all supplementary tables
tables <- listPackRatTables("my_project")
print(tables)

# Use to check what's available before generating output
available <- listPackRatTables("my_project")
if ("expression" \%in\% available) {
  generateGeneSheet(include_supplementary = "expression",
                    project_dir = "my_project")
}
}
\seealso{
  \code{\link{addRatTable}}, \code{\link{generateGeneSheet}}
}
\author{
  Brian Gural, Todd Kimball, Anh Luu, Christoph D. Rau
}