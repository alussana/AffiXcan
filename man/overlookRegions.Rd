% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AffiXcan.R
\name{overlookRegions}
\alias{overlookRegions}
\title{Count the number of genomic regions on which the TBA was computed}
\usage{
overlookRegions(tbaPaths)
}
\arguments{
\item{tbaPaths, }{A vector of strings, which are the paths to
MultiAssayExperiment RDS files containing the tba values}
}
\value{
An integer, that is the summation of length(assays()) of every
MultiAssayExperiment RDS object indicated in the param tbaPaths
}
\description{
Count the number of genomic regions on which the TBA was computed
}
\examples{
if (interactive()) {
testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
package="AffiXcan")

regionsCount <- overlookRegions(tbaPaths=testingTbaPaths)
}
}
