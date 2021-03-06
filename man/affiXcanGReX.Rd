% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AffiXcan.R
\name{affiXcanGReX}
\alias{affiXcanGReX}
\title{Compute a GReX from variables and their coefficients for a set of genes}
\usage{
affiXcanGReX(affiXcanTraining, pcs, BPPARAM = bpparam())
}
\arguments{
\item{affiXcanTraining}{The returning object from affiXcanTrain()}

\item{pcs}{A list, which is the returning object from affiXcanPcs()}

\item{BPPARAM}{A BiocParallelParam object. Default is bpparam(). For
details on BiocParallelParam virtual base class see 
browseVignettes("BiocParallel")}
}
\value{
A SummarizedExperiment object containing the imputed GReX values
}
\description{
Compute a GReX from variables and their coefficients for a set of genes
}
\examples{
if (interactive()) {
trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
package="AffiXcan")

data(exprMatrix)
data(regionAssoc)
data(trainingCovariates)

assay <- "values"

training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
varExplained=80, scale=TRUE)

testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
package="AffiXcan")

pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, affiXcanTraining=training,
scale=TRUE)

exprmatrix <- affiXcanGReX(affiXcanTraining=training, pcs=pcs)
}
}
