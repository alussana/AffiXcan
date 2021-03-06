% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AffiXcan.R
\name{computePcs}
\alias{computePcs}
\title{Compute a matrix product between variables and eigenvectors}
\usage{
computePcs(region, tbaMatrix, scale, pca)
}
\arguments{
\item{region}{A string, which is the name of the object in the list
MultiAssayExperiment::experiments(tbaMatrix) that contains the TBA values of
the genomic region of interest (see the param tbaMatrix)}

\item{tbaMatrix}{A MultiAssayExperiment object containing the TBA values}

\item{scale}{A logical; if scale=FALSE the TBA values will be only centered,
not scaled before performing PCA}

\item{pca}{The returningObject$pca from affiXcanTrain()}
}
\value{
A data.frame containing the principal components values of the TBA
in a certain genomic region, as the result of the matrix product between the
TBA values and the eigenvectors
}
\description{
Compute a matrix product between variables and eigenvectors
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

region <- "ENSG00000256377.1"

tbaMatrixMAE <- readRDS(system.file("extdata","training.tba.toydata.rds",
package="AffiXcan"))

pca <- training$pca
pcs <- computePcs(region=region, tbaMatrix=tbaMatrixMAE, scale=TRUE, pca=pca)
}
}
