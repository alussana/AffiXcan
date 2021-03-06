% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasetsDoc.R
\docType{data}
\name{trainingCovariates}
\alias{trainingCovariates}
\title{Covariates of the population structure for 229 individuals}
\format{
An object of class \code{data.frame}
}
\usage{
data(trainingCovariates)
}
\description{
Toy data used in examples to describe affiXcanTrain() function.
}
\details{
This object consists in a data.frame where columns are the first three
principal components of the population genetic structure and rows are
individuals' IDs. These individuals are the same whom expression values are
stored in the expression matrix (see help(exprMatrix) )

Genotypes of the individuals were downloaded from the GEUVADIS public dataset
(\url{https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/}) in vcf format.
Following L. Price et al.
(\url{https://www.sciencedirect.com/science/article/pii/S0002929708003534}),
long range linkage disequilibrium (LRLD) regions were first filtered out
with vcf-tools. Then, following J. Novembre et al.
(\url{www.nature.com/articles/nature07331}), non-common alleles (MAF < 0.05)
were filtered out with vcftools and LD pruning was performed with plink.
Finally, principal components were computed with eigenstrat.
}
\examples{
data(trainingCovariates)
head(trainingCovariates)
}
\keyword{datasets}
