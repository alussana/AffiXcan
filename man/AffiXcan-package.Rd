\name{AffiXcan-package}
\alias{AffiXcan-package}
\alias{AffiXcan}
\docType{package}
\title{
A functional approach to impute GReX
}
\description{
Impute a GReX (Genetically Regulated Expression) for a set of genes in a sample 
of individuals, using a method based on the Total Binding Affinity (TBA) score.
Statistical models to impute GReX can be trained on a training dataset where
the real total expression values are known.
}
\details{
For every gene a linear regression model can be fitted on a training set of
individuals where the real total expression values, the Total Binding Affinity
(TBA) values for a set of genomic regions, and the covariates of the population
genetic structure, are known. AffiXcan performs a principal component analysis
(PCA) on the TBA values to fit a linear model using the formula: GReX ~ PC1 +
PC2 + PC3 + ... + COV1 + COV2 + COV3. Associations between the expressed genes
and the regulatory regions, on which the TBA values have to be computed, are
needed. TBA can be computed using "vcf_rider" software (see references)
}
\author{
Alessandro Lussana <alessandro.lussana@protonmail.com>
Maintainer: Alessandro Lussana <alessandro.lussana@protonmail.com>
}
\references{
\url{https://github.com/vodkatad/vcf_rider}
}
\keyword{ package }
\examples{
trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
package="AffiXcan")

data(exprMatrix)
data(regionAssoc)
data(trainingCovariates)

assay <- "values"

training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
varExplained=80, scale=TRUE, cores=1)

testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
package="AffiXcan")

exprmatrix <- affiXcanImpute(tbaPaths=testingTbaPaths, pca=training$pca,
bs=training$bs, scale=TRUE, cores=1)
}
