% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AffiXcan.R
\name{affiXcanTrain}
\alias{affiXcanTrain}
\title{Train the model needed to impute a GReX for each gene}
\usage{
affiXcanTrain(
  exprMatrix,
  assay,
  tbaPaths,
  regionAssoc,
  cov = NULL,
  varExplained = 80,
  scale = TRUE,
  BPPARAM = bpparam(),
  kfold = 1
)
}
\arguments{
\item{exprMatrix}{A SummarizedExperiment object containing expression data}

\item{assay}{A string with the name of the object in
SummarizedExperiment::assays(exprMatrix) that contains expression values}

\item{tbaPaths}{A vector of strings, which are the paths to
MultiAssayExperiment RDS files containing the tba values}

\item{regionAssoc}{A data.frame with the association between regulatory
regions and expressed genes and with colnames = c("REGULATORY_REGION",
"EXPRESSED_REGION")}

\item{cov}{Optional. A data.frame with covariates values for the population 
structure where the columns are the PCs and the rows are the individual IIDs
Default is NULL}

\item{varExplained}{An integer between 0 and 100; varExplained=80 means that
the principal components selected to fit the models must explain at least
80 percent of variation of TBA values; default is 80}

\item{scale}{A logical; if scale=FALSE the TBA values will be only centered,
not scaled before performing PCA; default is TRUE}

\item{BPPARAM}{A BiocParallelParam object. Default is bpparam(). For
details on BiocParallelParam virtual base class see 
browseVignettes("BiocParallel")}

\item{kfold}{An integer. The k definition of k-fold cross-validation on the
training dataset. Default is 1. This argument controls the behavior of the
function in the following way:
     
\itemize{
 \item kfold<2: train the models on the whole dataset, not performing the
  cross-validation
 \item kfold>=2: perform the k-fold cross-validation
}}
}
\value{
The output depends on the parameter kfold

If kfold<2: a list containing three objects: pca, bs, regionsCount
\itemize{
 \item pca: A list containing lists named as the
 MultiAssayExperiment::experiments() found in the MultiAssayExperiment
 objects listed in the param tbaPaths. Each of these lists contain two
 objects:
 \itemize{
     \item eigenvectors: A matrix containing eigenvectors for those principal
     components of the TBA selected according to the param varExplained
     \item pcs: A matrix containing the principal components values of the 
     TBA selected according to the param varExplained
     \item eigenvalues: A vector containing eigenvalues for those principal
     components of the TBA selected according to the param varExplained
 }
 \item bs: A list containing lists named as the EXPRESSED_REGIONS found in
 the param regionAssoc that have a correspondent rowname in the expression
 values stored SummarizedExperiment::assays(exprMatrix)$assay.
 Each of the lists in bs contains three objects:
 \itemize{
     \item coefficients: The coefficients of the principal components used in the
     model, completely similar to the "coefficients" from the results of lm()
     \item p.val: The uncorrected anova pvalue of the model
     \item r.sq: The coefficient of determination between the real total expression
     values and the imputed GReX, retrived from summary(model)$r.squared
     \item corrected.p.val: The p-value of the model, corrected for multiple
     testing with benjamini-hochberg procedure
 }
 \item regionsCount: An integer, that is the number of genomic regions taken into
 account during the training phase
}

If kfold>=2: a list containing k-fold objects, named from 1 to kfold and
corresponding to the different cross-validations [i]; each one of these
objects is a list containing lists named as the expressed gene IDs [y] (i.e.
the rownames() of the object in SummarizedExperiment::assays(exprMatrix)
containing the expression values), for which a GReX could be imputed. Each of
these inner lists contain two objects:
\itemize{     
 \item rho: the pearson's correlation coefficient (R) between the real
 expression values and the imputed GReX for the cross-validation i on
 the expressed gene y, computed with cor()
 \item rho.sq: the coefficient of determination (R^2) between the real
 expression values and the imputed GReX for the cross-validation i on
 the expressed gene y, computed as pearson^2
 \item cor.test.p.val: the p-value of the cor.test() between the real expression values
 and the imputed GReX for the cross-validation i on the expressed gene y
 \item model.p.val: The uncorrected anova pvalue of the model
 \item model.corrected.p.val: The p-value of the model, corrected for
 multiple testing with benjamini-hochberg procedure
 \item model.r.sq: the model's coefficient of determination (R^2) on the
 training data
}
}
\description{
Train the model needed to impute a GReX for each gene
}
\examples{
if(interactive()) {
trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
package="AffiXcan")

data(exprMatrix)
data(regionAssoc)
data(trainingCovariates)

assay <- "values"

training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
varExplained=80, scale=TRUE, kfold=3)
}
}
