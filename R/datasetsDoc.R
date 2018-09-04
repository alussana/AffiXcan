#' Expression data of four genes for 230 individuals
#'
#' Toy dataset used in examples to describe affiXcanTrain() function
#'
#' @docType data
#'
#' @usage data(exprMatrix)
#'
#' @format An object of class \code{"SummarizedExperiment"}.
#'
#' @keywords datasets
#'
#' @examples
#' library(SummarizedExperiment)
#' data(exprMatrix)
#' toyExpressionMatrix <- assays(exprMatrix)$values
"exprMatrix"

#' Covariates of the population structure for 230 individuals
#'
#' Toy data used in examples to describe affiXcanTrain() function;
#' Columns are the first three principal components of the population structure
#' Rows are individials
#'
#' @docType data
#'
#' @usage data(trainingCovariates)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @keywords datasets
#'
#' @examples
#' data(trainingCovariates)
#' head(trainingCovariates)
"trainingCovariates"

#' Associations between regulatory regions and expressed genes
#'
#' Toy data used in examples to describe affiXcanTrain() and affiXcanImpute() 
#' functions; it has two columns: for every EXPRESSED_REGION are listed the
#' associated REGULATORY_REGION(s)
#'
#' @docType data
#'
#' @usage data(regionAssoc)
#'
#' @format An object of class \code{"data.frame"}
#'
#' @keywords datasets
#'
#' @examples
#' data(regionAssoc)
#' head(regionAssoc)
"regionAssoc"
