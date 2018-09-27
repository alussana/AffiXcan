#' Expression data of two genes for 229 individuals
#'
#' Toy dataset used in examples to describe affiXcanTrain() function.
#'
#' The data consist in a SummarizedExperiment object that contains a matrix of
#' expression values (RPKM) of two randomly chosen genes ("ENSG00000139269.2"
#' and "ENSG00000256377.1") for 229 individuals of european descent. The data
#' was retrieved subsetting the RNA-sequencing data of EBV-transformed
#' lymphocytes from the GEUVADIS public dataset (see 'source')
#'
#' @docType data
#'
#' @usage data(exprMatrix)
#'
#' @format An object of class \code{SummarizedExperiment}.
#'
#' @keywords datasets
#'
#' @source \url{https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-3/}
#'
#' @examples
#' library(SummarizedExperiment)
#' data(exprMatrix)
#' toyExpressionMatrix <- assays(exprMatrix)$values
"exprMatrix"

#' Covariates of the population structure for 229 individuals
#'
#' Toy data used in examples to describe affiXcanTrain() function.
#'
#' This object consists in a data.frame where columns are the first three
#' principal components of the population genetic structure and rows are
#' individuals' IDs. These individuals are the same whom expression values are
#' stored in the expression matrix (see help(exprMatrix) )
#'
#' Genotypes of the individuals were downloaded from the GEUVADIS public dataset
#' (\url{https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/}) in vcf format.
#' Following L. Price et al.
#' (\url{https://www.sciencedirect.com/science/article/pii/S0002929708003534}),
#' long range linkage disequilibrium (LRLD) regions were first filtered out
#' with vcf-tools. Then, following J. Novembre et al.
#' (\url{www.nature.com/articles/nature07331}), non-common alleles (MAF < 0.05)
#' were filtered out with vcftools and LD pruning was performed with plink.
#' Finally, principal components were computed with eigenstrat. 
#'
#' @docType data
#'
#' @usage data(trainingCovariates)
#'
#' @format An object of class \code{data.frame}
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
#' functions.
#'
#' The object consists in a data.frame with two columns: for every
#' "EXPRESSED_REGION" are listed the associated "REGULATORY_REGION"(s).
#' For the illustrative purpose there are only two expressed genes:
#' "ENSG00000139269.2" and "ENSG00000256377.1" (the same names are used in the
#' SummarizedExperiment object containing the expression matrix, see
#' help(exprMatrix) ). 
#'
#' For every gene, the associated regulatory regions were retrieved among the
#' enhancers predicted by preSTIGE, a tool developed by O. Corradin et al.
#' (\url{https://genome.cshlp.org/content/24/1/1.full}), in EBV-transformed
#' lymphocytes cell lines. The name of the regulatory region being identical to
#' the name of the expressed gene means that the regulatory region refers to
#' a genomic window that spans at least for 2 Kbp and includes the most upstream
#' TSS of the gene.
#'
#' @docType data
#'
#' @usage data(regionAssoc)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
#' @examples
#' data(regionAssoc)
#' head(regionAssoc)
"regionAssoc"
