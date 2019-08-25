#' Train the model needed to impute a GReX for each gene 
#'
#' @param exprMatrix A SummarizedExperiment object containing expression data
#' @param assay A string with the name of the object in
#' SummarizedExperiment::assays(exprMatrix) that contains expression values 
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param regionAssoc A data.frame with the association between regulatory
#' regions and expressed genes and with colnames = c("REGULATORY_REGION",
#' "EXPRESSED_REGION")
#' @param cov Optional. A data.frame with covariates values for the population 
#' structure where the columns are the PCs and the rows are the individual IIDs
#' Default is NULL
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values; default is 80
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA; default is TRUE
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#' @param kfold An integer. The k definition of k-fold cross-validation on the
#' training dataset. Default is 1. This argument controls the behavior of the
#' function in the following way:
#'      
#' \itemize{
#'  \item kfold<2: train the models on the whole dataset, not performing the
#'   cross-validation
#'  \item kfold>=2: perform the k-fold cross-validation
#' }
#'
#' @return The output depends on the parameter kfold
#' 
#' If kfold<2: a list containing three objects: pca, bs, regionsCount
#' \itemize{
#'  \item pca: A list containing lists named as the
#'  MultiAssayExperiment::experiments() found in the MultiAssayExperiment
#'  objects listed in the param tbaPaths. Each of these lists contain two
#'  objects:
#'  \itemize{
#'      \item eigenvectors: A matrix containing eigenvectors for those principal
#'      components selected according to the param varExplained
#'      \item pcs: A matrix containing the principal components values selected
#'      according to the param varExplained
#'  }
#'  \item bs: A list containing lists named as the REGULATORY_REGIONS found in
#'  the param regionAssoc that have a correspondent colname in the experiments
#'  stored in MultiAssayExperiment objects listed in the param tbaPaths.
#'  Each of the lists in bs contains three objects:
#'  \itemize{
#'      \item coefficients: The coefficients of the principal components used in the
#'      model, completely similar to the "coefficients" from the results of lm()
#'      \item pval: The uncorrected anova pvalue of the model, retrieved from
#'      anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'      If cov==NULL, i.e. no covariates for the population structure have
#'      been provided, pval will be NA
#'      \item r.sq: The coefficient of determination between the real total expression
#'      values and the imputed GReX, retrived from summary(model)$r.squared
#'  }
#'  \item regionsCount: An integer, that is the number of genomic regions taken into
#'  account during the training phase
#' }
#'
#' If kfold>=2: a list containing k-fold objects, named from 1 to kfold and
#' corresponding to the different cross-validations [i]; each one of these
#' objects is a list containing lists named as the expressed gene IDs [y] (i.e.
#' the rownames() of the object in SummarizedExperiment::assays(exprMatrix)
#' containing the expression values), for which a GReX could be imputed. Each of
#' these inner lists contain two objects:
#' \itemize{     
#'  \item pearson: the pearson's correlation coefficient (R) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed with cor()
#'  \item rSquared: the coefficient of determination (R^2) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed as pearson^2
#' }
#' @import MultiAssayExperiment SummarizedExperiment BiocParallel crayon
#' @export
#'
#' @examples
#' if(interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, kfold=3)
#' }
affiXcanTrain <- function(exprMatrix, assay, tbaPaths, regionAssoc, cov=NULL,
    varExplained=80, scale=TRUE, BPPARAM=bpparam(), kfold=1) {
    regionsCount <- overlookRegions(tbaPaths)
    sampleNames <- colnames(exprMatrix)
    nSamples <- length(sampleNames)
    sampGroups <- subsetKFold(kfold, nSamples)
    for (i in seq(length(sampGroups))) {
        sampGroups[[i]] <- colnames(exprMatrix)[sampGroups[[i]]]
    }
    trainingOutput <- vector('list', length=kfold)
    
    for (i in seq(kfold)) {
        if (kfold > 1) {
            cat(green(bold("\nAffiXcan")),": ",cyan("Performing Training With",
                kfold, "Fold Cross-Validation (",i,"/",kfold,")\n"), sep="")
            testingSamples <- sampGroups[[i]]
            trainingSamples <- sampleNames[!sampleNames %in% testingSamples]
        }
        else {
            cat(green(bold("\nAffiXcan")),": ", cyan("Training The Models\n"),
                sep="")
            trainingSamples <- sampleNames
        }

        cat(bold(cyan("\t-->")),
            green("Performing Principal Components Analysis\n"))
        pca <- affiXcanPca(tbaPaths, varExplained, scale, regionsCount,
                           BPPARAM, trainingSamples)
        cat(bold(cyan("\t-->")), green("Training Coefficients\n"))
        bs <- affiXcanBs(exprMatrix, assay, regionAssoc, pca, cov, BPPARAM,
                         trainingSamples)
        affiXcanTraining <- list(pca=pca, bs=bs, regionsCount=regionsCount)

        if (kfold > 1) {
            cat(bold(cyan("\t-->")),
                green("Computing Principal Components On Testing Set\n"))
            pcs <- affiXcanPcs(tbaPaths, affiXcanTraining, scale, BPPARAM,
                               testingSamples)
            cat(bold(cyan("\t-->")), green("Imputing GReX Values\n"))
            imputedExpr <- affiXcanGReX(affiXcanTraining, pcs, BPPARAM)
            cat(bold(cyan("\t-->")), green("Computing Cross-Validated R^2\n"))
            correlation <- computeRSquared(exprMatrix, imputedExpr, assay,
                                         testingSamples, BPPARAM)
            trainingOutput[[i]] <- correlation
        }
        else {
            trainingOutput <- affiXcanTraining
        }
        cat(bold(cyan("\tDone!\n")))
    }
    
    return(trainingOutput)
}

#' Split the numbers from 1 to n in k equally-sized lists
#'
#' @param k An integer. The number of groups in which the first n natural
#' numbers are to be splitted
#' @param n An integer. Defines the interval 1..n 
#'
#' @return A list of lists; the first natural n numbers equally distributed 
#' across k lists
#'
#' @export
#'
#' @examples
#' sampGroups <- subsetKFold(k=5, n=23)
subsetKFold <- function(k, n) {
    sampGroups <- as.vector(seq(n))
    sampGroups <- split(sampGroups, seq_along(sampGroups)%%k + 1)
    return(sampGroups)
}

#' Perform a PCA on each experiment found in MultiAssayExperiment objects
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values; default is 80
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA; default is TRUE
#' @param regionsCount An integer, that is the summation of length(assays()) of
#' every MultiAssayExperiment RDS object indicated in the param tbaPaths; it is
#' the returning value from overlookRegions()
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#' @param trainingSamples A vector of strings. The identifiers (e.g. row names
#' of MultiAssayExperiment objects from tbaPaths) of the samples that have to
#' be considered in the training phase, and not used for the cross-validation 
#'
#' @return pca: A list containing lists named as the 
#' MultiAssayExperiment::experiments() found in the MultiAssayExperiment objects
#' listed in the param tbaPaths. Each of these lists contain two objects:
#' \itemize{
#'  \item eigenvectors: A matrix containing eigenvectors for those principal
#'  components selected according to the param varExplained
#'  \item pcs: A matrix containing the principal components values selected
#'  according to the param varExplained
#' }
#' @import MultiAssayExperiment BiocParallel 
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(exprMatrix)
#'
#' tbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' regionsCount <- overlookRegions(tbaPaths)
#'
#' sampleNames <- colnames(exprMatrix)
#' nSamples <- length(sampleNames)
#' sampGroups <- subsetKFold(k=3, n=nSamples)
#' for (i in seq(length(sampGroups))) {
#'      sampGroups[[i]] <- colnames(exprMatrix)[sampGroups[[i]]]
#' }
#' 
#' testingSamples <- sampGroups[[1]]
#' trainingSamples <- sampleNames[!sampleNames %in% testingSamples]
#'
#' pca <- affiXcanPca(tbaPaths=tbaPaths, varExplained=80, scale=TRUE,
#' regionsCount=regionsCount, trainingSamples=trainingSamples)
#' }
affiXcanPca <- function(tbaPaths, varExplained=80, scale=TRUE, regionsCount,
                        BPPARAM=bpparam(), trainingSamples) {

    pca <- vector("list", regionsCount)
    index <- 0

    for(i in seq(1,length(tbaPaths))) {

        tbaMatrixMAE <- readRDS(tbaPaths[i])
        tbaMatrixMAE <- MultiAssayExperiment::subsetByRow(tbaMatrixMAE,
                                                          trainingSamples)
        tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)

        newPca <- BiocParallel::bplapply(X=tbaMatrix, FUN=computePca,
        varExplained, scale, BPPARAM=BPPARAM)

        for(l in seq(1,length(newPca))) {
            pca[index + l] <- newPca[l]
            names(pca)[index + l] <- names(newPca)[l]
        }

        index <- index + length(newPca)
    
    }

    return(pca)

}

#' Perform a PCA on a matrix where columns are variables
#'
#' @param data A matrix containing the TBA values for a certain genomic region;
#' columns are PWMs, rows are individuals IIDs
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values; default is 80
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA; default is TRUE
#'
#' @return A list containing two objects:
#' \itemize{ 
#'  \item eigenvectors: a matrix containing eigenvectors for those principal
#'  components selected according to the param varExplained
#'  \item pcs: a matrix containing the principal components values selected
#'  according to the param varExplained
#' }
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(exprMatrix)
#' sampleNames <- colnames(exprMatrix)
#' nSamples <- length(sampleNames)
#' sampGroups <- subsetKFold(k=3, n=nSamples)
#' for (i in seq(length(sampGroups))) {
#'      sampGroups[[i]] <- colnames(exprMatrix)[sampGroups[[i]]]
#' }
#' testingSamples <- sampGroups[[i]]
#' trainingSamples <- sampleNames[!sampleNames %in% testingSamples]
#'
#' tbaMatrixMAE <- readRDS(system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan"))
#' tbaMatrixMAE <- MultiAssayExperiment::subsetByRow(tbaMatrixMAE,
#'                                                   trainingSamples)
#' tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)
#' tba <- tbaMatrix$ENSG00000256377.1
#'
#' pca <- computePca(data=tba, varExplained=80, scale=TRUE)
#' 
#' }
computePca <- function(data, varExplained=80, scale=TRUE) {

    tryCatch (
        {
            notScalablePwm <- vector('character')
            for(pwm in colnames(data)) {
                if(TRUE %in% is.na(data[[pwm]]) || var(data[[pwm]])==0) {        
                    notScalablePwm <- c(notScalablePwm, pwm)
                    data[[pwm]] <- NULL
                }
            }

            dataPrcomp <- prcomp(data, center = TRUE, scale. = scale)
            dataPcs <- dataPrcomp$x
            pcVar <- dataPrcomp$sdev^2
            totVar <- sum(pcVar)
            eigenvalues <- pcVar * 100 / totVar 
            eigenvectors <- dataPrcomp$rotation
            nWantedPcs <- which.max(cumsum(eigenvalues) >= varExplained)
            wantedPcs <- subset(dataPcs, select=colnames(dataPcs)[1:nWantedPcs])
            wantedEigenvectors <- subset(eigenvectors,
                select=colnames(eigenvectors)[1:nWantedPcs])

            for(pwm in notScalablePwm) {
                newLine <- matrix(rep(0, nWantedPcs), ncol=1)
                wantedEigenvectors <- rbind(wantedEigenvectors, newLine)
                rownames(wantedEigenvectors)[[nrow(wantedEigenvectors)]] <- pwm
            }

            wantedEigenvectors <- as.data.frame(wantedEigenvectors)
            wantedEigenvectors <- 
                wantedEigenvectors[order(rownames(wantedEigenvectors)),]
            wantedEigenvectors <- as.matrix(wantedEigenvectors)

            return(list(eigenvectors=wantedEigenvectors, pcs=wantedPcs))

        },
        error=function(e) {
            wantedEigenvectors <- NA
            wantedPcs <- NA
            return(list(eigenvectors=wantedEigenvectors, pcs=wantedPcs))
        }
    )
}

#' Fit a linear model and compute ANOVA p value
#'
#' @param exprMatrix A SummarizedExperiment object containing expression data
#' @param assay A string with the name of the object in
#' SummarizedExperiment::assays(exprMatrix) that contains expression values
#' @param regionAssoc A data.frame with the associations between regulatory
#' regions and expressed genes, and with colnames = c("REGULATORY_REGION",
#' "EXPRESSED_REGION")
#' @param pca A list, which is the returningObject$pca from affiXcanPca()
#' @param cov Optional; a data.frame with covariates values for the population
#' structure where the columns are the PCs and the rows are the individual IDs;
#' default is NULL
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#' @param trainingSamples A vector of strings. The identifiers (e.g. row names
#' of MultiAssayExperiment objects from tbaPaths) of the samples that have to
#' be considered in the training phase, and not used for the cross-validation
#'
#' @return A list containing lists named as the REGULATORY_REGIONS found in the
#' param regionAssoc that have a correspondent name in the param pca.
#' Each of these lists contain three objects:
#' \itemize{
#'  \item coefficients: An object containing the coefficients of the principal
#'  components used in the model, completely similar to the "coefficients"
#'  from the results of lm()
#'  \item pval: The uncorrected anova pvalue of the model, retrieved from
#'  anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'  If cov==NULL, i.e. no covariates for the population structure have
#'  been provided, pval will be NA
#'  \itwm r.sq: The coefficient of determination between the real total expression
#'  values and the imputed GReX, retrived from summary(model)$r.squared
#' }
#' @import SummarizedExperiment BiocParallel
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(exprMatrix)
#' data(trainingCovariates)
#' data(regionAssoc)
#'
#' tbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' regionsCount <- overlookRegions(tbaPaths)
#' assay <- "values"
#'
#' sampleNames <- colnames(exprMatrix)
#' nSamples <- length(sampleNames)
#' sampGroups <- subsetKFold(k=5, n=nSamples)
#' for (i in seq(length(sampGroups))) {
#'      sampGroups[[i]] <- colnames(exprMatrix)[sampGroups[[i]]]
#' }
#' 
#' testingSamples <- sampGroups[[1]]
#' trainingSamples <- sampleNames[!sampleNames %in% testingSamples]
#'
#' pca <- affiXcanPca(tbaPaths=tbaPaths, varExplained=80, scale=TRUE,
#' regionsCount=regionsCount, trainingSamples=trainingSamples)
#'
#' bs <- affiXcanBs(exprMatrix=exprMatrix, assay=assay, regionAssoc=regionAssoc,
#' pca=pca, cov=trainingCovariates, trainingSamples=trainingSamples)
#' }
affiXcanBs <- function(exprMatrix, assay, regionAssoc, pca, cov=NULL,
                       BPPARAM=bpparam(), trainingSamples) {

    if (is.null(cov)==FALSE) {
        cov <- cov[rownames(cov) %in% trainingSamples,]
    }
    
    expr <- SummarizedExperiment::assays(exprMatrix)[[assay]]
    expr <- expr[,colnames(expr) %in% trainingSamples]
    tDfExprMatrix<-t(as.data.frame(expr))

    expressedRegions <- as.vector(unique(regionAssoc$EXPRESSED_REGION))
    assocList <- BiocParallel::bplapply(X=expressedRegions, FUN=assoc2list,
        regionAssoc, BPPARAM=BPPARAM)
    names(assocList) <- expressedRegions

    bs <- BiocParallel::bplapply(X=assocList, FUN=computeBs, pca,
        tDfExprMatrix, cov, BPPARAM=BPPARAM)

    ## Remove ensg for which the model is NA
    for(i in names(bs)) {
        if(is.na(bs[[i]]$r.sq)) {
            bs[[i]] <- NULL
        }
    }

    ## ANOVA is not performed anymore (changes in version 1.3.1)

    ## Compute multitest pvalue correction (benjamini-hochberg)
    #pvalues <- vector("numeric", length=length(bs))
    #
    #for(i in seq(1:length(pvalues))) {
    #    pvalues[[i]] <- bs[[i]]$pval
    #}
    #
    #multi.test.pvals <- p.adjust(pvalues, method="BH")
    #
    #for(i in seq(1:length(bs))) {
    #    bs[[i]]$correctedP <- multi.test.pvals[[i]]
    #}
    
    return(bs)
}

#' Reorganize associations table in a list
#'
#' @param gene A string; the name of an expressed gene
#' @param regionAssoc A data.frame with the associations between regulatory
#' regions and expressed genes, and with colnames = c("REGULATORY_REGION",
#' "EXPRESSED_REGION")
#'
#' @return A list of data frames, each of which has the same structure of the
#' param regionAssoc, except that contains the information relative to one
#' expressed gene
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(regionAssoc)
#' expressedRegions <- as.list(as.vector(unique(regionAssoc$EXPRESSED_REGION)))
#' gene <- expressedRegions[[1]]
#' assocList <- assoc2list(gene, regionAssoc)
#' }
assoc2list <- function(gene, regionAssoc) {
    return(regionAssoc[regionAssoc$EXPRESSED_REGION==gene,])
}

#' Fit a linear model to impute a GReX for a certain gene
#'
#' @param assocRegions A data.frame with the associations between regulatory
#' regions and one expressed gene, and with colnames = c("REGULATORY_REGION",
#' "EXPRESSED_REGION")
#' @param pca The returningObject$pca from affiXcanTrain()
#' @param expr A matrix containing the real total expression values, where the
#' columns are genes and the rows are individual IIDs
#' @param covariates Optrional; a data.frame with covariates values for the
#' population structure where the columns are the PCs and the rows are the
#' individual IIDs; default is NULL
#'
#' @return A list containing three objects:
#' \itemize{
#'  \item coefficients: An object containing the coefficients of the principal
#'  components used in the model, completely similar to the "coefficients"
#'  from the results of lm()
#'  \item pval: The uncorrected anova pvalue of the model, retrieved from
#'  anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'  If cov==NULL, i.e. no covariates for the population structure have
#'  been provided, pval will be NA
#'  \item r.sq: The coefficient of determination between the real total expression
#'  values and the imputed GReX, retrived from summary(model)$r.squared
#' }
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(exprMatrix)
#' data(trainingCovariates)
#' data(regionAssoc)
#'
#' tbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' regionsCount <- overlookRegions(tbaPaths)
#' assay <- "values"
#'
#' sampleNames <- colnames(exprMatrix)
#' nSamples <- length(sampleNames)
#' sampGroups <- subsetKFold(k=5, n=nSamples)
#' for (i in seq(length(sampGroups))) {
#'      sampGroups[[i]] <- colnames(exprMatrix)[sampGroups[[i]]]
#' }
#' 
#' testingSamples <- sampGroups[[1]]
#' trainingSamples <- sampleNames[!sampleNames %in% testingSamples]
#'
#' pca <- affiXcanPca(tbaPaths=tbaPaths, varExplained=80, scale=TRUE,
#' regionsCount=regionsCount, trainingSamples=trainingSamples)
#'
#' cov <- trainingCovariates
#' cov <- cov[rownames(cov) %in% trainingSamples,]
#'
#' expr <- SummarizedExperiment::assays(exprMatrix)[[assay]]
#' expr <- expr[,colnames(expr) %in% trainingSamples]
#' expr <- t(as.data.frame(expr))
#'
#' expressedRegions <- as.vector(unique(regionAssoc$EXPRESSED_REGION))
#' assocList <- BiocParallel::bplapply(X=expressedRegions, FUN=assoc2list,
#' regionAssoc)
#' names(assocList) <- expressedRegions
#' assocRegions <- assocList[[1]]
#'
#' bs <- computeBs(assocRegions=assocRegions, pca=pca, expr=expr,
#' covariates=cov)
#' }
computeBs <- function(assocRegions, pca, expr, covariates=NULL) {
    expressedRegion <- as.vector(unique(assocRegions$EXPRESSED_REGION))
    expression <- as.vector(as.data.frame(expr)[[expressedRegion]])
    tryCatch (
        {
            myPcs <- as.data.frame(matrix(nrow = nrow(expr)))
            myPcsColnames <- vector("character")

            for(regulRegion in assocRegions$REGULATORY_REGION) {

                ## Discards the pcs for those regions that are NA in the model
                if(is.null(pca[[regulRegion]]$pcs[[1]]) ||
                    is.na(pca[[regulRegion]]$pcs[[1]])) {

                    next

                } else {

                    nNames <- length(colnames(pca[[regulRegion]]$pcs))
                    regionColnames <- colnames(pca[[regulRegion]]$pcs)
                    regionColnames <- paste(regionColnames, regulRegion,
                                                                    sep="@")
                    myPcsColnames <- c(myPcsColnames, regionColnames)
                    regulRegionPcs <- pca[[regulRegion]]$pcs
                    myPcs <- cbind(myPcs, regulRegionPcs)

                }
            }

            if (is.null(covariates)==FALSE) {
                myPcsColnames <- c(myPcsColnames, names(covariates))
                myPcs <- cbind(myPcs, covariates)
            }

            myPcs <- as.matrix(subset(myPcs, select=-V1))
            colnames(myPcs) = myPcsColnames

            model <- lm(expression~myPcs)
            
            if (is.null(covariates)==FALSE) {
                modelReduced <- lm(expression~as.matrix(covariates))
                a <- anova(model, modelReduced, test="F")
                pval <- a$'Pr(>F)'[2]
            } else {
                pval <- NA
            }

            if (is.null(covariates)==FALSE) {
                betas <- model$coefficients[seq(1,
                    (length(model$coefficients)-length(names(covariates))))]
            } else {
                betas <- model$coefficients
            }

            r.sq <- summary(model)$r.squared
            results <- list(coefficients=betas, pval=pval, r.sq=r.sq)
            return(results)
        },
        error = function(e) {
            ## TODO: throw warning
            betas <- NA
            pval <- NA
            r.sq <- NA
            results <- list(coefficients=betas, pval=pval, r.sq=r.sq)
            return(results)
        }
    )
}

#' Count the number of genomic regions on which the TBA was computed
#'
#' @param tbaPaths, A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#'
#' @return An integer, that is the summation of length(assays()) of every
#' MultiAssayExperiment RDS object indicated in the param tbaPaths
#' @import MultiAssayExperiment
#' @export
#'
#' @examples
#' if (interactive()) {
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#'
#' regionsCount <- overlookRegions(tbaPaths=testingTbaPaths)
#' }
overlookRegions <- function(tbaPaths) {

    regionsCount <- 0 
    
    for(path in tbaPaths) {
        tba <- readRDS(path)
        assaysNum <- length(assays(tba))
        regionsCount <- regionsCount + assaysNum
    }

    return(regionsCount)
    
}

#' Compute PCs in MultiAssayExperiment objects using eigenvectors given by user
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param affiXcanTraining The returning object from affiXcanTrain()
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#' @param testingSamples A vector of strings. The identifiers (e.g. row names
#' of MultiAssayExperiment objects from tbaPaths) of the samples that have not
#' been considered in the training phase, to be used in the cross-validation;
#' default is NULL; if is.null(testingSamples)==TRUE then no filtering is
#' performed
#'
#' @return A list of matrices containing the principal components values of TBA
#' for each region; each object of the list is named after the
#' MultiAssayExperiment object from which it derives
#' @import MultiAssayExperiment BiocParallel
#' @export
#'
#' @examples
#' if (interactive()) {
#' data(exprMatrix)
#' data(trainingCovariates)
#' data(regionAssoc)
#'
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#'
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, affiXcanTraining=training,
#' scale=TRUE)
#' }
affiXcanPcs <- function(tbaPaths, affiXcanTraining, scale, BPPARAM=bpparam(),
                        testingSamples=NULL) {

    pcs <- vector("list", affiXcanTraining$regionsCount)
    index <- 0

    for(i in seq(1,length(tbaPaths))) {

        pca <- affiXcanTraining$pca
        tbaMatrixMAE <- readRDS(tbaPaths[i])
        
        if (is.null(testingSamples)==FALSE) {
            tbaMatrixMAE <- MultiAssayExperiment::subsetByRow(tbaMatrixMAE,
                                                              testingSamples)
        }
        
        tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)
        regionsList <- setNames(as.list(names(tbaMatrix)), names(tbaMatrix))
        
        newPcs <- BiocParallel::bplapply(X=regionsList, FUN=computePcs,
            tbaMatrix, scale, pca, BPPARAM=BPPARAM)
        
        for(l in seq(1,length(newPcs))) {
            pcs[index + l] <- newPcs[l]
            names(pcs)[index + l] <- names(newPcs)[l]
        }

        index <- index + length(newPcs)
    
    }

    return(pcs)
}

#' Compute a matrix product between variables and eigenvectors
#'
#' @param region A string, which is the name of the object in the list
#' MultiAssayExperiment::experiments(tbaMatrix) that contains the TBA values of
#' the genomic region of interest (see the param tbaMatrix)
#' @param tbaMatrix A MultiAssayExperiment object containing the TBA values
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param pca The returningObject$pca from affiXcanTrain()
#'
#' @return A data.frame containing the principal components values of the TBA
#' in a certain genomic region, as the result of the matrix product between the
#' TBA values and the eigenvectors
#' @import MultiAssayExperiment
#' @export
#'
#' @examples
#' if (interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#'
#' region <- "ENSG00000256377.1"
#'
#' tbaMatrixMAE <- readRDS(system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan"))
#'
#' pca <- training$pca
#' pcs <- computePcs(region=region, tbaMatrix=tbaMatrixMAE, scale=TRUE, pca=pca)
#' }
computePcs <- function(region, tbaMatrix, scale, pca) {

    myRegion <- region[[1]]
    myTbaMatrix <- tbaMatrix[[myRegion]]

    tryCatch (
        {
            myEigenvectors <- as.matrix(pca[[myRegion]]$eigenvectors)
            dataScaled <- as.matrix(scale(myTbaMatrix, center = TRUE,
                scale = scale))
            myPcs <- as.data.frame(dataScaled %*% myEigenvectors)
            return(pcs=myPcs)
        },
        error = function(e) {
            myPcs <- NA
            return(pcs=myPcs)
        }
    )
}

#' Compute a GReX from variables and their coefficients for a set of genes
#'
#' @param affiXcanTraining The returning object from affiXcanTrain()
#' @param pcs A list, which is the returning object from affiXcanPcs()
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#'
#' @return A SummarizedExperiment object containing the imputed GReX values
#' @import SummarizedExperiment BiocParallel 
#' @export
#'
#' @examples
#' if (interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#'
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, affiXcanTraining=training,
#' scale=TRUE)
#'
#' exprmatrix <- affiXcanGReX(affiXcanTraining=training, pcs=pcs)
#' }
affiXcanGReX <- function(affiXcanTraining, pcs, BPPARAM=bpparam()) {
    bs <- affiXcanTraining$bs
    
    expr <- BiocParallel::bplapply(X=bs, FUN=computeExpr, pcs, BPPARAM=BPPARAM)

    exprmatrix <- matrix(nrow = nrow(pcs[[1]]), ncol=0)
    rownames(exprmatrix) <- rownames(pcs[[1]])
    
    for(i in seq(1:length(expr))) {
        if(is.na(expr[[i]][[1]])==FALSE) {
            exprMatrixPart <- matrix(expr[[i]], ncol=1)
            colnames(exprMatrixPart)[[1]] <- ls(expr[i])
            exprMatrixPart <- as.matrix(exprMatrixPart, ncol=1)
            exprmatrix <- cbind(exprmatrix, exprMatrixPart)
        }
    }
    
    ## TODO: possibly include a colData given by the user in the SE object
    Summarized.GReX <- SummarizedExperiment(assays=list(GReX=t(exprmatrix)))
    
    return(Summarized.GReX)
}

#' Compute the imputed GReX for a certain gene on a set of individuals
#'
#' @param bs A list containing three objects:
#' \itemize{
#'  \item coefficients: An object containing the coefficients of the principal
#'  components used in the model, completely similar to the "coefficients"
#'  object from the results of lm()
#'  \item pval: The uncorrected anova pvalue of the model, retrieved from
#'  anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'  If cov==NULL, i.e. no covariates for the population structure have
#'  been provided, pval will be NA
#'  \item r.sq: The coefficient of determination between the real total expression
#'  values and the imputed GReX, retrived from summary(model)$r.squared
#' }
#' @param pcs A list, which is the returning object of affiXcanPcs()
#'
#' @return A vector of imputed GReX values
#' @export
#'
#' @examples
#' if (interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, affiXcanTraining=training,
#' scale=TRUE)
#'
#' region <- "ENSG00000256377.1"
#' bs <- training$bs[[region]]
#'
#' exprmatrix <- computeExpr(bs=bs, pcs=pcs)
#' }
computeExpr <- function(bs, pcs) {

    tryCatch (
        {   
            ## ANOVA p-value is not used anymore (changes in version 1.3.1)
            #if(is.na(bs$correctedP) == TRUE || bs$correctedP > 0.05) {
            #    exprVector <- NA
            #    return(exprVector)
            #} else {
            
            ## create vector with all regulatory region names associated
            ## with the expressed region
            regulRegions <- names(bs$coefficients[-1])
            for(i in seq(1:length(regulRegions))) {

                name <- strsplit(regulRegions[[i]], split="@", fixed=TRUE)
                regulRegions[[i]] <- name[[1]][[2]]

            }
            regulRegions <- unique(regulRegions)
            ## create data frame with all the predictors 
            pcVals <- as.data.frame(matrix(nrow = nrow(pcs[[regulRegions]]),
                ncol=0))
            rownames(pcVals) <- as.vector(rownames(pcs[[regulRegions]]))
            for(region in regulRegions) {
                pcValsPart <- pcs[[region]]
                colnames(pcValsPart) <- paste(colnames(pcValsPart), "@",
                    region, sep="")
                pcVals <- cbind(pcVals, pcValsPart)
            }
            ## create vector with imputed expression values
            exprVector <- vector(mode="numeric", length=nrow(pcVals))
            intercept <- bs$coefficients[[1]]
            coefficients <- as.vector(bs$coefficients[-1])
            for(j in seq(1, nrow(pcVals))) {
                expr <- sum(coefficients * as.vector(pcVals[j, ]))
                    + intercept
                exprVector[j] <- expr
            }

            return(exprVector) 
        },
        error = function(e) {
            exprVector <- NA
            return(exprVector)
        }
    )
}

#' Compute R and R^2 between rows of two SummarizedExperiment assays
#'
#' @param realExpr A SummarizedExperiment object containing expression data
#' @param imputedExpr The returning object of affiXcanImpute()
#' @param assay A string with the name of the object in
#' SummarizedExperiment::assays(realExpr) that contains expression values
#' @param testingSamples A vector of strings. The identifiers of the samples
#' that have to be considered by the function; default is NULL; if 
#' is.null(testingSamples)==TRUE then no filtering is performed
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#'
#' @return A list of lists; inner lists are named after the rows for which the
#' correlation between realExpr and imputedExpr have been computed; inner
#' lists contain two objects:
#' \itemize{     
#'  \item pearson: the pearson's correlation coefficient (R) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed with cor()
#'  \item rSquared: the coefficient of determination (R^2) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed as pearson^2
#' }
#' @import SummarizedExperiment BiocParallel
#' @export
#'
#' @examples
#' if (interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#' 
#' assay <- "values"
#' 
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#' 
#' imputedExpr <- affiXcanImpute(tbaPaths=trainingTbaPaths,
#' affiXcanTraining=training, scale=TRUE)
#' realExpr <- exprMatrix
#' 
#' correlation <- computeRSquared(realExpr, imputedExpr, assay)
#' }
computeRSquared <- function(realExpr, imputedExpr, assay,
                            testingSamples=NULL, BPPARAM=bpparam()) {
    geneNames <- rownames(imputedExpr)
    names(geneNames) <- geneNames

    if (is.null(testingSamples)==FALSE) {
        realExpr <- realExpr[,colnames(exprMatrix) %in% testingSamples]
    }
   
    ## check/correct the consistency of the order of the samples
    if(is.na(table(colnames(imputedExpr)==colnames(realExpr))['FALSE'])==FALSE){
        imputedExpr <- imputedExpr[,colnames(realExpr)]
    }
    
    imputedExpr <- SummarizedExperiment::assays(imputedExpr)$GReX
    realExpr <- SummarizedExperiment::assays(realExpr)[[assay]] 
    correlation <- BiocParallel::bplapply(X=geneNames, FUN=computeCorrelation,
                                       realExpr, imputedExpr, BPPARAM=BPPARAM)
    return(correlation)
}

#' Compute R and R^2 on a particular row of two SummarizedExperiment assays
#'
#' @param geneName A string. The row name in realExpr and imputedExpr objects
#' that identifies the vectors between which R and R^2 have to be computed
#' @param realExpr A SummarizedExperiment object containing expression data
#' @param imputedExpr The returning object of affiXcanImpute()
#'
#' @return A list of two objects:
#' \itemize{     
#'  \item pearson: the pearson's correlation coefficient (R) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed with cor()
#'  \item rSquared: the coefficient of determination (R^2) between the real
#'  expression values and the imputed GReX for the cross-validation i on
#'  the expressed gene y, computed as pearson^2
#' }
#' @export
#'
#' @examples
#' if (interactive()) {
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#' 
#' assay <- "values"
#' 
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE)
#' 
#' imputedExpr <- affiXcanImpute(tbaPaths=trainingTbaPaths,
#' affiXcanTraining=training, scale=TRUE)
#' realExpr <- exprMatrix
#' 
#' geneName <- "ENSG00000256377.1"
#' imputedExpr <- SummarizedExperiment::assays(imputedExpr)$GReX
#' realExpr <- SummarizedExperiment::assays(realExpr)[[assay]]
#'
#' correlation <- computeCorrelation(geneName, realExpr, imputedExpr) 
#' }
computeCorrelation <- function(geneName, realExpr, imputedExpr) {
    realExpr <- realExpr[rownames(realExpr)==geneName]
    imputedExpr <- imputedExpr[rownames(imputedExpr)==geneName]
    pearson <- cor(imputedExpr, realExpr)
    rSquared <- pearson^2
    correlation <- list("pearson"=pearson, "rSquared"=rSquared)
    return(correlation)
}

#' Impute a GReX for each gene for which a model was generated
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param affiXcanTraining The returning object from affiXcanTrain()
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA; default is TRUE
#' @param BPPARAM A BiocParallelParam object. Default is bpparam(). For
#' details on BiocParallelParam virtual base class see 
#' browseVignettes("BiocParallel")
#'
#' @return A SummarizedExperiment object containing imputed GReX values
#' @import MultiAssayExperiment SummarizedExperiment BiocParallel crayon
#' @export
#'
#' @examples
#' trainingTbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' data(exprMatrix)
#' data(regionAssoc)
#' data(trainingCovariates)
#'
#' assay <- "values"
#'
#' training <- affiXcanTrain(exprMatrix=exprMatrix, assay=assay,
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc,
#' varExplained=80, scale=TRUE)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' exprmatrix <- affiXcanImpute(tbaPaths=testingTbaPaths,
#' affiXcanTraining=training, scale=TRUE)
affiXcanImpute <- function(tbaPaths, affiXcanTraining, scale=TRUE,
                           BPPARAM=bpparam()) {
    regionsCount <- overlookRegions(tbaPaths)
    if(regionsCount!=affiXcanTraining$regionsCount) {
        warning(cat("The amount of genomic regions included in the training
                    phase is different from those found in the TBA matrices:\n
                    TBA in the training set was computed on ",
                    affiXcanTraining$regionsCount, " regions\n TBA provided in
                    tbaPaths refers to ", regionsCount, " regions\n"))
    }
    cat(bold(green("\nAffiXcan")),": ",
        cyan("Imputing Genetically Regulated Expression (GReX)\n"), sep="")
    cat(bold(cyan("\t--> ")), green("Computing Principal Components\n"), sep="")
    pcs <- affiXcanPcs(tbaPaths, affiXcanTraining, scale, BPPARAM)
    cat(bold(cyan("\t--> ")), green("Imputing GReX values\n"), sep="")
    exprmatrix <- affiXcanGReX(affiXcanTraining, pcs, BPPARAM)
    cat(bold(cyan("\tDone!\n")))
    return(exprmatrix)
}
