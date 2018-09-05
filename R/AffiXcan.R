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
#' @param cov A data.frame with covariates values for the population structure
#' where the columns are the PCs and the rows are the individual IIDs
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return A list containing two objects: pca and bs
#'
#' pca: A list containing lists named as the MultiAssayExperiment::experiments()
#' found in the MultiAssayExperiment objects listed in the param tbaPaths. Each
#' of these lists contains two objects:
#'
#'    eigenvectors: A matrix containing eigenvectors for those principal
#'    components selected according to the param varExplained
#'
#'    pcs: A matrix containing the principal components values selected
#'    according to the param varExplained
#'
#' bs: A list containing lists named as the REGULATORY_REGIONS found in the
#' param regionAssoc that have a correspondent colname in the experiments
#' stored in MultiAssayExperiment objects listed in the param tbaPaths.
#' Each of the lists in bs contains four objects:
#'
#'    coefficients: The coefficients of the principal components used in the
#'    model, completely similar to the "coefficients" from the results of lm()
#'
#'    pval: The uncorrected anova pvalue of the model, retrieved from
#'    anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'
#'    r.sq: The coefficient of determination between the real total expression
#'    values and the imputed GReX, retrived from summary(model)$r.squared
#'
#'    correctedP: The p value after the benjamini-hochberg correction for
#'    multiple testing, retrived using p.adjust(pvalues, method="BH")
#' @import MultiAssayExperiment SummarizedExperiment doParallel plyr
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
affiXcanTrain <- function(exprMatrix, assay, tbaPaths, regionAssoc, cov,
    varExplained, scale, cores) {
    cat("AffiXcan: Performing PCA ...", "\n")
    pca <- affiXcanPca(tbaPaths, varExplained, scale, cores)
    cat("AffiXcan: Training coefficients ...", "\n")
    bs <- affiXcanBs(exprMatrix, assay, regionAssoc, pca, cov, cores)
    affiXcanTraining <- list(pca=pca, bs=bs)
    cat("AffiXcan: Done!", "\n")
    return(affiXcanTraining)
}

#' Impute a GReX for each gene for whom a model was generated
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param pca A list, which is the returningObject$pca from affiXcanTrain()
#' @param bs A list, which is the returningObject$bs from affiXcanTrain()
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return A SummarizedExperiment object containing imputed GReX values
#' @import MultiAssayExperiment SummarizedExperiment doParallel plyr
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' exprmatrix <- affiXcanImpute(tbaPaths=testingTbaPaths, pca=training$pca,
#' bs=training$bs, scale=TRUE, cores=1)
affiXcanImpute <- function(tbaPaths, pca, bs, scale, cores) {
    cat("AffiXcan: Computing principal components ...", "\n")
    pcs <- affiXcanPcs(tbaPaths, pca, scale, cores)
    cat("AffiXcan: Imputing GReX values ...", "\n")
    exprmatrix <- affiXcanGReX(bs, pcs, cores)
    cat("AffiXcan: Done!", "\n")
    return(exprmatrix)
}

#' Perform a PCA on each experiment found in MultiAssayExperiment objects
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return pca: A list containing lists named as the 
#' MultiAssayExperiment::experiments() found in the MultiAssayExperiment objects
#' listed in the param tbaPaths. Each of these lists contains two objects:
#'
#'    eigenvectors: A matrix containing eigenvectors for those principal
#'    components selected according to the param varExplained
#'
#'    pcs: A matrix containing the principal components values selected
#'    according to the param varExplained
#' @import MultiAssayExperiment doParallel plyr
#' @export
#'
#' @examples
#' tbaPaths <- system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan")
#'
#' pca <- affiXcanPca(tbaPaths=tbaPaths, varExplained=80, scale=TRUE, cores=1)
affiXcanPca <- function(tbaPaths, varExplained, scale, cores) {

    pca <- list()

    for(path in tbaPaths) {

        tbaMatrixMAE <- readRDS(path)
    	updateObject(tbaMatrixMAE)
        tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)
        rm(tbaMatrixMAE)
        gc()

        if (as.numeric(cores) > 1) {
            doParallel::registerDoParallel(cores)
            newPca <- plyr::llply(.data=tbaMatrix, .fun=computePca,
                varExplained, scale, .parallel = TRUE)
        } else {
            newPca <- plyr::llply(.data=tbaMatrix, .fun=computePca,
                varExplained, scale, .parallel = FALSE)
        }

        pca <- append(pca, newPca)
        rm(tbaMatrix)
        gc()
    }

    return(pca)

}

#' Compute PCs in MultiAssayExperiment objects using eigenvectors given by user
#'
#' @param tbaPaths A vector of strings, which are the paths to
#' MultiAssayExperiment RDS files containing the tba values
#' @param pca A list, which is the returningObject$pca from affiXcanTrain()
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return A list of matrices containing the principal components values of TBA
#' for each region; each object of the list is named after the
#' MultiAssayExperiment object from which it derives
#' @import MultiAssayExperiment doParallel plyr
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, pca=training$pca, scale=TRUE,
#' cores=1)
affiXcanPcs <- function(tbaPaths, pca, scale, cores) {

    pcs <- list()

    for(path in tbaPaths) {

        tbaMatrixMAE <- readRDS(path)
    	updateObject(tbaMatrixMAE)
        tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)
        regionsList <- setNames(as.list(names(tbaMatrix)), names(tbaMatrix))
        rm(tbaMatrixMAE)
        gc()

        if (as.numeric(cores) > 1) {
            doParallel::registerDoParallel(cores)
            newPcs <- plyr::llply(.data=regionsList, .fun=computePcs, tbaMatrix,
                scale, pca, .parallel = TRUE)
        } else {
            newPcs <- plyr::llply(.data=regionsList, .fun=computePcs, tbaMatrix,
                scale, pca, .parallel = FALSE)
        }

        pcs <- append(pcs, newPcs)
        rm(tbaMatrix)
        rm(regionsList)
        gc()
    }

    return(pcs)
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
#' @param cov A data.frame with covariates values for the population structure
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return A list containing lists named as the REGULATORY_REGIONS found in the
#' param regionAssoc that have a correspondent name in the param pca.
#' Each of these lists contains three objects:
#'
#'    coefficients: An object containing the coefficients of the principal
#'    components used in the model, completely similar to the "coefficients"
#'    from the results of lm()
#'
#'    pval: The uncorrected anova pvalue of the model, retrieved from
#'    anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'
#'    r.sq: The coefficient of determination between the real total expression
#'    values and the imputed GReX, retrived from summary(model)$r.squared
#'
#'    correctedP: The p value after the benjamini-hochberg correction for
#'    multiple testing, retrived using p.adjust(pvalues, method="BH")
#' @import SummarizedExperiment doParallel plyr
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
#'
#' pca <- training$pca 
#'
#' bs <- affiXcanBs(exprMatrix=exprMatrix, assay=assay, regionAssoc=regionAssoc,
#' pca=pca, cov=trainingCovariates, cores=1)
affiXcanBs <- function(exprMatrix, assay, regionAssoc, pca, cov, cores) {
    expr <- SummarizedExperiment::assays(exprMatrix)[[assay]]
    tDfExprMatrix<-t(as.data.frame(expr))
    rm(expr)
    gc()
    if (as.numeric(cores) > 1) {
        doParallel::registerDoParallel(cores)
        bs <- dlply(.data=regionAssoc, .(EXPRESSED_REGION), .fun=computeBs, pca,
            tDfExprMatrix, cov, .parallel = TRUE)
    } else {
        bs <- dlply(.data=regionAssoc, .(EXPRESSED_REGION), .fun=computeBs, pca,
            tDfExprMatrix, cov, .parallel = FALSE)
    }

    ## Remove ensg for which the model is NA
    for(i in names(bs)) {
        if(is.na(bs[[i]]$r.sq)) {
            bs[[i]] <- NULL
        }
    }

    ## Compute multitest pvalue correction (benjamini-hochberg)
    pvalues <- vector("numeric", length=length(bs))
    for(i in seq(1:length(pvalues))) {
        pvalues[[i]] <- bs[[i]]$pval
    }
    multi.test.pvals <- p.adjust(pvalues, method="BH")
    for(i in seq(1:length(bs))) {
        bs[[i]]$correctedP <- multi.test.pvals[[i]]
    }
    rm(exprMatrix)
    rm(regionAssoc)
    gc()
    return(bs)
}

#' Compute a GReX from variables and their coefficients for a set of genes
#'
#' @param bs A list. which is the returningObject$bs from affiXcanTrain()
#' @param pcs A list, which is the returningObject from affiXcanPcs()
#' @param cores An integer >0; if cores=1 processes will not be parallelized
#'
#' @return A SummarizedExperiment object containing the imputed GReX values
#' @import SummarizedExperiment doParallel plyr
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#'
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, pca=training$pca, scale=TRUE,
#' cores=2)
#'
#' bs <- training$bs
#'
#' exprmatrix <- affiXcanGReX(bs=bs, pcs=pcs, cores=1)
affiXcanGReX <- function(bs, pcs, cores) {
    if (as.numeric(cores > 1)) {
        doParallel::registerDoParallel(cores)
        expr <- plyr::llply(.data=bs, .fun=computeExpr, pcs, .parallel = TRUE)
    } else {
        expr <- plyr::llply(.data=bs, .fun=computeExpr, pcs, .parallel = FALSE)
    }

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
    ## TODO: eventually include a colData given by the user in the SE object
    Summarized.GReX <- SummarizedExperiment(assays=list(GReX=t(exprmatrix)))
    rm(exprmatrix)
    gc()
    return(Summarized.GReX)
}

#' Perform a PCA on a matrix where columns are variables
#'
#' @param data A matrix containing the TBA values for a certain genomic region;
#' columns are PWMs, rows are individuals IIDs
#' @param varExplained An integer between 0 and 100; varExplained=80 means that
#' the principal components selected to fit the models must explain at least
#' 80 percent of variation of TBA values
#' @param scale A logical; if scale=FALSE the TBA values will be only centered,
#' not scaled before performing PCA
#'
#' @return A list containing two objects:
#' 
#'     eigenvectors: a matrix containing eigenvectors for those principal
#'     components selected according to the param varExplained
#' 
#'     pcs: a matrix containing the principal components values selected
#'     according to the param varExplained
#' @export
#'
#' @examples 
#' tbaMatrixMAE <- readRDS(system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan"))
#'
#' tbaMatrix <- MultiAssayExperiment::experiments(tbaMatrixMAE)
#' 
#' tba <- tbaMatrix$ENSG00000256377.1
#'
#' pca <- computePca(data=tba, varExplained=80, scale=TRUE)
computePca <- function(data, varExplained, scale) {

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
#' varExplained=80, scale=TRUE, cores=1)
#'
#' region <- "ENSG00000256377.1"
#'
#' tbaMatrixMAE <- readRDS(system.file("extdata","training.tba.toydata.rds",
#' package="AffiXcan"))
#'
#' pca <- training$pca
#' pcs <- computePcs(region=region, tbaMatrix=tbaMatrixMAE, scale=TRUE, pca=pca)
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

#' Fit a linear model to impute a GReX for a certain gene
#'
#' @param assocRegions A data.frame with the associations between regulatory
#' regions and one expressed gene, and with colnames = c("REGULATORY_REGION",
#' "EXPRESSED_REGION")
#' @param pca The returningObject$pca from affiXcanTrain()
#' @param expr A matrix containing the real total expression values, where the
#' columns are genes and the rows are individual IIDs
#' @param covariates A data.frame with covariates values for the population
#' structure where the columns are the PCs and the rows are the individual IIDs
#'
#' @return A list containing three objects:
#'
#'    coefficients: An object containing the coefficients of the principal
#'    components used in the model, completely similar to the "coefficients"
#'    from the results of lm()
#'
#'    pval: The uncorrected anova pvalue of the model, retrieved from
#'    anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'
#'    r.sq: The coefficient of determination between the real total expression
#'    values and the imputed GReX, retrived from summary(model)$r.squared
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
#' assocRegions <- regionAssoc[regionAssoc$EXPRESSED_REGION==
#' "ENSG00000256377.1",]
#'
#' pca <- affiXcanPca(tbaPaths=trainingTbaPaths, varExplained=80, scale=TRUE,
#' cores=1)
#'
#' expr <- SummarizedExperiment::assays(exprMatrix)[[assay]]
#' expr <- t(as.data.frame(expr))
#'
#' bs <- computeBs(assocRegions=assocRegions, pca=pca, expr=expr,
#' covariates=trainingCovariates)
computeBs <- function(assocRegions, pca, expr, covariates) {
    expressedRegion <- as.vector(unique(assocRegions$EXPRESSED_REGION))
    expression <- as.vector(as.data.frame(expr)[[expressedRegion]])
    tryCatch (
        {
            myPcs <- as.data.frame(matrix(nrow = nrow(covariates)))
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

            myPcsColnames <- c(myPcsColnames, names(covariates))
            myPcs <- cbind(myPcs, covariates)
            myPcs <- as.matrix(subset(myPcs, select=-V1))
            colnames(myPcs) = myPcsColnames

            model <- lm(expression~myPcs)
            modelReduced <- lm(expression~as.matrix(covariates))
            a <- anova(model, modelReduced, test="F")
            pval <- a$'Pr(>F)'[2]

            betas <- model$coefficients[seq(1,
                (length(model$coefficients)-length(names(covariates))))]
            r.sq <- summary(model)$r.squared
            results <- list(coefficients=betas, pval=pval, r.sq=r.sq)
            return(results)
        },
        error = function(e) {
            betas <- NA
            pval <- NA
            r.sq <- NA
            results <- list(coefficients=betas, pval=pval, r.sq=r.sq)
            return(results)
        }
    )
}

#' Compute the imputed GReX for a certain gene on a set of individuals
#'
#' @param bs A list containing three objects:
#'
#'    coefficients: An object containing the coefficients of the principal
#'    components used in the model, completely similar to the "coefficients"
#'    object from the results of lm()
#'
#'    pval: The uncorrected anova pvalue of the model, retrieved from
#'    anova(model, modelReduced, test="F")$'Pr(>F)'[2]
#'
#'    r.sq: The coefficient of determination between the real total expression
#'    values and the imputed GReX, retrived from summary(model)$r.squared
#' @param pcs A list, which is the returning object of affiXcanPcs()
#'
#' @return A summarizedExperiment object; 
#' SummarizedExperiment::assays(returningObject)$GReX is a matrix containing the
#' imputed GReX values
#' @import SummarizedExperiment
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
#' tbaPaths=trainingTbaPaths, regionAssoc=regionAssoc, cov=trainingCovariates,
#' varExplained=80, scale=TRUE, cores=1)
#'
#' testingTbaPaths <- system.file("extdata","testing.tba.toydata.rds",
#' package="AffiXcan")
#' 
#' pcs <- affiXcanPcs(tbaPaths=testingTbaPaths, pca=training$pca, scale=TRUE,
#' cores=1)
#'
#' bs <- training$bs
#'
#' exprmatrix <- computeExpr(bs=bs, pcs=pcs)
computeExpr <- function(bs, pcs) {

    tryCatch (
        {
            if(is.na(bs$correctedP) == TRUE || bs$correctedP > 0.05) {

                exprVector <- NA
                return(exprVector)

            } else {

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
            }
        },
        error = function(e) {
            exprVector <- NA
            return(exprVector)
        }
    )
}
