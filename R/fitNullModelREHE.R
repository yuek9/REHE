#############################################
##
##      Fit null model for GWAS analysis via REHE
## 
##      Note: Currently only support Gaussian family (linear mixed models)
##
##
##      Code built based on R package GENESIS (v2.14.3) 
##      Reference: Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using the GENESIS R/Bioconductor package. Bioinformatics. doi:10.1093/bioinformatics/btz567.
## 
##      Author:           Kun Yue
##      Date modified:    02/21/2020
#############################################

## required library: GENESIS(v2.14.3 works), quadprog (v1.5-7 works) , 

setGeneric("fitNullModelREHE", function(x, ...) standardGeneric("fitNullModelREHE"))

setMethod("fitNullModelREHE",
          "data.frame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   family = "gaussian",
                   start = NULL,
                   AIREML.tol = 1e-4,
                   max.iter = 100,
                   EM.iter = 0,
                   drop.zeros = TRUE,
                   verbose = TRUE, computeCI=F) {
            
              if(!is.null(group.var)){
                stop('Heterogeneous residual variance currently not supported. Please leave group.var as NULL')
              }
              
              desmat <- createDesignMatrix(x, outcome, covars, group.var)

              # if there was missing data, need to subset cov.mat
              if (!is.null(cov.mat)) {
                  .checkRownames(cov.mat, x)
                  if (nrow(desmat$X) < nrow(x)) {
                      ind <- which(rownames(x) %in% rownames(desmat$X))
                      cov.mat <- .covMatSubset(cov.mat, ind)
                  }
              }
              
              .fitNullModelREHE(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                            group.idx=desmat$group.idx, family=family,
                            start=start, AIREML.tol=AIREML.tol,
                            max.iter=max.iter, EM.iter=EM.iter,
                            drop.zeros=drop.zeros, verbose=verbose, computeCI=computeCI)
          })

setMethod("fitNullModelREHE",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   ...) {
              
              x <- pData(x)
              rownames(x) <- x$sample.id

              if (!is.null(cov.mat)) {
                  .checkSampleId(cov.mat, x)
              }
              
              ## subset data.frame and cov.mat for selected samples
              if (!is.null(sample.id)) {
                  stopifnot(all(sample.id %in% x$sample.id))
                  ind <- x$sample.id %in% sample.id
                  x <- x[ind,]
                  if (!is.null(cov.mat)) {
                      ind <- which(.covMatNames(cov.mat) %in% sample.id)
                      cov.mat <- .covMatSubset(cov.mat, ind)
                  }
              }

              ## reorder data.frame to match cov.mat
              if (!is.null(cov.mat)) {
                  ind <- match(.covMatNames(cov.mat), rownames(x))
                  x <- x[ind,]
              }
              
              nm <- fitNullModelREHE(x, outcome, covars, cov.mat, group.var, ...)
              nm$sample.id <- rownames(nm$model.matrix)
              nm
          })

setMethod("fitNullModelREHE",
          "SeqVarData",
          function(x, ...) {
              fitNullModelREHE(sampleData(x), ...)
          })

setMethod("fitNullModelREHE",
          "ScanAnnotationDataFrame",
          function(x, ...) {
              class(x) <- "AnnotatedDataFrame"
              varLabels(x)[varLabels(x) == "scanID"] <- "sample.id"
              fitNullModelREHE(x, ...)
          })

setMethod("fitNullModelREHE",
          "GenotypeData",
          function(x, ...) {
              fitNullModelREHE(getScanAnnotation(x), ...)
          })


