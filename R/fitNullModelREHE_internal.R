### Non-gaussian families were not supported. 
### Sparse correlation matrices were regarded as 'dense' in the computation; possible to adapt the code for faster computation by using sparse matrices characteristics
### The following comments from GENESIS developer were kept for reference

    ### A function that fits a null model, assuming complete data - will have to be wrapped using a different function that performs checks. 
    ## or checks will be added... for a start - complete valid data. nrow(X) == length(y), dimensions of covMatList similarly appropriate. We
    ## do not deal with IDs, just indices, as everything is assumed to match. 
    ## X is assumed to have an intercept. 

    ## y - outcome vector
    ## X - data.frame or model.matrix

    ## computeCI - if Ture, will compute confidence interval for variance component estimates via bootstrap; will be slow for large samples
    ## VConly - if True, will only focus on computing variance component point estimations (will not be able to use fitted model for subsequent GWAS analysis). Fast.

.fitNullModelREHE <- function(y, X, covMatList = NULL, group.idx = NULL, family = "gaussian", start = NULL,
                          AIREML.tol = 1e-4, max.iter = 100, EM.iter = 0,
                          drop.zeros = TRUE, verbose = TRUE, computeCI=F){
    
    ### checks
    if(!is.null(covMatList)){
        if (!is.list(covMatList)){
            covMatList <- list(A = covMatList)
        }
        # if any Matrix objects; coerce all to Matrix objects (coerced ones are not sparse)
        covMatList <- .checkMatrixType(covMatList)
    }

    if (is.null(colnames(X))){
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }
    if (!is.element(family$family, c("gaussian"))){
        stop("family must be gaussian for estimation via REHE")
    }

    ### Gaussian family
  
    if (is.null(covMatList) & is.null(group.idx)) {
      # linear regression
      mod <- lm(y ~ -1 + X)
      out <- .nullModOutReg(y, X, mod, family)
    }
  
    if (is.null(covMatList) & !is.null(group.idx)){
      vc.mod <- .runWLSgaussian(y, X, group.idx = group.idx, start = start, 
                                AIREML.tol = AIREML.tol, max.iter = max.iter, 
                                EM.iter = EM.iter, verbose = verbose)
      out <- .nullModOutWLS(y, X, vc.mod = vc.mod, family = family, group.idx = group.idx)
    }
    if (!is.null(covMatList)){
      if (is.null(group.idx)) group.idx <- list(resid.var = 1:length(y))
        
      vc.mod <- fitREHE(y, X, K = covMatList, group.idx = group.idx, computeCI = computeCI)
    

      # use paramtric bootstrap to get covariance/inverse covariance of VCcomp, 
      nK = length(covMatList)
      n = length(y)
    
      # compute confidence interval for VC
      varCI = NA
    
      if(computeCI){
        cat('compute VC confidence interval \n')
        varCI = REHE_CI(n, K=covMatList, vc_rehe_full = c(vc.mod$varComp[nK+1], vc.mod$varComp[-c(nK+1)]), In=vc.mod$In, cholSigma = vc.mod$cholSigma)
        cat('compute CI finished \n')
        
      }
    
    
      out <- .nullModOutMM(y = y, workingY = y, X = X, vc.mod = vc.mod, 
                           family = family, covMatList = covMatList, 
                           group.idx = group.idx, drop.zeros = drop.zeros)
    
      out$varCompCI <- varCI
    }
  

    out.class <- class(out)
    nullprep <- nullModelTestPrep(out)
    out <- c(out, nullprep)
    class(out) <- out.class
    
    return(out)
    
}





