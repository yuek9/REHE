### A function that fits a null model, assuming complete data - will have to be wrapped using a different function that performs checks. 
## or checks will be added... for a start - complete valid data. nrow(X) == length(y), dimensions of covMatList similarly appropriate. We
## do not deal with IDs, just indices, as everything is assumed to match. 
## X is assumed to have an intercept. 
## non-gaussian families can only be binomial and poisson. 

## y - outcome vector
## X - data.frame or model.matrix
.fitNullModelREHE <- function(y, X, covMatList = NULL, group.idx = NULL, family = "gaussian", start = NULL,
                          AIREML.tol = 1e-4, max.iter = 100, EM.iter = 0,
                          drop.zeros = TRUE, verbose = TRUE, computeCI=F,
                          subsample=F, ratio=0.05, nsub = 50, option='mean'
                          ){
    
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
    if (!is.element(family$family, c("gaussian", "binomial", "poisson"))){
        stop("family must be one of gaussian, binomial, or poisson")
    }

    ### Gaussian family
    if (family$family == "gaussian"){
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
            # vc.mod <- .runAIREMLgaussian(y, X, start = start, covMatList = covMatList, 
            #                              group.idx = group.idx, AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,  
            #                              max.iter = max.iter, EM.iter = EM.iter, verbose = verbose)
            
            # need to modify to incoporate multiple variance components
            vc.mod <- fitREHE(y, X, K = covMatList, group.idx = group.idx, subsample=subsample, ratio=ratio, nsub = nsub, option=option)
            
            
            
            # use paramtric bootstrap to get covariance/inverse covariance of VCcomp, 
            # here can be 'exact' since we have the decomposition
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
    } 

    ### Non-Gaussian family
    if (family$family != "gaussian"){
        # initial fit with glm
        mod <- glm(y ~ X, family = family)
        
        if (!is.null(covMatList)){ ## iterate between computing workingY and estimating VCs. 
            iterate.out <- .iterateAIREMLworkingY(glm.mod = mod, X = X, family = family, 
                                                  start = start, covMatList = covMatList, AIREML.tol = AIREML.tol,
                                                  drop.zeros = drop.zeros, max.iter = max.iter, EM.iter = EM.iter,
                                                  verbose = verbose)
            
            vc.mod <- iterate.out$vc.mod
            working.y <- iterate.out$working.y
            
      	    ## check whether all variance components were estimated as zero:
            if (vc.mod$allZero == TRUE){
                out <- .nullModOutReg(y, X, mod, family)
                out$zeroFLAG <- TRUE
            } else{
                out <- .nullModOutMM(y = y, workingY = working.y$Y, X = X, vc.mod = vc.mod, 
                                     family = family, covMatList = covMatList, 
                                     vmu = working.y$vmu, gmuinv = working.y$gmuinv, drop.zeros = drop.zeros)
            }	
        } else{
            out <- .nullModOutReg(y, X, mod, family)
        }
    }

    out.class <- class(out)
    nullprep <- nullModelTestPrep(out)
    out <- c(out, nullprep)
    class(out) <- out.class
    
    return(out)
    
}



# update this version to compute CI based on several covariance matrices
REHE_CI = function(n, K, # provide sparse Kin_sparse if do not provide cholesky decomposition; or can provide dense K with cholSigma
                   vc_rehe_full, In, cholSigma=NULL){
    
    # vc_rehe_full should be in order of c(e2, k2...)
    # output confidence interval for both sd based CI (wald type) and empirical CI (percentile of difference)
    
    nK = length(K)

    # get the sd estimation from 100 replications; could also use quantile based
    n_boots = 100
    
    
    if(is.null(cholSigma)){ # will report issue if covariance matrix cannot have chol; for now do not really this this step
    
    tmp = lapply(K,Matrix::chol) # this is fast, gives R, where R'R = Sigma
    
    sqrt.Sigma_approx_Kin = Reduce('+', lapply(1:nK, function(i) t(as.matrix(tmp[[i]]))*sqrt(vc_rehe_full[i+1])))
    
    tmp_sample = matrix(rnorm(2*n*n_boots), nrow=n, ncol=2*n_boots)
    
    residual_boots = sqrt.Sigma_approx_Kin %*% tmp_sample[,1:n_boots] + # this is adding Kin and In parts together
        sqrt(vc_rehe_full[1])* tmp_sample[,-c(1:n_boots)] #this is already ordered to Kin_sparse
    
    
    # not fully developed, do not use this for now
    if(class(Kin_sparse)=='dsCMatrix'){
        colindex = do.call(c,mapply(1:n, each=diff(Kin_sparse@p), FUN=rep)) # this will get nonzero entry column index
        rowindex = Kin_sparse@i+1 # this will give the nonzero row index
    }
    
    if(class(Kin_sparse)=='dsTMatrix'){
        colindex = Kin_sparse@j+1 # this will get nonzero entry column index
        rowindex = Kin_sparse@i+1 # this will give the nonzero row index
    }
    
    vec_K = Kin_sparse[cbind(rowindex, colindex)] # this is effective nonzero entries
    vec_I = In[cbind(rowindex, colindex)]

    
    vec_y_boots = residual_boots[rowindex,]*residual_boots[colindex,]

    
    }else{
        
    # if provided with cholesky decomposition of Sigma, directly use it for bootstrap
        sqrt.Sigma_approx = t(cholSigma)
        residual_boots = sqrt.Sigma_approx%*%matrix(rnorm(n*n_boots), nrow=n, ncol=n_boots) #this is already ordered to Kin_sparse
        
        all_K = list(In)
        all_K[2:(nK+1)] = K
    
        vec_K = list( as.vector(In))
        vec_K[2:4]= lapply(K, as.vector)
        #dim(residual_boots)
        
        # vec_y_boots = residual_boots[row(In),]*residual_boots[col(In),] # this is too large to store; do this one by one
        # list_y_boots = lapply(1:n_boots, function(j)tcrossprod(residual_boots[,j])) #same size
        
        #dim(vec_y_boots)
        XTY_boots = sapply(1:n_boots, function(j){  # this seems slower, but memory efficient
            tmp = tcrossprod(residual_boots[,j])
            sapply(1:(nK+1), function(i) sum(tmp*all_K[[i]]))
            })
        
        }
    
    # construct XTX, XTY matrix
    tmp = as.matrix(expand.grid(1:(nK+1), 1:(nK+1)))
    XTX_approx = matrix(sapply(1:nrow(tmp), function(l)sum(all_K[[tmp[l,1]]]*all_K[[tmp[l,2]]])), nK+1, nK+1)
    
    # XTY_boots = t(sapply(1:(nK+1), function(i)vec_K[[i]]%*%vec_y_boots))
    

    vc_boots = sapply(1:n_boots, function(i){res_qp_boots = solve.QP(Dmat = XTX_approx, dvec=XTY_boots[,i], Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
    res_qp_boots$solution})
    
    
    heri_boots = t(t(vc_boots)/colSums(vc_boots))
    
    heri_full = vc_rehe_full/sum(vc_rehe_full)
    
    # the following based on Wald type SD
    CIvc_wald = t(sapply(1:(nK+1), function(i) vc_rehe_full[i] + c(-1, 1)*qnorm(0.975)*sd(vc_boots[i,])))
    CIvc_quan = t(sapply(1:(nK+1), function(i) vc_rehe_full[i] + quantile(vc_boots[i,]-vc_rehe_full[i],c(0.025, 0.975))))
    
    CIheri_wald = t(sapply(1:(nK+1), function(i) heri_full[i] + c(-1, 1)*qnorm(0.975)*sd(heri_boots[i,])))
    CIheri_quan = t(sapply(1:(nK+1), function(i) heri_full[i] + quantile(heri_boots[i,]-heri_full[i],c(0.025, 0.975))))

        
    
    return(list(CIvc = list(wald = CIvc_wald, quantile = CIvc_quan),
                CIheri = list(wald = CIheri_wald, quantile = CIheri_quan)))  
    
    
}


