## main function for estimating linear mixed models via REHE
## Code is adapted from R package GENESIS (v2.14.3)

## For faster computation, set VConly to TRUE if only computing point estimations for variance components (OLS estimation for beta will be computed)
## If subsequent GWAS analysis will be performed, need to set VConly=F to obtain the required intermediate quantities (WLS estimation for beta will be computed)

## computeCI=T if would like to compute confidence intervals for variance components. Fixed effects will use weighted lease square.

fitREHE = function(y, X, K, group.idx=NULL, computeCI = F, VConly = F){ 
  
  # y is the response vector, X is the fixed covariate design matrix (the first column being 1), K is the list of correlation matrices
  
  
  if (is.null(group.idx)){group.idx <- list(resid.var = 1:length(y))}else{stop("Heterogeneous residual variance feature is not developed.")}
  
  XTXXT = solve(t(X)%*%X, t(X))
  beta = XTXXT%*%y
  # P = X%*%XTXXT
  n = length(y)
  In = diag(1, n)
  nK = length(K)
  
  residual = y - X%*%beta

  y_matrix = as.vector(tcrossprod(residual))
  X_matrix = cbind(as.vector(In), sapply(1:nK, function(w)as.vector(K[[w]])))
  
  XTX = crossprod(X_matrix)
  XTY = crossprod(X_matrix, y_matrix)
  
  
  res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
  
  estimates = res_qp$solution
  estimates_he = res_qp$unconstrained.solution
  

  # numerical correction for REHE estimations (due to machine error)
  for(s in 1:length(estimates)){
    estimates[s] = ifelse(estimates[s]<0, 0, estimates[s])}
  
  estimates = c(estimates[2:(nK+1)], residual=estimates[1])
  names(estimates) = c(paste('K', 1:nK), 'residual')
  
  if(!VConly || computeCI){ 
    # functions to prepare compatible quantities for GWAS analysis under GENESIS package
    sq = .computeSigmaQuantities(varComp = estimates, covMatList = K, group.idx = group.idx)
    lq = .calcLikelihoodQuantities(Y = y, X = X, Sigma.inv = sq$Sigma.inv, cholSigma.diag = sq$cholSigma.diag)
    output = list(varComp = estimates, AI = NA, converged = T, zeroFLAG = 0, niter = 0,
                  Sigma.inv = sq$Sigma.inv, beta = lq$beta, residM = lq$residM, fits = lq$fits, eta = 0, 
                  logLikR = lq$logLikR, logLik = lq$logLik, RSS = lq$RSS, In=In, cholSigma = sq$cholSigma)
  }else{ 
    sq=list(NULL)
    lq=list(NULL)
    output = list(varComp = estimates, beta = beta)
  }
  
  
  return(output)
  
  
}




