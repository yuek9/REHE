fitREHE = function(y, X, K, group.idx,
                   subsample=F, ratio=0.05, nsub = 50, option='mean'){
  
  # y is the response, X is the fixed covariate design matrix
  # need to regress out the fixed effects in advance
  # use LS estimation for beta
  XTXXT = solve(t(X)%*%X, t(X))
  beta = XTXXT%*%y
  P = X%*%XTXXT
  n = length(y)
  In = diag(1, n)
  nK = length(K)
  
  
  cat('projection of covariance calculation (do not project in this implementation) \n')
  # K[[1]] = (In-P)%*%K[[1]]%*%(In-P)   ## this projection step is very time consuming
  cat('projection Done \n')
  
  residual = y - X%*%beta

  if(!subsample){
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
    
    estimates = c(estimates[2:(nK+1)], estimates[1]) # need to put the residual at the last
    
  }else{
    
    yyT = tcrossprod(residual)
    
    reREHE_res = replicate(n=nsub, {
      
      # for one subsample
      sub_index = sample(n, size=round(ratio*n), replace=T) # sample the responses
      
      # system.time({
      #   index_list = as.matrix(expand.grid(sub_index, sub_index))
      # 
      # y_matrix = as.vector(residual[sub_index] %*% t(residual[sub_index]))
      # X_matrix = cbind(In[index_list], sapply(1:nK, function(w)K[[w]][index_list]))
      # })
      
      
      y_matrix = as.vector(yyT[sub_index, sub_index])
      X_matrix = cbind(as.vector(In[sub_index, sub_index]), sapply(1:nK, function(w)as.vector(K[[w]][sub_index, sub_index])))
      
      
      XTX = crossprod(X_matrix)
      XTY = crossprod(X_matrix, y_matrix)
      
#      print(X_matrix)
      
      res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
      
      estimates = res_qp$solution
      estimates_he = res_qp$unconstrained.solution
      
      # non-negative correction for HE estimations
      for(s in 1:length(estimates_he)){
        estimates_he[s] = ifelse(estimates_he[s]<0, 0, estimates_he[s])}
      
      # numerical correction for REHE estimations (due to machine error)
      for(s in 1:length(estimates)){
        estimates[s] = ifelse(estimates[s]<0, 0, estimates[s])}
      
      return(estimates)
    })
    
#    print(reREHE_res)
    if(option=='mean'){
      estimates = apply(reREHE_res, 1, mean)
    }
    
    if(option=='median'){
      estimates = apply(reREHE_res, 1, median)
    }
    
    estimates = c(estimates[2:(nK+1)], estimates[1]) # need to put the residual at the last
    
  }
  

  
  sq = .computeSigmaQuantities(varComp = estimates, covMatList = K, group.idx = group.idx)
  lq = .calcLikelihoodQuantities(Y = y, X = X, Sigma.inv = sq$Sigma.inv, cholSigma.diag = sq$cholSigma.diag)
  
  
  
  return(list(varComp = estimates, AI = 0, converged = T, zeroFLAG = 0, niter = 0,
       Sigma.inv = sq$Sigma.inv, beta = lq$beta, residM = lq$residM, fits = lq$fits, eta = 0, 
       logLikR = lq$logLikR, logLik = lq$logLik, RSS = lq$RSS, In=In, cholSigma = sq$cholSigma))
  
  
}




