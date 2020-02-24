REHE_CI = function(n, 
                   K,              # provide sparse Kin_sparse if do not provide cholesky decomposition; or can provide dense K with cholSigma. 
                   vc_rehe_full,   # vc_rehe_full should be in order of c(residual, K[[1]], K[[2]]...)
                   In, 
                   cholSigma=NULL,
                   n_boots = 100   # by default use 100 bootstrap samples
){ 
  
  
  if(is.null(K)){
    stop('K missing. Need to provide a list of correlation matrix.')
  }
  
  if(class(K)!='list'){
    stop('K must be a list.')
  }
  # output confidence interval for both SD based CI (wald type) and quantile type CI
  
  nK = length(K)
  
  
  if(is.null(cholSigma)){ # if not provided with decomposition, need to compute it 
    
    Sigma = Reduce('+', lapply(1:nK, function(i) K[[i]]*(vc_rehe_full[i+1])))+In*vc_rehe_full[1]
    
    sqrt.Sigma_approx = Matrix::chol(Sigma) # this is fast, gives R, where R'R = Sigma
    
    tmp_sample = matrix(rnorm(2*n*n_boots), nrow=n, ncol=n_boots)
    
    residual_boots = sqrt.Sigma_approx %*% tmp_sample 
    
  }else{    # if provided with cholesky decomposition of Sigma, directly use it for bootstrap
    sqrt.Sigma_approx = t(cholSigma)
    residual_boots = sqrt.Sigma_approx%*%matrix(rnorm(n*n_boots), nrow=n, ncol=n_boots) #this is already ordered to Kin_sparse
  }
  
  
  if(nK==1 & class(K[[1]]) %in% c('dsCMatrix', 'dsTMatrix')){    # faster operations with only one sparse correlation matrix
    
    if(class(K[[1]])=='dsCMatrix'){
      colindex = do.call(c,mapply(1:n, each=diff(K[[1]]@p), FUN=rep)) # this will get nonzero entry column index
      rowindex = K[[1]]@i+1 # this will give the nonzero row index
    }
    
    if(class(K[[1]])=='dsTMatrix'){
      colindex = K[[1]]@j+1 # this will get nonzero entry column index
      rowindex = K[[1]]@i+1 # this will give the nonzero row index
    }
    
    vec_K = K[[1]][cbind(rowindex, colindex)] # this is effective nonzero entries
    vec_I = In[cbind(rowindex, colindex)]
    
    
    vec_y_boots = residual_boots[rowindex,]*residual_boots[colindex,]
    
    tmp = sum(vec_I*vec_K)
    XTX_approx = matrix(c(sum(vec_I), tmp, tmp, sum(vec_K^2)), 2, 2)
    
    XTY_boots = matrix(c(vec_I, vec_K), nrow=2, byrow=T)%*% vec_y_boots
    
  }else{                  # for more than one correlation matrices, or for dense matrices
    
    all_K = list(In)
    all_K[2:(nK+1)] = K
    
    #   vec_K= lapply(all_K, as.vector)
    #dim(residual_boots)
    
    # vec_y_boots = residual_boots[row(In),]*residual_boots[col(In),] # this is too large to store; do this one by one
    # list_y_boots = lapply(1:n_boots, function(j)tcrossprod(residual_boots[,j])) #same size
    
    
    XTY_boots = sapply(1:n_boots, function(j){  # this is slower, but memory efficient
      tmp = tcrossprod(residual_boots[,j])
      sapply(1:(nK+1), function(i) sum(tmp*all_K[[i]]))
    })
    
    
    # construct XTX, XTY matrix
    tmp = as.matrix(expand.grid(1:(nK+1), 1:(nK+1)))
    XTX_approx = matrix(sapply(1:nrow(tmp), function(l)sum(all_K[[tmp[l,1]]]*all_K[[tmp[l,2]]])), nK+1, nK+1)
    
    # XTY_boots = t(sapply(1:(nK+1), function(i)vec_K[[i]]%*%vec_y_boots))
    
  }
  
  
  
  vc_boots = sapply(1:n_boots, function(i){
    res_qp_boots = solve.QP(Dmat = XTX_approx, dvec=XTY_boots[,i], Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
    res_qp_boots$solution
  })
  
  
  heri_boots = t(t(vc_boots)/colSums(vc_boots))
  
  heri_full = vc_rehe_full/sum(vc_rehe_full)
  
  # the following based on Wald type SD
  CIvc_wald = t(sapply(1:(nK+1), function(i) vc_rehe_full[i] + c(-1, 1)*qnorm(0.975)*sd(vc_boots[i,])))
  rownames(CIvc_wald)<-c('residual', paste('K', 1:nK))
  CIvc_quan = t(sapply(1:(nK+1), function(i) vc_rehe_full[i] + quantile(vc_boots[i,]-vc_rehe_full[i],c(0.025, 0.975))))
  rownames(CIvc_quan)<-c('residual', paste('K', 1:nK))
  
  CIheri_wald = t(sapply(1:(nK+1), function(i) heri_full[i] + c(-1, 1)*qnorm(0.975)*sd(heri_boots[i,])))
  rownames(CIheri_wald)<-c('residual', paste('K', 1:nK))
  
  CIheri_quan = t(sapply(1:(nK+1), function(i) heri_full[i] + quantile(heri_boots[i,]-heri_full[i],c(0.025, 0.975))))
  rownames(CIheri_quan)<-c('residual', paste('K', 1:nK))
  
  
  
  return(list(CIvc = list(wald = CIvc_wald, quantile = CIvc_quan),
              CIheri = list(wald = CIheri_wald, quantile = CIheri_quan)))  
  
  
}
