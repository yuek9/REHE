## Kun Yue, yuek@uw.edu
## 02/24/2020


##### a collection of functions to be used in simulation
library(quadprog)


# K can take sparse matrix and possibility improve speed
# need to provide the reordered residuals with ordering matching the sparsified matrix (sparsified matrix will have rows/columns reordered)
fitREHE <- function(n, residual, K, In){
  nK = length(K)
  estimates=NULL
  estimates_he = NULL
  print(class(K[[1]]))
  
  if(class(K[[1]])!='matrix'){
    
  if(class(K[[1]])=='dsyMatrix'){
    return('error: the threshold too high and dense matrix resulted')
  }
    
  if(class(K[[1]])=='dsCMatrix'){
    Kin_sparse = K[[1]]
    colindex = do.call(c,mapply(1:n, each=diff(Kin_sparse@p), FUN=rep)) # this will get nonzero entry column index
    rowindex = Kin_sparse@i+1 # this will give the nonzero row index
  }
  
  if(class(K[[1]])=='dsTMatrix'){
    Kin_sparse = K[[1]]
    colindex = Kin_sparse@j+1 # this will get nonzero entry column index
    rowindex = Kin_sparse@i+1 # this will give the nonzero row index
  }

    vec_K = Kin_sparse[cbind(rowindex, colindex)] # this is effective nonzero entries
    vec_I = In[cbind(rowindex, colindex)]
    
    vec_y = residual[rowindex]*residual[colindex]
    
    x11 = sum(vec_I*vec_I)
    x12 = sum(vec_I*vec_K)
    x22 = sum(vec_K*vec_K)
    XTX = matrix(c(x11, x12, x12, x22), nrow=2)

    XTY  = matrix(c(sum(vec_y*vec_I), sum(vec_y*vec_K)), ncol=1)
    

  }else{
    
    index=1:n
    y_matrix = as.vector(tcrossprod(residual))
    X_matrix = cbind(as.vector(In), sapply(1:nK, function(w)as.vector(K[[w]])))
    
    XTX = crossprod(X_matrix)
    XTY = crossprod(X_matrix, y_matrix)
  
  }
  
  res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
  
  estimates = res_qp$solution
  estimates_he = res_qp$unconstrained.solution
  
  # non-negative correction for HE estimations
  for(s in 1:length(estimates_he)){
    estimates_he[s] = ifelse(estimates_he[s]<0, 0, estimates_he[s])}
  
  # numerical correction for REHE estimations (due to machine error)
  for(s in 1:length(estimates)){
    estimates[s] = ifelse(estimates[s]<0, 0, estimates[s])}
  
  return(list(rehe = estimates, he = estimates_he))
}

REHE_CI = function(n, K,
                   vc_rehe_full, In, cholSigma=NULL){
  
  # vc_rehe_full should be in order of c(residual, K1, K2...)
  # output confidence interval for both sd based CI (wald type) and empirical CI (percentile of difference)
  
  nK = length(K)
  
  # get the sd estimation from 100 replications; could also use quantile based
  n_boots = 100
  
  
  if(is.null(cholSigma)){ # will report issue if covariance matrix cannot have chol; for now do not really this this step
    
    tmp = list(Matrix::chol(K[[1]])) # this is fast, gives R, where R'R = Sigma
    
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
    
    all_K = list(In) # always put In as the first one, so can speed up a little bit
    all_K[2:(nK+1)] = K
    
    vec_K = list( as.vector(In))
    vec_K[2:4]= lapply(K, as.vector)
    #dim(residual_boots)
    
    # vec_y_boots = residual_boots[row(In),]*residual_boots[col(In),] # this is too large to store; do this one by one
    # list_y_boots = lapply(1:n_boots, function(j)tcrossprod(residual_boots[,j])) #same size
    
    #dim(vec_y_boots)
    XTY_boots = sapply(1:n_boots, function(j){  # this seems slower, but memory efficient
      tmp = tcrossprod(residual_boots[,j])
      c(sum(diag(tmp)), sapply(2:(nK+1), function(i) sum(tmp*all_K[[i]])))
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


check_coverage = function(CI, truth){
  if(length(truth)>1){
    sapply(1:length(truth), function(i)CI[i,1]<=truth[i] & CI[i,2]>=truth[i])
  }else{
    CI[1]<=truth & CI[2]>=truth
  }
}
