######
# resample REHE function
#####

fitreREHE <- function(n, residual, K, In, sample_ratio = 0.1, reps=50, option='median'){
  nK = length(K)
  estimates=NULL
  estimates_he = NULL
  print(class(K[[1]]))
  
  if(class(K[[1]])!='matrix'){
    stop('Correlation matrix type should be matrix. Do not use reREHE if the correlation matrix is already sparsified.')
  }
  
  
  # if(class(K[[1]])!='matrix'){
  #   
  #   
  #   if(class(K[[1]])=='dsyMatrix'){
  #     return('error: the threshold too high and dense matrix resulted')
  #   }
  #   
  #   if(class(K[[1]])=='dsCMatrix'){
  #     Kin_sparse = K[[1]]
  #     colindex = do.call(c,mapply(1:n, each=diff(Kin_sparse@p), FUN=rep)) # this will get nonzero entry column index
  #     rowindex = Kin_sparse@i+1 # this will give the nonzero row index
  #   }
  #   
  #   if(class(K[[1]])=='dsTMatrix'){
  #     Kin_sparse = K[[1]]
  #     colindex = Kin_sparse@j+1 # this will get nonzero entry column index
  #     rowindex = Kin_sparse@i+1 # this will give the nonzero row index
  #   }
  #   
  #   
  #   
  #   vec_K = Kin_sparse[cbind(rowindex, colindex)] # this is effective nonzero entries
  #   vec_I = In[cbind(rowindex, colindex)]
  #   
  #   vec_y = residual[rowindex]*residual[colindex]
  #   
  #   x11 = sum(vec_I*vec_I)
  #   x12 = sum(vec_I*vec_K)
  #   x22 = sum(vec_K*vec_K)
  #   XTX = matrix(c(x11, x12, x12, x22), nrow=2)
  #   
  #   XTY  = matrix(c(sum(vec_y*vec_I), sum(vec_y*vec_K)), ncol=1)
  #   
  #   
  # }else{
  #   
  #   index=1:n
  #   
  #   y_matrix = as.vector(tcrossprod(residual))
  #   X_matrix = cbind(as.vector(In), sapply(1:nK, function(w)as.vector(K[[w]])))
  #   
  #   XTX = crossprod(X_matrix)
  #   XTY = crossprod(X_matrix, y_matrix)
  #   
  # }
  
  
  yyT = tcrossprod(residual)
  
  reREHE_res = replicate(n=reps, {
    
  # for one subsample
  sub_index = sample(n, size=round(sample_ratio*n), replace=T) # sample the responses
  
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
  
  res_qp = solve.QP(Dmat = XTX, dvec=XTY, Amat = diag(1,nK+1), bvec = rep(0, nK+1), meq=0, factorized=FALSE)
  
  estimates = res_qp$solution
  estimates_he = res_qp$unconstrained.solution
  
  # non-negative correction for HE estimations
  for(s in 1:length(estimates_he)){
    estimates_he[s] = ifelse(estimates_he[s]<0, 0, estimates_he[s])}
  
  # numerical correction for REHE estimations (due to machine error)
  for(s in 1:length(estimates)){
    estimates[s] = ifelse(estimates[s]<0, 0, estimates[s])}
  
  estimates
  })
  
  if(option=='mean'){
    estimates = apply(reREHE_res, 1, mean)
  }
  
  if(option=='median'){
    estimates = apply(reREHE_res, 1, median)
  }
  

  return(list(rerehe = estimates))
}
