##### this is the final simulation code for heritability study
##  Kun Yue, yuek@uw.edu, 12/14/2020

filepath = '..' #Pearson


library(GENESIS)
library(gdsfmt)
library(GWASTools) 
library(quadprog)
library(MASS)
library(foreach)
library(parallel)
library(doParallel)
sessionInfo()

setwd(filepath)
source('library.R')
source('reREHE.R')

input = read.table(paste0('input_simulation.txt'), header=T)
input_arg = commandArgs(T)
input_row = as.numeric(input_arg[1])
argv = as.matrix(input[as.integer(input_row),])
#argv[1] <-9000
n = as.integer(argv[1])
sig.e2 = as.numeric(argv[2])
sig.k2 = as.numeric(argv[3])
condition = as.integer(argv[4]) # condition for different kinship choice
# ratio = argv[5]
print(argv)

sig.h2=0  # for this simulation study do not consider household effect

set.seed(200)
nreps=200

In = diag(1, n)


## compute this file to save time for later replications
load(paste0('sqrt.Sigma_e_k_h', paste0(c(sig.e2, sig.k2, sig.h2), collapse ='_'), 'n', n,'c', condition, '.RData'))

if(F){
  #Kinship matrix
  Kin = diag(1, n)
  Kh1 = matrix(c(1, 0.8, 0.2, 0.8, 1, 0.4, 0.2,0.4, 1), 3, 3 )
  Kh2 = matrix(c(1, 0.05, 0.05, 0.05, 1, 0.1, 0.05, 0.1, 1), 3, 3)
  Kh3 = matrix(c(1, 0.1, 0.05, 0.1, 1, 0.3, 0.05, 0.3, 1), 3, 3)
  
  Hh = matrix(1, 3, 3)
  
  Kh = list(Kh1 ,Kh2)
  
  if(condition%in% 1:2){
    Kh = Kh[[condition]]
    for (i in 1:(n/3)){
      Kin[(i*3-2):(i*3), (i*3-2):(i*3)] = Kh
    }
    
    #Household correlation matrix
    H = diag(1, n)
    Hh= matrix(1, 3, 3)
    for (i in 1:(n/3)){
      H[(i*3-2):(i*3), (i*3-2):(i*3)] = Hh
    }
  }
  if(condition==3){
    Kin = matrix(runif(n*n, min=-0.001, max=0.001), n, n)
    Kin = (Kin+t(Kin))/2
    Kh = Kh3
    for (i in 1:(n/3)){
      Kin[(i*3-2):(i*3), (i*3-2):(i*3)] = Kh
    }
    
    #Household correlation matrix
    H = diag(1, n)
    Hh= matrix(1, 3, 3)
    for (i in 1:(n/3)){
      H[(i*3-2):(i*3), (i*3-2):(i*3)] = Hh
    }
  }
  if(condition==4){
    Kin <- getobj(paste0(filepath,'pcreap_phi_matrix_without_asian_outliers.RData'))
    dim(Kin)
    Kin.ids <- colnames(Kin)
    
    H <- matrix(0, n, n)
    # H <- getobj(paste0(filepath,"covar_matrix_HH_ID.RData"))
    # dim(H) # [1] 12803 12803
    # H = H[Kin.ids, Kin.ids]
    
    ids = 1:dim(Kin)[1]
    
    Kin = Kin[ids[1:n], ids[1:n]]
    H = H[ids[1:n],ids[1:n]]
  }
  
  
  #residual coavriance matrix for simulation
  Sigma = sig.e2*In + sig.k2*Kin + sig.h2*H
  colnames(Kin)<-1:n
  rownames(Kin)<-1:n
  colnames(H)<-1:n
  rownames(H)<-1:n
  
  svd.Sigma <- svd(Sigma)
  sqrt.Sigma <- svd.Sigma$v %*% diag(sqrt(svd.Sigma$d)) %*% t(svd.Sigma$v)
  
  save(file=paste0('sqrt.Sigma_e_k_h', paste0(c(sig.e2, sig.k2, sig.h2), collapse ='_'), 'n', n,'c', condition, '.RData'),
       list=c('sqrt.Sigma', 'Kin', 'H'))
  
}

cat('SVD preperation for simulation finished \n')

# the following generate fixed effects, but in this simulation we do not include fixed effects
w = matrix(rnorm(n)*10, ncol=1)
#W = cbind(rep(1, n), w)
#true.beta = matrix(c(2, 3), nrow=2, ncol=1)
true.beta=c(0, 0)
W = matrix(0, n, 2)

rm('Sigma')
rm('svd.Sigma')


#######################
# in the main paper do not include sparsification results


simulation_iteration=function(iter){
  set.seed(iter)
  print(iter)
  
vc_rehe_full<-vc_rehe_full_sparse<-vc_he_full<-vc_he_full_sparse<-vc_rehe_sub<-vc_rehe_sub_mean<-heri_rehe_full<-heri_rehe_full_sparse<-heri_he_full<-heri_he_full_sparse<-heri_rehe_sub<-
heri_rehe_sub_mean<-vc_CI<-vc_CI_q<-heri_CI<-heri_CI_q<-vc_reml<-vc_reml_CI<- coverage_reml<- heri_reml<- heri_reml_CI<-vc_reml_sparse<-
vc_reml_CI_sparse<- coverage_reml_sparse<- heri_reml_sparse<- heri_reml_CI_sparse<-time_CI_sparse<-time_full<- time_full_sparse_fit<- time_reml<- 
time_reml_sparse<- time_reml_sparse_fit<-time_sparsify<-time_sub<- NA


  # separately look at sparcification time, notice the starting iterations have long warm up time
  time_sparsify = c(0,0,0)
  Kin_sparse = Kin
  time_sparsify = system.time({Kin_sparse = makeSparseMatrix(Kin, thresh = 2^(-6.5)*2, verbose=F)})
  
  #generate residuals
  ind.err <- rnorm(n)
  residual <- sqrt.Sigma %*% ind.err
  
  #generate corresponding observation values
  y = W%*%true.beta+residual
  residual = y ## notice if we fit fixed coefficients, then the residual distribution Sigma is NOT directly K+I, but need projection P-x (K+I) P-x onto orthogonal column space of X
  
  if(sig.h2!=0){
    K = list(Kin_sparse, H)
    true.sig2 = c(sig.e2, sig.k2, sig.h2)
  }else{
    K = list(Kin_sparse)
    true.sig2 = c(sig.e2, sig.k2)
  }
  
  nK = length(K)
  n=length(residual)

  
  
  ##################
  ### Full REHE
  ###################
  
  if(T){
    
    cat('Full REHE.. \n')
    
    # directly perform full REHE with dense covariance matrix
    time_full = system.time({
      vc_full = fitREHE(n, residual, Kin, In)
      vc_rehe_full = vc_full$rehe
      heri_rehe_full = vc_rehe_full[2] / sum(vc_rehe_full)
      vc_he_full = vc_full$he
      heri_he_full = vc_he_full[2] / sum(vc_he_full)
      
    })
    
    cat('resample REHE .. \n')
    if(F){
    time_sub = system.time({
      vc_sub = fitreREHE(n, residual, K, In, sample_ratio = ratio, reps=50, option='median')
      vc_rehe_sub = vc_sub$rerehe
      heri_rehe_sub = vc_rehe_sub[2] / sum(vc_rehe_sub)
    })
    
    time_sub = system.time({
      vc_sub_mean = fitreREHE(n, residual, K, In, sample_ratio = ratio, reps=50, option='mean')
      vc_rehe_sub_mean = vc_sub_mean$rerehe
      heri_rehe_sub_mean = vc_rehe_sub_mean[2] / sum(vc_rehe_sub_mean)
    })
    }

    
    cat('Sparse Full REHE... \n')
    
    # perform sparsify and then rehe, little time needed for fitting REHE alone
    time_full_sparse_fit = rep(0,3)
    vc_rehe_full_sparse <-heri_rehe_full_sparse <-vc_he_full_sparse<-heri_he_full_sparse<-NA
    
    
    # time_full_sparse_fit = system.time({
    #   n_ordering = as.numeric(Kin_sparse@Dimnames[[1]])
    #   vc_full = fitREHE(n, residual[n_ordering], K=list(Kin_sparse), In)
    #   vc_rehe_full_sparse = vc_full$rehe
    #   heri_rehe_full_sparse = vc_rehe_full_sparse[2] / sum(vc_rehe_full_sparse)
    #   
    #   vc_he_full_sparse = vc_full$he
    #   heri_he_full_sparse = vc_he_full_sparse[2] / sum(vc_he_full_sparse)
    #   
    # })
    
  
    # sparified parametric bootstrapping for confidence interval construction
          
    call('\n Get Full REHE CI with sparse matrix... \n')

    if(T){

      time_CI_sparse = system.time({
        
        
        for(s in 1:length(vc_rehe_full)){
          vc_rehe_full[s] = ifelse(vc_rehe_full[s]<0, 0, vc_rehe_full[s])}
        
        varComp =  c(vc_rehe_full[2:(nK+1)], vc_rehe_full[1]) # c(random effect vc, residual vc)
        
        {
          
          Vre <- Reduce("+", mapply("*", K, varComp[1:nK], SIMPLIFY=FALSE))
          diag(Vre) <- diag(Vre) + varComp[nK+1]  
          # cholesky decomposition
          cholSigma <- chol(Vre)
        }

        
        res = REHE_CI(n, K, vc_rehe_full, In, cholSigma = cholSigma)
      })


      vc_CI = res$CIvc$wald
      vc_CI_q = res$CIvc$quantile
      heri_CI = res$CIheri$wald
      heri_CI_q = res$CIheri$quantile
      
    }
    
    
  }
  
  
  
  
  
  
  
  ########################################################################
  ##########
  ##REML
  #########
  #refer to the paper, package GENESIS is used, function fitNullMM is used to estimate VC based on AIREML (average information REML)
  if(F){
    cat('REML...')
    # do not used covariate in this simulation
    scanData = data.frame(list(scanID=1:n, outcome = y))
    time_reml = system.time({
      res = fitNullModel(scanData, outcome = 'outcome',  
                         cov.mat = K, verbose = F)
    })
    
    
    vc_reml = c(res$varComp[nK+1], res$varComp[-(nK+1)])
    
    vc_reml_CI =  varCompCI(res, prop=F)[c(nK+1, 1:nK),2:3]
    
    coverage_reml = sapply(1:(nK+1), function(i) true.sig2[i] >=vc_reml_CI[i,1] & true.sig2[i] <=vc_reml_CI[i,2])
    
    heri_reml = vc_reml[2]/sum(vc_reml)
    
    heri_reml_CI =  varCompCI(res, prop=T)[1,2:3]
    
    
  }
  
  time_reml_sparse_fit<-time_reml_sparse<-vc_reml_sparse <-vc_reml_CI_sparse<- coverage_reml_sparse<-heri_reml_sparse<-heri_reml_CI_sparse <-NA
  
  if(T){
    cat('REML Sparse...')
 
    n_ordering = as.numeric(Kin_sparse@Dimnames[[1]])
    
    scanData = data.frame(list(scanID=n_ordering, outcome = y[n_ordering]))
    rownames(scanData)<-n_ordering
    
    time_reml_sparse_fit = system.time({res = fitNullModel(scanData, outcome = 'outcome',
                                                    cov.mat = list(Kin_sparse) , verbose = F)})
    
    time_reml_sparse  = time_sparsify+time_reml_sparse_fit
    
    
    vc_reml_sparse = c(res$varComp[nK+1], res$varComp[-(nK+1)])
    
    vc_reml_CI_sparse =  varCompCI(res, prop=F)[c(nK+1, 1:nK),2:3]
    
    coverage_reml_sparse = sapply(1:(nK+1), function(i) true.sig2[i] >=vc_reml_CI_sparse[i,1] & true.sig2[i] <=vc_reml_CI_sparse[i,2])
    
    heri_reml_sparse = vc_reml_sparse[2]/sum(vc_reml_sparse)
    
    heri_reml_CI_sparse =  varCompCI(res, prop=T)[1,2:3]
    
    
  }
  
  
  


  
  return(list(
    rehe = list(
      
      vc_rehe_full = vc_rehe_full,
      vc_rehe_full_sparse = vc_rehe_full_sparse,
      vc_he_full = vc_he_full,
      vc_he_full_sparse = vc_he_full_sparse,
      vc_rehe_sub = vc_rehe_sub,
      vc_rehe_sub_mean = vc_rehe_sub_mean,
      
      
      heri_rehe_full = heri_rehe_full,
      heri_rehe_full_sparse = heri_rehe_full_sparse,
      heri_he_full = heri_he_full,
      heri_he_full_sparse = heri_he_full_sparse,
      heri_rehe_sub = heri_rehe_sub,
      heri_rehe_sub_mean = heri_rehe_sub_mean,
      
      vc_CI= vc_CI,
      vc_CI_q = vc_CI_q,
      heri_CI = heri_CI,
      heri_CI_q = heri_CI_q
                ),
    reml = list(
      est = vc_reml,
      CI = vc_reml_CI, 
      coverage = coverage_reml, 
      heri = heri_reml, 
      heri_CI = heri_reml_CI),
  
    reml_sparse = list(
      est = vc_reml_sparse,
      CI = vc_reml_CI_sparse, 
      coverage = coverage_reml_sparse, 
      heri = heri_reml_sparse, 
      heri_CI = heri_reml_CI_sparse),
    
    time = list(
      rehe_CI = time_CI_sparse,
      full = time_full, 
      full_sparse = time_full_sparse_fit, 
      reml = time_reml, 
      reml_sparse = time_reml_sparse, 
      reml_sparse_fit_only = time_reml_sparse_fit,
      make_sparse = time_sparsify,
      sub = time_sub
      )
    ))
  
}

# simulation_iteration(1)

cat('\n start parallel \n')



cl <- makeCluster(5)

registerDoParallel(cl)

res <- foreach(iter = 1:nreps,
               .packages = c("MASS", 'quadprog', 'parallel', 'GENESIS'),
               .verbose = T,
               .errorhandling = 'pass',
               .export = c('true.beta', 'w', 'W', 'n', 'Kin')) %dopar% simulation_iteration(iter)

stopCluster(cl)

cat('\n finish \n')
res[[1]]


rm('Kin')
rm('H')
rm('sqrt.Sigma')
rm('In')



save.image(paste0("sparse_parametric_bootstrap_",paste(argv, collapse='_'), "_condition_", condition, ".RData"))






