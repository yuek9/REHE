############ 
## Kun Yue
## yuek@uw.edu
## 12/14/2020 


## use estimated network as the true network, and generate data based on the linear mixed model:
## X_i = Lambda * mu + Lambda * gamma_i + epsilon_i
## the influence matrix is Lambda, and the network information matrix has 1-weighted_Adj = {Lambda %*% t(Lambda)}^(-1)

## use the following command in the terminal
## Rscript semi_syhthetic_simulation.R <s2e> <s2g> <muvals> <iter> <target>

## s2e = noise variance
## s2g = random effect variance
## muvals = mean signal difference
## iter = the i^th iteration to run (different iter for a different random seed)
## target = 'est' for estimation, and 'truepower' for computing true power only

# s2e = 0.5
# s2g = 0.1
# muvals = 0.2
# iter = 1

library(Matrix)


inputs = commandArgs(T)
s2e = as.numeric(inputs[1])
s2g = as.numeric(inputs[2])
muvals = as.numeric(inputs[3])
iter = as.numeric(inputs[4])
target = as.character(inputs[5])

print(c(s2e, s2g, muvals, iter, target))



set.seed(1)
p=2598
n=403+117

today <- '20200827'

filepath = '...'
setwd(filepath)

out_path = paste0('Output/semi_synthetic/e',s2e, 'g', s2g, 'mu',muvals)
setwd('Output')
dir.create('semi_synthetic', showWarnings = F)
setwd('semi_synthetic')
dir.create(paste0('e',s2e, 'g', s2g, 'mu',muvals), showWarnings = F)
setwd(filepath)

load("../data/breastcancer2012_ready.rda") # use this to obtain index of differential genes
load('../data/breastcancer2012.rda')
## load the estimated network influence matrices
load(paste0(today, 'permute_breastcancer2012_network.rda'))

prec_list1 = network_info$Adj[[1]]
for(i in 1:length(prec_list1)){prec_list1[[i]] = diag(1,dim(prec_list1[[i]])[1]) - prec_list1[[i]]}
prec_list1[[1]][1:3, 1:3]
prec_list2 = network_info$Adj[[2]]
for(i in 1:length(prec_list2)){prec_list2[[i]] = diag(1,dim(prec_list2[[i]])[1]) - prec_list2[[i]]}

for(i in 1:length(prec_list1)){
  if(dim(prec_list1[[i]])[1]==1 & prec_list1[[i]][1]==0){
    prec_list1[[i]][1]<-1
    print(i)
  }
}

for(i in 1:length(prec_list2)){
  if(dim(prec_list2[[i]])[1]==1 & prec_list2[[i]][1]==0){
    print(i)
    prec_list2[[i]][1]<-1
  }
}

  
Sigma_list1 = lapply(prec_list1, function(x){solve(x)})
Sigma_list2 = lapply(prec_list2, function(x) solve(x))
D_list1 = lapply(Sigma_list1, function(x)t(chol(x)))
D_list2 = lapply(Sigma_list2, function(x)t(chol(x)))
  
DInv_list1 = lapply(D_list1, function(x)(solve(x)))
DInv_list2 = lapply(D_list2, function(x)(solve(x)))

# range(diag(bdiag(Sigma_list1)))
# range(diag(bdiag(Sigma_list2)))
  
n1=403
n2=117
n = n1+n2
  

set.seed(iter) # each data iteration reproducible
  
x1 = as.matrix(bdiag(D_list1)) %*% matrix(rnorm(n1*p), p, n1)*sqrt(s2g)+ 
    matrix(rnorm(n1*p), p, n1)*sqrt(s2e) # error term
  
x2 = as.matrix(bdiag(D_list2)) %*% matrix(rnorm(n2*p), p, n2)*sqrt(s2g)+ 
    matrix(rnorm(n2*p), p, n2)*sqrt(s2e)
  
## make sure the order of the genes are the same in both groups
if(!all(do.call(c, sapply(prec_list1, rownames)) == do.call(c, sapply(prec_list2, rownames)))){
  stop('order of genes in two groups are not the same')
}
## use this order to generate data, and set the pathways_mat
rownames(x1) <- rownames(x2) <- do.call(c, sapply(prec_list1, rownames))
  
pathways_mat = pathways_mat[, rownames(x2)]

## condition 1 has beta=0, condition 2 beta is added by the specified value 0.1~0.5 to certain gene entries
mu2 <- muvals * (match(rownames(x1), genes2affect_btw, nomatch = 0)>0) # use Jing's file, only change genes2affect_btw
x2 <- x2 + as.matrix(bdiag(D_list2)) %*% matrix(rep(mu2, n2), p, n2)

x = cbind(x1, x2)


group = c(rep(1, n1), rep(2, n2))


if(F){  # if mannually compute it
  
  x1_demean =  x1 - replicate(n1, rowMeans(x1))
  x2_demean =  x2 - replicate(n2, rowMeans(x2))
  
  XXXY = function(x, Sigma){
      Xi1 = as.vector(Sigma)
      Xi2 = as.vector(diag(1, p))
      n = ncol(x)
      XTX = n*matrix(c(sum(Xi1^2), sum(Xi1*Xi2), sum(Xi1*Xi2), sum(Xi2^2)), 2, 2)
      XTY = sapply(1:n, function(i){
        Yi = as.vector(as.vector(x[,i])%*% t(as.vector(x[,i])))
        c(sum(Xi1*Yi), sum(Xi2*Yi))
      })
      
      XTY = rowSums(XTY)
      return(list(XTX, XTY))
  }
  
  t1 = XXXY(x1_demean, as.matrix(bdiag(Sigma_list1)))
  t2 = XXXY(x2_demean, as.matrix(bdiag(Sigma_list2)))
    
    
  XTX = t1[[1]]+t2[[1]]
  XTY = t1[[2]]+t2[[2]]
  
  
  raw_res =  solve(XTX,XTY) ### problem here with NetGSA function: results not the same
  
  raw_res

}
  
library(netgsa)

time <- matrix(0, nrow=1, ncol=3, dimnames = list(NULL, c('REHE', 'reREHE', 'REML')))

# setwd('netgsa-master/netgsa-master/R_box')
# sapply(list.files(), source)
# setwd(filepath)
# library(corpcor)
# library(quadprog)
# library(genefilter)

if(target == 'est'){
  cat('run REHE \n')
  set.seed(100)
  time[1,'REHE'] = sum(system.time({
    out_REHE <- NetGSA(network_info[["Adj"]], 
                       x = x, 
                       group ,
                       pathways_mat,
                       lklMethod = "REHE", sampling = F)})[1:2])
  
  
  
  cat('run reREHE \n')
  set.seed(100)
  time[1,'reREHE'] = sum(system.time({
    out_reREHE <- NetGSA(network_info[["Adj"]], 
                         x, 
                         group, pathways_mat,
                         lklMethod = "REHE", sampling = T, 
                         sample_n = 0.1, sample_p=0.5)})[1:2])
  
  
  cat('run REML \n')
  out_REML=NULL
  set.seed(100)
  time[1,'REML'] = sum(system.time({
    out_REML <- NetGSA(network_info[["Adj"]],
                       x, group, pathways_mat,
                       lklMethod = "REML")})[1:2])
  
  out = list(REHE = out_REHE, reREHE = out_reREHE, REML = out_REML)
  
  
  
  save(list= c('time', 'out'), file = paste0(out_path, '/', today, '_iter_', iter, '.rda'))
  
}

if(target == 'truepower'){
  ## also need to compute the true power
  source('netgsa_true_power.R')
  truepower = netgsa_true_power(list(D_list1, D_list2),
                                x,
                                y=group,
                                B= pathways_mat,
                                trueVar = c(sqrt(s2g), sqrt(s2e)),
                                trueBeta = list(rep(0, p), mu2))
  
  save('truepower', file = paste0(out_path, '/truepower.rda'))
}


