## this code compute the true power with given network influence matrix, muvals, and variance components
## ref:    1) Shojaie A and Michailidis G (2009), JCB 
##	 	     2) Shojaie A and Michailidis G (2010), SAGMB



netgsa_true_power<-function(
  InfMat,			    #influence matrix 
  x, 			      	#the p x n data matrix
  y, 			      	#vector of class indicators of length n (for treatment or control indicator)
  B, 		    	  	#indicator matrix for pathways (npath x p)
  varEstMethod = "Newton", #method for estimating variance components
  verbose = FALSE,     #Whether to display the estimated values for s2g,s2e and dof, default to be False.
  trueVar = NULL,      #Input the true variance components, as c(sqrt(s2g), sqrt(s2e))
  trueBeta = NULL     # true beta vectors, a list as list(beta for x1, beta for x2)
){
  
  ##-----------------
  p = dim(x)[1] #No. of genes
  n = length(y) #No. of samples in total
  npath = dim(B)[1] #No. of gene sets to consider 
  
  y = as.integer(as.factor(y))
  n1 = sum(y == 1)
  n2 = sum(y == 2)
  
  X1 = x[, (y == 1)]
  X2 = x[, (y == 2)]
  
  
  
  ##-----------------
  ##setting up control parameters for the var estimation procedures
  varEstCntrl = list(lklMethod = "REML", #method for likelihood, REML, ML
                     s2profile = "se",   #parameter to profile out (se/sg)
                     lb = 0.5,          #lower bound for the line search method
                     ub = 100,           #upper bound for the line search method
                     tol = 0.01,         # tolerance level
                     useGrad = FALSE     #whether to use gradient (for DER method)
  )     

  
  ##--------------------
  #pre-examination of adjancency matrix and influence matrix property
  if (dim(x)[2] != n) {
    stop("The dimensions of the data matrix and class vector don't match.")
  }
  
  if (dim(B)[2] != p) {
    stop("The dimensions of the data matrix and indicator matrix don't match.")
  }
  
  if (length(unique(y)) != 2) {
    stop("There should be 2 unique classes in the class indicator!")
  }
  

  

  
  ##Matrices Needed 
  Ip = diag(rep(1, p))
  
  D1_list = InfMat[[1]]            ####### InfMat is the infunence matrix Lambda
  D2_list = InfMat[[2]]
  
  D1Inv_list = lapply(D1_list, solve)
  D2Inv_list = lapply(D2_list, solve)
  D1D1_list = lapply(D1_list, function(x) x %*% t(x))
  D2D2_list = lapply(D2_list, function(x) x %*% t(x))
  
  D1D1Inv_list = lapply(D1_list, function(x)solve(t(x) %*% x))
  D2D2Inv_list = lapply(D2_list, function(x)solve(t(x) %*% x))

  D1 = as.matrix(bdiag(D1_list))
  D2 = as.matrix(bdiag(D2_list))
  
  D1Inv = as.matrix(bdiag(D1Inv_list))
  D2Inv = as.matrix(bdiag(D2Inv_list))
  
  D1D1 = as.matrix(bdiag(D1D1_list))
  D2D2 = as.matrix(bdiag(D2D2_list))
  
  
  D1D1Inv = as.matrix(bdiag(D1D1Inv_list))
  D2D2Inv = as.matrix(bdiag(D2D2Inv_list))

  
  ##--------------------------
  ##--------------------------
  ##calculate test related statistics, per simulation
  
  ##Building the "contrast" matrix L, see Result in the paper
  L1 = (B %*% D1) * B
  L2 = (B %*% D2) * B
  LN = cbind(-L1, L2)
  

  
  
  ##--------------------------------------------------------------------
  ## calculate the true power when the true population mean and variance components are known.
  ## ncp is the non-centerality parameter
  ncp = matrix(0, npath, 1)
  df.alt = matrix(0, npath, 1)
  num.alt = matrix(0, npath, 1)
  if ((1-is.null(trueVar))*(1-is.null(trueBeta)) == 1) {
    ## Get the population parameters
    beta1.true = trueBeta[[1]]
    beta2.true = trueBeta[[2]]
    s2g.true = trueVar[1]^2
    s2e.true = trueVar[2]^2
    mctildi.true = s2e.true * D1D1Inv + s2g.true * Ip
    mttildi.true = s2e.true * D2D2Inv + s2g.true * Ip
    
    ## matrices needed to calculate the true degree of freedom.
    mc = chol2inv(chol(s2g.true * D1D1 + s2e.true * Ip))
    mca = mc %*% D1D1
    mt = chol2inv(chol(s2g.true * D2D2 + s2e.true * Ip))
    mta = mt %*% D2D2
    
    #trace of matrices needed for H (and therefore KK)
    t11c = sum(diag(mca %*% mca))
    t11t = sum(diag(mta %*% mta))
    t22c = sum(diag(mc %*% mc))
    t22t = sum(diag(mt %*% mt))
    t12c = sum(diag(mca %*% mc))
    t12t = sum(diag(mta %*% mt))
    
    EH11 = (1/2) * (t11c + t11t)
    EH22 = (1/2) * (t22c + t22t)
    EH12 = (1/2) * (t12c + t12t) #(-1/2) * (tt1 + tt2) + (1/2) * (tt3 + tt4)
    
    Kmat = matrix(c(EH11, EH12, EH12, EH22), 2, 2, byrow = TRUE)
    KmatInv = solve(Kmat)
    
    # Calculate the non-centrality parameter and the df under the alternative hypothesis
    for (rr in 1:npath){
      Lrow = t(as.matrix(LN[rr, ])) #single row of L3
      Lrow1 = t(as.matrix(Lrow[, 1:p]))
      Lrow2 = t(as.matrix(Lrow[, (p + 1):(2 * p)]))
      LC11Lprime.true = (1/n2) * Lrow2 %*% mttildi.true %*% t(Lrow2) + (1/n1) * Lrow1 %*% mctildi.true %*% t(Lrow1)
      g1 = (1/n2) * (Lrow2 %*% t(Lrow2)) + (1/n1) * (Lrow1 %*% t(Lrow1))
      g2 = (1/n2) * Lrow2 %*% D2D2Inv %*% t(Lrow2) + (1/n1) * Lrow1 %*% D1D1Inv %*% t(Lrow1)
      g = matrix(c(g1, g2), 2, 1)
      num.alt[rr] = Lrow2 %*% beta2.true + Lrow1 %*% beta1.true
      ncp[rr] = (Lrow2 %*% beta2.true + Lrow1 %*% beta1.true)/sqrt(LC11Lprime.true)	
      df.alt[rr] = 2 * (LC11Lprime.true)^2/(t(g) %*% KmatInv %*% g)
      
      if (df.alt[rr]<2) {df.alt[rr]=2}
    }
    
  }
  
  # Find the rejection region (cutoff in this case) from the null distribution; df should be based on the true vc and bias
  cutoff = qt(p = 1 - 0.025, df = df.alt, ncp = 0) 
  
  ##first row for the true power of each subnet test
  sigInd1 = pt(q = -cutoff, df = df.alt, ncp = ncp, lower.tail = TRUE) + 
    pt(q = cutoff, df = df.alt, ncp = ncp, lower.tail = FALSE)
  
  output = sigInd1
  
  
  
  
  return(output)
  
}
