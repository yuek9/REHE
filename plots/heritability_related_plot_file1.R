## summarize wrap up simulation results (part 1: produce the R images that contain necessary data matrices for plots)
library(ggplot2)
library(ggpubr)
library(scales)
library(cowplot)
filepath = '\\\\gcc-fs4\\gcc-fs4'
setwd(paste0(filepath, '/thornton/Kun_Yue_HCHS_SOL/Code/wrap_up_files'))
theme_set(theme_bw())

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
con4 = seq(4, 80, 4)
label_condition = function(i){
  paste('setting: ', i)
}
label_setting = c(
  '1'=paste0('e2: ', 0.1, ', k2: ', 0.1),
  '2'=paste0('e2: ', 0.04, ', k2: ', 0.1),
  '3'=paste0('e2: ', 0.1, ', k2: ', 0.04),
  '4'=paste0('e2: ', 0.01, ', k2: ', 0.1),
  '5'=paste0('e2: ', 0.1, ', k2: ', 0.01)
  
)


## heritability study


#########
# for main paper results
#########

{
  
  # results for simulation with sparsification when constructing confidence interval
  savename = paste0('sparse_parametric_bootstrap_')
  
  {
  table_point_est = matrix(NA, nrow=4*5*4, ncol=9*6+4)
  colnames(table_point_est) <- c('n', 'e2','k2','condition', 
                                      'Bias rehe e', 'Var rehe e', 'MSE rehe e','Bias rehe k', 'Var rehe k', 'MSE rehe k','Bias rehe heri', 'Var rehe heri', 'MSE rehe heri', #9
                                      'Bias rehe sparse e', 'Var rehe sparse e', 'MSE rehe sparse e','Bias rehe sparse k', 'Var rehe sparse k', 'MSE rehe sparse k','Bias rehe sparse heri', 'Var rehe sparse heri', 'MSE rehe sparse heri', #9
                                      
                                      'Bias he e', 'Var he e', 'MSE he e','Bias he k', 'Var he k', 'MSE he k','Bias he heri', 'Var he heri', 'MSE he heri', #9
                                      'Bias he sparse e', 'Var he sparse e', 'MSE he sparse e','Bias he sparse k', 'Var he sparse k', 'MSE he sparse k','Bias he sparse heri', 'Var he sparse heri', 'MSE he sparse heri', #9
                                      
                                      'Bias reml e', 'Var reml e', 'MSE reml e','Bias reml k', 'Var reml k', 'MSE reml k','Bias reml heri', 'Var reml heri', 'MSE reml heri',
                                      'Bias reml e', 'Var reml sparse e', 'MSE reml sparse e','Bias reml sparse k', 'Var reml sparse k', 'MSE reml sparse k','Bias reml sparse heri', 'Var reml sparse heri', 'MSE reml sparse heri'
  )
  
  table_CI = matrix(NA, nrow=4*5*4, ncol=4+3*2+4*2+3)
  colnames(table_CI) = c('n', 'e2','k2','condition', #4
                         
                         
  
                        'e rehe est', 'k rehe est', 'heri rehe est', #3
                        'e rehe sparse est', 'k rehe sparse est', 'heri rehe sparse est', #3
                        'e rehe quantile', 'k rehe quantile', 'heri rehe quantile', #3
                        'e reml', 'k reml', 'heri reml', 'REML NA prop', #4
                        'e reml sparse', 'k reml sparse', 'heri reml sparse', 'REML NA prop sparse' #4
  )
  
  
  # this records CPU time, which adds up system.time first two entries
  # notice that for MakeSpraseMatrix, it is sometimes very strange with huge computation time (typically first several iterations), I will summarize with median time
  
  table_time = matrix(NA, nrow=4*5*4, ncol=4+3+3)
  colnames(table_time) = c('n', 'e2','k2','condition', #4
                           'bootstrap REHE CI','REHE sparse fit', 'REHE fit', #3
                            'reml sparse fit', 'reml', 'sparsify time'
  )
  
  
  he_est = list(NULL)
  
  rehe_full_est = list(NULL)
  
  reml_est = list(NULL)
  
  
  
  counter = 0
  for(ww in 1:80){
    counter=counter+1 
    
    try({
      cat(ww, '\n')
      
      filepath = '..'
      input = read.table(paste0(filepath, '/input_simulation.txt'), header=T)
      input_row=ww
      
      argv = as.matrix(input[as.integer(input_row),])
      nK=1

      n = as.integer(argv[1])
      sig.e2 = as.numeric(argv[2])
      sig.k2 = as.numeric(argv[3])
      condition = as.integer(argv[4]) 
      
      true.heri = sig.k2 / sum(sig.e2, sig.k2)
      
      table_CI[counter,1:4]<-table_point_est[counter,1:4]<-table_time[counter, 1:4]<-argv

      load(paste0(filepath, "output/" ,savename,paste(argv, collapse='_'), "_condition_", condition, ".RData"))
      
      gather = function(tmp){
        if(class(tmp)=='list'){
          do.call(cbind, tmp)
        }else{
          tmp
        }
      }


      REHE =  list(e_k = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_full)),
                   e_k_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_full_sparse)),
                   
                   heri = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_full)),
                   heri_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_full_sparse)),
                   
                   e_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[,1])), # these based on dense full point + sparse sd
                   k_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[,2])),
                   heri_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI)),
                   
                   # quantile based confidence interval
                   e_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI_q[,1])), # these based on dense full point + sparse sd
                   k_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI_q[,2])),
                   heri_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI_q)),
                   
                   # sparse REHE centered confidence interval (Wald type only)
                   e_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[,1] -  res[[i]]$rehe$vc_rehe_full[1] +  res[[i]]$rehe$vc_rehe_full_sparse[1])), # these based on sparse full point + sparse sd
                   k_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[,2] -  res[[i]]$rehe$vc_rehe_full[2] +  res[[i]]$rehe$vc_rehe_full_sparse[2])),
                   heri_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI - res[[i]]$rehe$heri_rehe_full +  res[[i]]$rehe$heri_rehe_full_sparse))
      )   
      
      HE = list(
        e_k = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_he_full)),
        e_k_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_he_full_sparse)),
        
        heri = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_he_full)),
        heri_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_he_full_sparse))

      )
      
      reml = list(e_k=gather(sapply(1:nreps, function(i) res[[i]]$reml$est)),
                  e_k_sparse=gather(sapply(1:nreps, function(i) res[[i]]$reml_sparse$est)),
                  
                  e_CI = sapply(1:nreps, function(i) res[[i]]$reml$CI[1, ]),
                  k_CI = sapply(1:nreps, function(i) res[[i]]$reml$CI[2, ]),
                  heri=gather(sapply(1:nreps, function(i) res[[i]]$reml$heri)),
                  heri_CI = sapply(1:nreps, function(i) res[[i]]$reml$heri_CI),
                  eCI_cov =gather(sapply(1:nreps, function(i) res[[i]]$reml$coverage[1])),
                  kCI_cov=gather(sapply(1:nreps, function(i) res[[i]]$reml$coverage[2])),
                  heriCI_cov = gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml$heri_CI, true.heri))),
                  
                  
                  
                  e_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$CI[1, ]),
                  k_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$CI[2, ]),
                  heri_sparse =gather(sapply(1:nreps, function(i) res[[i]]$reml_sparse$heri)),
                  heri_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$heri_CI),
                  eCI_cov_sparse =gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$CI[1, ], sig.e2 ))),
                  kCI_cov_sparse =gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$CI[2, ], sig.k2 ))),
                  heriCI_cov_sparse = gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$heri_CI, true.heri)))
                  
                  
                  
                  
      )
      
      
      he_est[[counter]] = HE
      
      rehe_full_est[[counter]] = REHE
      
      reml_est[[counter]] = reml
      
      summ = function(xx, true){
        tmp = rep(0, 3)
        if(class(xx)=='matrix'){
          for(i in 1:nrow(xx)){
            x = as.vector(xx[i,])
            tmp = rbind(tmp, c(mean(x)-true[i], var(x), (mean(x)-true[i])^2+var(x)))
          }
          return(t(tmp[-1,]))
        }else{
          return(c(mean(xx)-true[1], var(xx), (mean(xx)-true[1])^2+var(xx)))
        }
      }
      
      table_point_est[counter, c(5:58)]<-c(
        
        c(summ(REHE$e_k, c(sig.e2, sig.k2))),
        summ(REHE$heri, true.heri),
        
        
        c(summ(REHE$e_k_sparse, c(sig.e2, sig.k2))),
        summ(REHE$heri_sparse,true.heri),
        
        c(summ(HE$e_k, c(sig.e2, sig.k2))),
        summ(HE$heri, true.heri),
        
        
        c(summ(HE$e_k_sparse, c(sig.e2, sig.k2))),
        summ(HE$heri_sparse,true.heri),
        
        c(summ(reml$e_k, c(sig.e2, sig.k2))),
        summ(reml$heri,true.heri),
        
        
        c(summ(reml$e_k_sparse, c(sig.e2, sig.k2))),
        summ(reml$heri_sparse,true.heri)

      )
      
      
      # bootstrap for CI, # full sparse , #full with dense Kin
      table_time[counter, 5:7]<-c(  median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$rehe_CI[1:2])))),
                                           median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$full_sparse[1:2])))),
                                           median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$full[1:2]))))
      )
      
      
      
      
      
      # reml sparse fit only , # reml based on dense Kin, # sparsify time
      table_time[counter, 8:10]<-c(   median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$reml_sparse_fit_only[1:2])))),
                                            median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$reml[1:2])))),
                                            median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$make_sparse[1:2]))))
      )
      
      
      
      REHE_sd = list(low = rbind(REHE$e_CI[1,], REHE$k_CI[1,]), up=rbind(REHE$e_CI[2,], REHE$k_CI[2,]), heri = t(REHE$heri_CI))
      
      REHE_quantile = list(low = rbind(REHE$e_CI_q[1,], REHE$k_CI_q[1,]), up=rbind(REHE$e_CI_q[2,], REHE$k_CI_q[2,]), heri = t(REHE$heri_CI_q))
      
      REHE_sd_sparse = list(low = rbind(REHE$e_CI_sparse[1,], REHE$k_CI_sparse[1,]), up=rbind(REHE$e_CI_sparse[2,], REHE$k_CI_sparse[2,]), heri = t(REHE$heri_CI_sparse))
      

      
      true_sig = c(sig.e2, sig.k2)
      true_heri = true_sig[2]/sum(true_sig)
      
      
      
      table_CI[counter, 5:13]<- c(
        sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_sd$low[i,], REHE_sd$up[i,]), rep(true_sig[i], length(REHE_sd$up[i,])))))
        ,mean(check_coverage(REHE_sd$heri, rep(true_heri, length(REHE$heri))))
        ,
        
        
        sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_sd_sparse$low[i,], REHE_sd_sparse$up[i,]), rep(true_sig[i], length(REHE_sd_sparse$up[i,])))))
        ,mean(check_coverage(REHE_sd_sparse$heri, rep(true_heri, length(REHE$heri_sparse)))),
        
        sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_quantile$low[i,], REHE_quantile$up[i,]), rep(true_sig[i], length(REHE_quantile$up[i,])))))
        ,mean(check_coverage(REHE_quantile$heri, rep(true_heri, length(REHE$heri))))
        
      )
      
      
      table_CI[counter,14:16]<-c(mean(reml$eCI_cov, na.rm=T), mean(reml$kCI_cov, na.rm=T), mean(reml$heriCI_cov, na.rm=T))
      table_CI[counter, 17] = mean(is.na(reml$heriCI_cov))
      
      table_CI[counter, 18:20]<-c(mean(reml$eCI_cov_sparse, na.rm=T), mean(reml$kCI_cov_sparse, na.rm=T), mean(reml$heriCI_cov_sparse, na.rm=T))
      table_CI[counter, 21] = mean(is.na(reml$heriCI_cov_sparse))
      
    })
  }
  
  rm('In')
  save.image('new_wrap_up_summary.RData')
  
  }
  
  
  
  # results without sparsification for confidence interval construction (anything about sparsification will be in the supplementary)
  savename = paste0('dense_parametric_bootstrap_')
  
  {
    table_point_est = matrix(NA, nrow=4*5*4, ncol=9*6+4)
    colnames(table_point_est) <- c('n', 'e2','k2','condition', 
                                   'Bias rehe e', 'Var rehe e', 'MSE rehe e','Bias rehe k', 'Var rehe k', 'MSE rehe k','Bias rehe heri', 'Var rehe heri', 'MSE rehe heri', #9
                                   'Bias rehe sparse e', 'Var rehe sparse e', 'MSE rehe sparse e','Bias rehe sparse k', 'Var rehe sparse k', 'MSE rehe sparse k','Bias rehe sparse heri', 'Var rehe sparse heri', 'MSE rehe sparse heri', #9
                                   
                                   'Bias he e', 'Var he e', 'MSE he e','Bias he k', 'Var he k', 'MSE he k','Bias he heri', 'Var he heri', 'MSE he heri', #9
                                   'Bias he sparse e', 'Var he sparse e', 'MSE he sparse e','Bias he sparse k', 'Var he sparse k', 'MSE he sparse k','Bias he sparse heri', 'Var he sparse heri', 'MSE he sparse heri', #9
                                   
                                   'Bias reml e', 'Var reml e', 'MSE reml e','Bias reml k', 'Var reml k', 'MSE reml k','Bias reml heri', 'Var reml heri', 'MSE reml heri',
                                   'Bias reml e', 'Var reml sparse e', 'MSE reml sparse e','Bias reml sparse k', 'Var reml sparse k', 'MSE reml sparse k','Bias reml sparse heri', 'Var reml sparse heri', 'MSE reml sparse heri'
    )
    
    table_CI = matrix(NA, nrow=4*5*4, ncol=4+3*2+4*2+3)
    colnames(table_CI) = c('n', 'e2','k2','condition', #4
                           
                           
                           
                           'e rehe est', 'k rehe est', 'heri rehe est', #3
                           'e rehe sparse est', 'k rehe sparse est', 'heri rehe sparse est', #3
                           'e rehe quantile', 'k rehe quantile', 'heri rehe quantile', #3
                           'e reml', 'k reml', 'heri reml', 'REML NA prop', #4
                           'e reml sparse', 'k reml sparse', 'heri reml sparse', 'REML NA prop sparse' #4
    )
    
    
    # this records CPU time, which adds up system.time first two entries
    # notice that for MakeSpraseMatrix, it is sometimes very strange with huge computation time (typically first several iterations), I will summarize with median time
    
    table_time = matrix(NA, nrow=4*5*4, ncol=4+3+3)
    colnames(table_time) = c('n', 'e2','k2','condition', #4
                             'bootstrap REHE CI','REHE sparse fit', 'REHE fit', #3
                             'reml sparse fit', 'reml', 'sparsify time'
    )
    
    
    he_est = list(NULL)
    
    rehe_full_est = list(NULL)
    
    reml_est = list(NULL)
    
    
    
    counter = 0
    for(ww in 1:80){
      counter=counter+1 
      
      try({
        cat(ww, '\n')
        
        filepath = '..'
        input = read.table(paste0(filepath, 'input_simulation.txt'), header=T)
        input_row=ww
        
        argv = as.matrix(input[as.integer(input_row),])
        nK=1
        
        n = as.integer(argv[1])
        sig.e2 = as.numeric(argv[2])
        sig.k2 = as.numeric(argv[3])
        condition = as.integer(argv[4]) 
        
        true.heri = sig.k2 / sum(sig.e2, sig.k2)
        
        table_CI[counter,1:4]<-table_point_est[counter,1:4]<-table_time[counter, 1:4]<-argv
        
        load(paste0(filepath, "output/" ,savename,paste(argv, collapse='_'), "_condition_", condition, ".RData"))
        
        gather = function(tmp){
          if(class(tmp)=='list'){
            do.call(cbind, tmp)
          }else{
            tmp
          }
        }
        
        
        REHE =  list(e_k = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_full)),
                     e_k_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_full_sparse)),
                     
                     heri = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_full)),
                     heri_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_full_sparse)),
                     
                     
                     ## attention some difference in output format here
                     e_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[1,])), # these based on dense full point + sparse sd
                     k_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[2,])),
                     heri_CI = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI[2,])),
                     
                     # quantile based confidence interval
                     e_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI_q[1,])), # these based on dense full point + sparse sd
                     k_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI_q[2,])),
                     heri_CI_q = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI_q[2,])),
                     
                     # sparse REHE centered confidence interval (Wald type only)
                     e_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[1,] -  res[[i]]$rehe$vc_rehe_full[1] +  res[[i]]$rehe$vc_rehe_full_sparse[1])), # these based on sparse full point + sparse sd
                     k_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_CI[2,] -  res[[i]]$rehe$vc_rehe_full[2] +  res[[i]]$rehe$vc_rehe_full_sparse[2])),
                     heri_CI_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_CI - res[[i]]$rehe$heri_rehe_full +  res[[i]]$rehe$heri_rehe_full_sparse))
        )   
        
        HE = list(
          e_k = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_he_full)),
          e_k_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_he_full_sparse)),
          
          heri = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_he_full)),
          heri_sparse = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_he_full_sparse))
          
        )
        
        reml = list(e_k=gather(sapply(1:nreps, function(i) res[[i]]$reml$est)),
                    e_k_sparse=gather(sapply(1:nreps, function(i) res[[i]]$reml_sparse$est)),
                    
                    e_CI = sapply(1:nreps, function(i) res[[i]]$reml$CI[1, ]),
                    k_CI = sapply(1:nreps, function(i) res[[i]]$reml$CI[2, ]),
                    heri=gather(sapply(1:nreps, function(i) res[[i]]$reml$heri)),
                    heri_CI = sapply(1:nreps, function(i) res[[i]]$reml$heri_CI),
                    eCI_cov =gather(sapply(1:nreps, function(i) res[[i]]$reml$coverage[1])),
                    kCI_cov=gather(sapply(1:nreps, function(i) res[[i]]$reml$coverage[2])),
                    heriCI_cov = gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml$heri_CI, true.heri))),
                    
                    # does not implement sparse
                    
                    e_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$CI),
                    k_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$CI),
                    heri_sparse =gather(sapply(1:nreps, function(i) res[[i]]$reml_sparse$heri)),
                    heri_CI_sparse = sapply(1:nreps, function(i) res[[i]]$reml_sparse$heri_CI),
                    eCI_cov_sparse =gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$CI, sig.e2 ))),
                    kCI_cov_sparse =gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$CI, sig.k2 ))),
                    heriCI_cov_sparse = gather(sapply(1:nreps, function(i) check_coverage(res[[i]]$reml_sparse$heri_CI, true.heri)))
                    
                    
                    
                    
        )
        
        
        he_est[[counter]] = HE
        
        rehe_full_est[[counter]] = REHE
        
        reml_est[[counter]] = reml
        
        summ = function(xx, true){
          tmp = rep(0, 3)
          if(class(xx)=='matrix'){
            for(i in 1:nrow(xx)){
              x = as.vector(xx[i,])
              tmp = rbind(tmp, c(mean(x)-true[i], var(x), (mean(x)-true[i])^2+var(x)))
            }
            return(t(tmp[-1,]))
          }else{
            return(c(mean(xx)-true[1], var(xx), (mean(xx)-true[1])^2+var(xx)))
          }
        }
        
        table_point_est[counter, c(5:58)]<-c(
          
          c(summ(REHE$e_k, c(sig.e2, sig.k2))),
          summ(REHE$heri, true.heri),
          
          
          # c(summ(REHE$e_k_sparse, c(sig.e2, sig.k2))),
          # summ(REHE$heri_sparse,true.heri),
          rep(NA, 9),
              
          
          c(summ(HE$e_k, c(sig.e2, sig.k2))),
          summ(HE$heri, true.heri),
          
          
          # c(summ(HE$e_k_sparse, c(sig.e2, sig.k2))),
          # summ(HE$heri_sparse,true.heri),
          rep(NA, 9),
          
          c(summ(reml$e_k, c(sig.e2, sig.k2))),
          summ(reml$heri,true.heri),
          
          
          # c(summ(reml$e_k_sparse, c(sig.e2, sig.k2))),
          # summ(reml$heri_sparse,true.heri)
          rep(NA, 9)
        )
        
        
        # bootstrap for CI, # full sparse , #full with dense Kin
        table_time[counter, 5:7]<-c(  median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$rehe_CI[1:2])))),
                                      median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$full_sparse[1:2])))),
                                      median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$full[1:2]))))
        )
        
        
        
        
        
        # reml sparse fit only , # reml based on dense Kin, # sparsify time
        table_time[counter, 8:10]<-c(   median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$reml_sparse_fit_only[1:2])))),
                                        median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$reml[1:2])))),
                                        median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$make_sparse[1:2]))))
        )
        
        
        
        REHE_sd = list(low = rbind(REHE$e_CI[1,], REHE$k_CI[1,]), up=rbind(REHE$e_CI[2,], REHE$k_CI[2,]), heri = t(REHE$heri_CI))
        
        REHE_quantile = list(low = rbind(REHE$e_CI_q[1,], REHE$k_CI_q[1,]), up=rbind(REHE$e_CI_q[2,], REHE$k_CI_q[2,]), heri = t(REHE$heri_CI_q))
        
        REHE_sd_sparse = list(low = rbind(REHE$e_CI_sparse[1,], REHE$k_CI_sparse[1,]), up=rbind(REHE$e_CI_sparse[2,], REHE$k_CI_sparse[2,]), heri = t(REHE$heri_CI_sparse))
        
        
        
        true_sig = c(sig.e2, sig.k2)
        true_heri = true_sig[2]/sum(true_sig)
        
        
        
        table_CI[counter, 5:13]<- c(
          sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_sd$low[i,], REHE_sd$up[i,]), rep(true_sig[i], length(REHE_sd$up[i,])))))
          ,mean(check_coverage(REHE_sd$heri, rep(true_heri, length(REHE$heri))))
          ,
          
          
          sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_sd_sparse$low[i,], REHE_sd_sparse$up[i,]), rep(true_sig[i], length(REHE_sd_sparse$up[i,])))))
          ,mean(check_coverage(REHE_sd_sparse$heri, rep(true_heri, length(REHE$heri_sparse)))),
          
          sapply(1:(nK+1), function(i)mean(check_coverage(cbind(REHE_quantile$low[i,], REHE_quantile$up[i,]), rep(true_sig[i], length(REHE_quantile$up[i,])))))
          ,mean(check_coverage(REHE_quantile$heri, rep(true_heri, length(REHE$heri))))
          
        )
        
        
        table_CI[counter,14:16]<-c(mean(reml$eCI_cov, na.rm=T), mean(reml$kCI_cov, na.rm=T), mean(reml$heriCI_cov, na.rm=T))
        table_CI[counter, 17] = mean(is.na(reml$heriCI_cov))
        
        table_CI[counter, 18:20]<-c(mean(reml$eCI_cov_sparse, na.rm=T), mean(reml$kCI_cov_sparse, na.rm=T), mean(reml$heriCI_cov_sparse, na.rm=T))
        table_CI[counter, 21] = mean(is.na(reml$heriCI_cov_sparse))
        
      })
      print( table_CI[counter, 5:13])
    }
    
    
    dense_table_time = table_time
    dense_table_CI = table_CI
    dense_table_point_est = table_point_est
    dense_he_est = he_est
    dense_rehe_full_est = rehe_full_est
    dense_reml_est = reml_est
    
    rm('In')
    
    save.image('new_dense_wrap_up_summary.RData')
    
  }
  
  
  # point estimation results based on subsample methods: reREHE (based on mean, median, and subsampling ratio 0.05, 0.1)
  savename = paste0('dense_reREHE_only_')
  
  {
    table_point_est = (matrix(NA, nrow=4*5*4*2*2, ncol=15))
    colnames(table_point_est) <- c('n', 'e2','k2','condition',  'ratio','mean1median2',
                                   'Bias rerehe e', 'Var rerehe e', 'MSE rerehe e','Bias rerehe k', 'Var rerehe k', 'MSE rerehe k','Bias rerehe heri', 'Var rerehe heri', 'MSE rerehe heri' #9
    )
    
    
    
    # this records CPU time, which adds up system.time first two entries

    table_time = (matrix(NA, nrow=4*5*4*2*2, ncol= 7))
    colnames(table_time) = c('n', 'e2','k2','condition','ratio','mean1median2', 
                             'reREHE fit' # only point estimation time
    )
    
    
    rehe_full_est = list(NULL)
    
    
    
    counter = 0
    
    for(ww in 1:160){
      counter=counter+1 
      
      try({
        cat(ww, '\n')
        
        filepath = '...'
        input = read.table(paste0(filepath, '/input_simulation_reREHE.txt'), header=T)
        input_row=ww
        
        argv = as.matrix(input[as.integer(input_row),])
        nK=1
        
        n = as.integer(argv[1])
        sig.e2 = as.numeric(argv[2])
        sig.k2 = as.numeric(argv[3])
        condition = as.integer(argv[4]) 
        ratio = argv[5]
        
        true.heri = sig.k2 / sum(sig.e2, sig.k2)
        
        table_point_est[counter,1:5]<-table_time[counter, 1:5]<-argv
        table_point_est[counter+160,1:5]<-table_time[counter+160, 1:5]<-argv
        
        load(paste0(filepath, "output/" ,savename,paste(argv, collapse='_'), '_ratio_', ratio, "_condition_", condition, ".RData"))
        
        gather = function(tmp){
          if(class(tmp)=='list'){
            do.call(cbind, tmp)
          }else{
            tmp
          }
        }
        
        
        REHE =  list(
          e_k_median = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_sub)),
          e_k_mean = gather(sapply(1:nreps, function(i) res[[i]]$rehe$vc_rehe_sub_mean)),
                     
          heri_median = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_sub)),
          heri_mean = gather(sapply(1:nreps, function(i) res[[i]]$rehe$heri_rehe_sub_mean))
                     
        )   
        
        
        rehe_full_est[[counter]] = REHE

        
        summ = function(xx, true){
          tmp = rep(0, 3)
          if(class(xx)=='matrix'){
            for(i in 1:nrow(xx)){
              x = as.vector(xx[i,])
              tmp = rbind(tmp, c(mean(x)-true[i], var(x), (mean(x)-true[i])^2+var(x)))
            }
            return(t(tmp[-1,]))
          }else{
            return(c(mean(xx)-true[1], var(xx), (mean(xx)-true[1])^2+var(xx)))
          }
        }
        
        table_point_est[counter, c(6:15)]<-c(
          2, 
          c(summ(REHE$e_k_median, c(sig.e2, sig.k2))),
          summ(REHE$heri_median, true.heri)
        )
        
        table_point_est[counter+160, c(6:15)]<-c(
          1, 
          c(summ(REHE$e_k_mean, c(sig.e2, sig.k2))),
          summ(REHE$heri_mean, true.heri)
        )
        
        
        table_time[counter, 6:7]<-c(12,  median(gather(sapply(1:nreps, function(i)sum(res[[i]]$time$sub[1:2]))))
        )
        
        

      })
    }
    
    
    reREHE_table_time = data.frame(table_time[1:160, ])
    # reREHE_table_time
    reREHE_table_point_est = data.frame(table_point_est)
    rehe_sub_est = rehe_full_est

    rm('In')
    
    save.image('new_reREHE_wrap_up_summary.RData')
    
  }
  
  load('new_reREHE_wrap_up_summary.RData')
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')

 
}





 
  
  
  
  
  