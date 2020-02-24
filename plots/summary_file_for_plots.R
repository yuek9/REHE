## summarize wrap up simulation results
library(ggplot2)
library(ggpubr)
library(scales)
library(cowplot)

filepath = '...'
setwd(paste0(filepath, '/wrap_up_files'))

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


########
# wrap up simulation results into tables
########
if(T){
  # results for simulation with sparsification when constructing confidence interval
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
    
    
    
    savename = paste0('sparse_parametric_bootstrap_')
    counter = 0
    for(ww in 1:80){
      counter=counter+1 
      
      try({
        cat(ww, '\n')
        
        filepath = '\\\\gcc-fs4\\gcc-fs4\\thornton\\Kun_Yue_HCHS_SOL'
        input = read.table(paste0(filepath, '/Code/wrap_up_files/input_simulation.txt'), header=T)
        input_row=ww
        
        argv = as.matrix(input[as.integer(input_row),])
        nK=1
        
        n = as.integer(argv[1])
        sig.e2 = as.numeric(argv[2])
        sig.k2 = as.numeric(argv[3])
        condition = as.integer(argv[4]) 
        
        true.heri = sig.k2 / sum(sig.e2, sig.k2)
        
        table_CI[counter,1:4]<-table_point_est[counter,1:4]<-table_time[counter, 1:4]<-argv
        
        load(paste0(filepath, "/Code/wrap_up_files/" ,savename,paste(argv, collapse='_'), "_condition_", condition, ".RData"))
        
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
    
    
    
    savename = paste0('dense_parametric_bootstrap_')
    counter = 0
    for(ww in 1:80){
      counter=counter+1 
      
      try({
        cat(ww, '\n')
        
        filepath = '\\\\gcc-fs4\\gcc-fs4\\thornton\\Kun_Yue_HCHS_SOL'
        input = read.table(paste0(filepath, '/Code/wrap_up_files/input_simulation.txt'), header=T)
        input_row=ww
        
        argv = as.matrix(input[as.integer(input_row),])
        nK=1
        
        n = as.integer(argv[1])
        sig.e2 = as.numeric(argv[2])
        sig.k2 = as.numeric(argv[3])
        condition = as.integer(argv[4]) 
        
        true.heri = sig.k2 / sum(sig.e2, sig.k2)
        
        table_CI[counter,1:4]<-table_point_est[counter,1:4]<-table_time[counter, 1:4]<-argv
        
        load(paste0(filepath, "/Code/wrap_up_files/" ,savename,paste(argv, collapse='_'), "_condition_", condition, ".RData"))
        
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
  
}

load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
load('new_wrap_up_summary.RData')



##############################
##results for condition 4 only (for main paper)##
##############################
if(T){

  
  #time
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  {
    con4 = seq(4, 80, 4)
    table_time = table_time[con4,]
    dense_table_time = dense_table_time[con4,]
    
    #################time
    ## time for sparsification
    { 
      tmp = round(table_time,3)
      tmp = data.frame(table_time)
      tmp$setting = rep(1:5, 4)
      # time under same sample size are very similar, so might just average over same sample size and show time for different methods
      time_agg = aggregate(tmp,by=list(tmp$n),function(x)mean(x, na.rm=T))[-c(1, 3, 4, 5, 12)]
      colnames(time_agg)<-c('n', 'sREHE CI', 'sREHE est', 'REHE est', 'sREML est', 'REML', 'Sparsification')
      time_agg$sREML = time_agg$`sREML est`+time_agg$Sparsification
      
      time_agg = reshape(time_agg, varying =list(names(time_agg)[2:8]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[2:8], timevar = 'method')
      
      # png('time_agg_method4.png',width=600, height=300, unit='px')
      # ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
      #   geom_line(size=2)
      # dev.off()
      
      time_agg$method = factor(time_agg$method, levels = c('REML', 'sREML', 'sREML est', 'REHE est', 'sREHE est' , 'sREHE CI', 'Sparsification' ))
      
      # time_plot = ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
      #   geom_line(size=2)
      
      
      time_plot = ggplot(time_agg[! time_agg$method %in% c('sREML', 'sREHE est'),],aes(x=n, y=log10(time), group=method, color=method))+
        geom_line(size=1.7, aes(linetype = method))+
        geom_point(size=2)+
        scale_linetype_manual(values=c(1, 3, 1, 3, 4))+
        scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 6)])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
    }
    
    
    ## time for dense only
    {
      tmp = data.frame(dense_table_time)
      tmp$setting = rep(1:5, 4)
      # time under same sample size are very similar, so might just average over same sample size and show time for different methods
      time_agg = aggregate(tmp,by=list(tmp$n),function(x)mean(x, na.rm=T))[-c(1, 3, 4, 5, 7, 9, 11,12)]
      colnames(time_agg)<-c('n', 'REHE CI', 'REHE est',  'REML')
      time_agg = reshape(time_agg, varying =list(names(time_agg)[2:ncol(time_agg)]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[2:ncol(time_agg)], timevar = 'method')
      time_agg$method = factor(time_agg$method, levels = c('REML', 'REHE est', 'REHE CI'))
      
      dense_time_plot = ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
        geom_line(size=1.7, aes(linetype = method))+
        geom_point(size=2)+
        scale_linetype_manual(values=c(1,  1, 3))+
        scale_color_manual(values=cbPalette[c(2, 3, 3 )])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
    }
    
    
  }
  
  # time_plot+ theme_bw()
  # dense_time_plot
  ##################
  #MSE
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  {
    # for condition 1 and 2, sparse and nonsparse will be the same
    table_point_est = data.frame(table_point_est[con4,])
    table_point_est$setting = rep(1:5, 4)
    
    # plot for: 7=e, 10=k, 13=heri
    
    
    start=13
    
    condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
    colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
    
    condition34 = reshape(condition34, varying =list(names(condition34)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='MSE', times=names(condition34)[4:9], timevar = 'method')
    
    
    
    # png(paste0(start,'point_condition4.png'),width=1250, height=800, unit='px')
    # ggplot(condition34, aes(x=n, y=MSE, group=method, color=method, linetype=method))+
    #   geom_line(size=1.5)+
    #   scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    #   facet_wrap(.~condition+setting, nrow=2, scales='free',
    #              labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    # 
    # ggplot(condition34[!condition34$method%in%c('HE-sparse', 'REHE-sparse', 'sREML'),], aes(x=n, y=MSE, group=method, color=method, linetype=method))+
    #   geom_line(size=1.5)+
    #   scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    #   facet_wrap(.~condition+setting, nrow=2, scales='free',
    #              labeller = labeller(condition = label_condition, setting = label_setting))
    
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse'))
    
    
    
    MSE_plot_h = 
      ggplot(condition34[condition34$setting==1& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('MSE')+ggtitle(label ='Heritability')+
      labs(caption = bquote(~ sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 4))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 4)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    dense_MSE_plot_h = 
      ggplot(condition34[condition34$setting==1& !(condition34$method %in% c('HE-sparse', 'REHE-sparse', 'sREML')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('MSE')+ggtitle(label ='Heritability')+
      labs(caption = bquote(~ sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 1, 4))+
      scale_color_manual(values=cbPalette[c(2,  3, 4)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    start=10
    
    condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
    colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
    
    condition34 = reshape(condition34, varying =list(names(condition34)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='MSE', times=names(condition34)[4:9], timevar = 'method')
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse'))
    
    
    
    MSE_plot_v = 
      ggplot(condition34[condition34$setting==4& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('MSE')+
      ggtitle(label =bquote(sigma[1]^2 ))+ 
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 4))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 4)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    
    dense_MSE_plot_v = 
      ggplot(condition34[condition34$setting==4& !(condition34$method %in% c('HE-sparse', 'REHE-sparse', 'sREML')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('MSE')+
      ggtitle(label =bquote(sigma[1]^2 ))+ 
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1,  1, 4))+
      scale_color_manual(values=cbPalette[c(2,  3, 4)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    
    
    # more for supplement
    {
      start=7
      
      condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
      colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
      
      condition34 = reshape(condition34, varying =list(names(condition34)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                            v.names='MSE', times=names(condition34)[4:9], timevar = 'method')
      condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse'))
      
      MSE_main_sup1 = ggplot(condition34[condition34$setting==5& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
        geom_line(size=1.5)+
        ylab('MSE')+
        ggtitle(label =bquote(sigma[0]^2 ))+ 
        labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
        geom_point(size=2)+
        scale_linetype_manual(values=c(1, 3, 1, 4))+
        scale_color_manual(values=cbPalette[c(2, 2, 3, 4)])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
        scale_y_continuous(labels=scientific)
      
      
      start=10
      
      condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
      colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
      
      condition34 = reshape(condition34, varying =list(names(condition34)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                            v.names='MSE', times=names(condition34)[4:9], timevar = 'method')
      condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse'))
      
      MSE_main_sup2 = ggplot(condition34[condition34$setting==5& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=(MSE), group=method, color=method, linetype=method))+
        geom_line(size=1.5)+
        ylab('MSE')+
        ggtitle(label =bquote(sigma[1]^2 ))+ 
        labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
        geom_point(size=2)+
        scale_linetype_manual(values=c(1, 3, 1, 4))+
        scale_color_manual(values=cbPalette[c(2, 2, 3, 4)])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
        scale_y_continuous(labels=scientific)
      
      }
    
    
    
  }
  # MSE_plot_h
  # MSE_plot_v
  # dense_MSE_plot_h
  # dense_MSE_plot_v
  
  
  ##### comparison of confidence interval coverage
  
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  # for sparsification
  {
    
    table_CI = data.frame(table_CI[con4,])
    table_CI$setting = rep(1:5, 4)
    CI_NA = table_CI[, c(17,21)]
    
    # cbind(table_CI[(CI_NA[,1]>0),c(1,2,3,4,20)], CI_NA[CI_NA[,1]>0,])
    
    
    table_CI = table_CI[, -c(17,21)]
    table_CI[,c(1,4,20, seq(5, 19, by=3))]
    
    
    start=7
    
    condition34 = table_CI[table_CI$condition%in%c(3, 4),c(1,4,20, seq(start, 19, by=3))]
    condition34 = cbind(condition34,   dense_table_CI[con4,-c(17,21)][, seq(start, 19, by=3)][,c(1, 3)])
    colnames(condition34)<-c( "n" ,"condition" ,"setting" , "sREHE Wald" ,    "REHE-spase Wald", "sREHE Quantile" ,"REML", "sREML" , 'REHE Wald', 'REHE Quantile')
    
    # do not show REHE-sparse CI, not useful
    condition34 = condition34[,-5]
    condition34 = reshape(condition34, varying =list(names(condition34)[4:ncol(condition34)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='Coverage', times=names(condition34)[4:ncol(condition34)], timevar = 'method')
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    # png(paste0(start,'CI_condition4.png'),width=1250, height=800, unit='px')
    # ggplot(condition34, aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
    #   geom_line(size=1.5)+
    #   geom_hline(yintercept = 0.95,linetype=10, color=1)+
    #   facet_wrap(.~condition+setting, nrow=2, scales='free',
    #              labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    
    
    
    CIcov_plot1 = ggplot(condition34[condition34$setting==1,], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      ylab('CI coverage')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 4, 4)])
    
    
    start=6
    condition34 = table_CI[table_CI$condition%in%c(3, 4),c(1,4,20, seq(start, 19, by=3))]
    condition34 = cbind(condition34,   dense_table_CI[con4,-c(17,21)][, seq(start, 19, by=3)][,c(1, 3)])
    colnames(condition34)<-c( "n" ,"condition" ,"setting" , "sREHE Wald" ,    "REHE-spase Wald", "sREHE Quantile" ,"REML", "sREML" , 'REHE Wald', 'REHE Quantile')
    
    # do not show REHE-sparse CI, not useful
    condition34 = condition34[,-5]
    condition34 = reshape(condition34, varying =list(names(condition34)[4:ncol(condition34)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='Coverage', times=names(condition34)[4:ncol(condition34)], timevar = 'method')
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    CIcov_plot2 = ggplot(condition34[condition34$setting==4,], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      ylab('CI coverage')+
      ggtitle(label=bquote( sigma[1]^2 ))+
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 4, 4)])
    
  }
  
  
  # for dense
  {
    con4 = seq(4, 80, 4)
    
    dense_table_CI = data.frame(dense_table_CI[con4,])
    dense_table_CI$setting = rep(1:5, 4)
    CI_NA = dense_table_CI[, c(17,21)]
    
    
    dense_table_CI = dense_table_CI[, -c(17,21)]
    
    start=7
    
    condition34 = dense_table_CI[dense_table_CI$condition%in%c(3, 4),c(1,4,20, seq(start, 19, by=3))]
    colnames(condition34)<-c( "n" ,"condition" ,"setting" , "REHE Wald" ,    "REHE-spase Wald", "REHE Quantile" ,"REML", "sREML" )
    
    # do not show REHE-sparse CI, not useful
    condition34 = condition34[,-5]
    condition34 = reshape(condition34, varying =list(names(condition34)[4:7]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='Coverage', times=names(condition34)[4:7], timevar = 'method')
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    # png(paste0(start,'CI_condition4.png'),width=1250, height=800, unit='px')
    # ggplot(condition34, aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
    #   geom_line(size=1.5)+
    #   geom_hline(yintercept = 0.95,linetype=10, color=1)+
    #   facet_wrap(.~condition+setting, nrow=2, scales='free',
    #              labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    # 
    # ggplot(condition34[condition34$method!='sREML',], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
    #   geom_line(size=1.5)+
    #   geom_hline(yintercept = 0.95,linetype=10, color=1)+
    #   facet_wrap(.~condition+setting, nrow=2, scales='free',
    #              labeller = labeller(condition = label_condition, setting = label_setting))
    
    
    dense_CIcov_plot1 = ggplot(condition34[condition34$setting==1&condition34$method!='sREML',], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      ylab('CI coverage')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1,  1, 3))+
      scale_color_manual(values=cbPalette[c(2,  3, 3)])
    
    
    start=6
    condition34 = dense_table_CI[dense_table_CI$condition%in%c(3, 4),c(1,4,20, seq(start, 19, by=3))]
    colnames(condition34)<-c( "n" ,"condition" ,"setting" , "REHE Wald" ,    "REHE-spase Wald", "REHE Quantile" ,"REML", "sREML" )
    
    # do not show REHE-sparse CI, not useful
    condition34 = condition34[,-5]
    condition34 = reshape(condition34, varying =list(names(condition34)[4:7]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='Coverage', times=names(condition34)[4:7], timevar = 'method')
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    dense_CIcov_plot2 = ggplot(condition34[condition34$setting==4&condition34$method!='sREML',], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      ylab('CI coverage')+
      ggtitle(label=bquote( sigma[1]^2 ))+
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1,  1, 3))+
      scale_color_manual(values=cbPalette[c(2,  3, 3)])
    
    # ggarrange(dense_CIcov_plot1, dense_CIcov_plot2, nrow=1, common.legend = T, legend='bottom')
    
  }
  
  # dense_CIcov_plot1
  # dense_CIcov_plot2
  
  #### confidence interval width comparison
  
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  # for sparsification;
  {
    table_CI = data.frame(table_CI[con4,])
    table_CI$setting = rep(1:5, 4)
    rehe_full_est =  rehe_full_est[con4]
    reml_est = reml_est[con4]
    dense_table_CI = data.frame(dense_table_CI[con4,])
    dense_table_CI$setting = rep(1:5, 4)
    dense_rehe_full_est =  dense_rehe_full_est[con4]
    dense_reml_est = dense_reml_est[con4]
    
    table_CI_width = matrix(NA, nrow=1, ncol = ncol(table_CI))
    
    tmp_plot = list(NULL)
    
    
    target = sig.k2
    target_name = 'k2' 
    plot_data = data.frame(width=NULL, method=NULL, condition=NULL, setting=NULL)
    
    get_plt_data =  function(target_name, i){   
      sig.e2 = table_CI[i,2]
      sig.k2 = table_CI[i,3]
      sig.heri = round(sig.k2/(sig.e2+sig.k2), 3)
      print(c(sig.e2, sig.k2, sig.heri))   
      
      if(is.null(rehe_full_est[[i]])){
        return(NULL)
      }
      
      
      if(target_name=='e2'){
        target = sig.e2
        tmp = rehe_full_est[[i]]$e_CI
        tmpq = rehe_full_est[[i]]$e_CI_q
        tmp2 = t(gather(reml_est[[i]]$e_CI))
        tmp3 = t(gather(reml_est[[i]]$e_CI_sparse))
      }else   if(target_name=='k2'){
        # 
        # 
        target=sig.k2
        tmp = rehe_full_est[[i]]$k_CI
        tmpq = rehe_full_est[[i]]$k_CI_q
        tmp2 = t(gather( reml_est[[i]]$k_CI))
        tmp3 = t(gather(reml_est[[i]]$k_CI_sparse))
      }else  if(target_name=='heri'){  
        target=sig.heri
        tmp = rehe_full_est[[i]]$heri_CI
        tmpq = rehe_full_est[[i]]$heri_CI_q
        tmp2 = t(gather(reml_est[[i]]$heri_CI))
        tmp3 = t(gather(reml_est[[i]]$heri_CI_sparse))
      }  else{
        print('error')
      }   
      
      plot_data = data.frame(width = c((tmp[2,]-tmp[1,])/2, (tmpq[2,]-tmpq[1,])/2, (gather(tmp2[,2])-gather(tmp2[,1]))/2, (gather(tmp3[,2])-gather(tmp3[,1]))/2), 
                             method = c(rep('sREHE Wald', ncol(tmp)),rep('sREHE Quantile',ncol(tmpq)), rep('REML', nrow(tmp2)), rep('sREML', nrow(tmp3))),
                             n=table_CI[i,1],
                             condition = table_CI[i,4],
                             setting = table_CI$setting[i])
      # coverage = c(mean(check_coverage(t(tmp), rep(target, ncol(tmp))), na.rm=T), 
      #              mean(check_coverage(t(tmpq), rep(target, ncol(tmpq))), na.rm=T),
      #              mean(check_coverage(tmp2, rep(target, nrow(tmp2))), na.rm=T), 
      #              mean(check_coverage(tmp3, rep(target, nrow(tmp3))), na.rm=T))
      # plt1 = ggplot(plot_data, aes(x=method, y=width, color=method))+
      #   geom_boxplot()+
      #   ggtitle(paste0('half CI width: ', target_name, '=',target))+
      #   geom_text(data=data.frame(), aes(x=1:length(coverage), y=rep( range(plot_data$width, na.rm=T)[1], length(coverage)),label=paste0('CI:', round(coverage, 3)), col='red', size=3))+
      #   theme(legend.position = 'none')
      
      return(plot_data)
    }
    
    get_plt_data_dense =  function(target_name, i){   
      sig.e2 = dense_table_CI[i,2]
      sig.k2 = dense_table_CI[i,3]
      sig.heri = round(sig.k2/(sig.e2+sig.k2), 3)
      print(c(sig.e2, sig.k2, sig.heri))   
      
      if(is.null(dense_rehe_full_est[[i]])){
        return(NULL)
      }
      
      
      if(target_name=='e2'){
        target = sig.e2
        tmp = dense_rehe_full_est[[i]]$e_CI
        tmpq = dense_rehe_full_est[[i]]$e_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$e_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$e_CI_sparse))
      }else   if(target_name=='k2'){
        # 
        # 
        target=sig.k2
        tmp = dense_rehe_full_est[[i]]$k_CI
        tmpq = dense_rehe_full_est[[i]]$k_CI_q
        tmp2 = t(gather( dense_reml_est[[i]]$k_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$k_CI_sparse))
      }else  if(target_name=='heri'){  
        target=sig.heri
        tmp = dense_rehe_full_est[[i]]$heri_CI
        tmpq = dense_rehe_full_est[[i]]$heri_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$heri_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$heri_CI_sparse))
      }  else{
        print('error')
      }   
      
      plot_data = data.frame(width = c((tmp[2,]-tmp[1,])/2, (tmpq[2,]-tmpq[1,])/2, (gather(tmp2[,2])-gather(tmp2[,1]))/2, (gather(tmp3[,2])-gather(tmp3[,1]))/2), 
                             method = c(rep('REHE Wald', ncol(tmp)),rep('REHE Quantile',ncol(tmpq)), rep('REML', nrow(tmp2)), rep('sREML', nrow(tmp3))),
                             n=dense_table_CI[i,1],
                             condition = dense_table_CI[i,4],
                             setting = dense_table_CI$setting[i])
      
      
      return(plot_data)
    }
    
    # tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
    #   print(i)
    #   get_plt_data('e2', i)}, error = function(e) NULL))
    # tmp = do.call(rbind, tmp)
    # tmp$n = as.factor(tmp$n)
    # 
    # png(paste0('e2CI_width4.png'),width=1500, height=1250, unit='px')
    # ggplot(tmp, aes(x=method, y=width, fill=n))+
    #   geom_boxplot()+
    #   facet_wrap(.~setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    # 
    # 
    # tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
    #   print(i)
    #   get_plt_data('k2', i)}, error = function(e) NULL))
    # tmp = do.call(rbind, tmp)
    # tmp$n = as.factor(tmp$n)
    # 
    # png(paste0('k2CI_width4.png'),width=1500, height=1250, unit='px')
    # ggplot(tmp, aes(x=method, y=width, fill=n))+
    #   geom_boxplot()+
    #   facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    # 
    # 
    # tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
    #   print(i)
    #   get_plt_data('heri', i)}, error = function(e) NULL))
    # tmp = do.call(rbind, tmp)
    # tmp$n = as.factor(tmp$n)
    # png(paste0('heriCI_width4.png'),width=1500, height=1250, unit='px')
    # ggplot(tmp, aes(x=method, y=width, fill=n))+
    #   geom_boxplot()+
    #   facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    
    
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('heri', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'sREHE Wald', 'sREHE Quantile'))
    
    tmp2 = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data_dense('heri', i)}, error = function(e) NULL))
    tmp2 = do.call(rbind, tmp2)
    tmp2$n = as.factor(tmp2$n)
    tmp2$method = factor(tmp2$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    tmp = rbind(tmp, tmp2[tmp2$method%in%c( 'REHE Wald', 'REHE Quantile'),])
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML','REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    
    # ggplot(tmp, aes(x=n, y=width, fill=method))+
    #   geom_boxplot()+facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    
    
    CIw_plot1  = 
      ggplot(tmp[tmp$n==12000&tmp$setting==1,], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent',cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
    
    tmp2 = aggregate(data=tmp, width~method+n+condition+setting, median)
    
    
    CIw_plot1_m = ggplot(tmp2[tmp2$setting==1, ], aes(x=n, y=width, group=method, color=method,linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 4, 4)])
    
    
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('k2', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'sREHE Wald', 'sREHE Quantile'))
    tmp2 = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data_dense('k2', i)}, error = function(e) NULL))
    tmp2 = do.call(rbind, tmp2)
    tmp2$n = as.factor(tmp2$n)
    tmp2$method = factor(tmp2$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    tmp = rbind(tmp, tmp2[tmp2$method%in%c( 'REHE Wald', 'REHE Quantile'),])
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML','REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    
    CIw_plot2   = 
      
      ggplot(tmp[tmp$n==12000&tmp$setting==4,], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption= bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent',cbPalette[3], 'transparent',cbPalette[4], 'transparent'))
    
    tmp2 = aggregate(data=tmp, width~method+n+condition+setting, median)
    
    
    CIw_plot2_m = ggplot(tmp2[tmp2$setting==4, ], aes(x=n, y=width, group=method, color=method,linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3, 1, 4))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 4, 4)])
    
    
    supp_CIw_con4 =     ggplot(tmp[tmp$n==12000&tmp$setting==5,], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption= bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent',cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
    
    
    supp_CIw_con4_m = ggplot(tmp2[tmp2$setting==5,], aes(x=n, y=width, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 4, 4)])
    
    
  }
  
  # for dense;  use line plots for main paper (half width median); pick some boxplot for supplementary
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  {
    dense_table_CI = data.frame(dense_table_CI[con4,])
    dense_table_CI$setting = rep(1:5, 4)
    dense_rehe_full_est =  dense_rehe_full_est[con4]
    dense_reml_est = dense_reml_est[con4]
    
    table_CI_width = matrix(NA, nrow=1, ncol = ncol(dense_table_CI))
    
    tmp_plot = list(NULL)
    
    
    target = sig.k2
    target_name = 'k2' 
    plot_data = data.frame(width=NULL, method=NULL, condition=NULL, setting=NULL)
    
    get_plt_data =  function(target_name, i){   
      sig.e2 = dense_table_CI[i,2]
      sig.k2 = dense_table_CI[i,3]
      sig.heri = round(sig.k2/(sig.e2+sig.k2), 3)
      print(c(sig.e2, sig.k2, sig.heri))   
      
      if(is.null(dense_rehe_full_est[[i]])){
        return(NULL)
      }
      
      
      if(target_name=='e2'){
        target = sig.e2
        tmp = dense_rehe_full_est[[i]]$e_CI
        tmpq = dense_rehe_full_est[[i]]$e_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$e_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$e_CI_sparse))
      }else   if(target_name=='k2'){
        # 
        # 
        target=sig.k2
        tmp = dense_rehe_full_est[[i]]$k_CI
        tmpq = dense_rehe_full_est[[i]]$k_CI_q
        tmp2 = t(gather( dense_reml_est[[i]]$k_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$k_CI_sparse))
      }else  if(target_name=='heri'){  
        target=sig.heri
        tmp = dense_rehe_full_est[[i]]$heri_CI
        tmpq = dense_rehe_full_est[[i]]$heri_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$heri_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$heri_CI_sparse))
      }  else{
        print('error')
      }   
      
      plot_data = data.frame(width = c((tmp[2,]-tmp[1,])/2, (tmpq[2,]-tmpq[1,])/2, (gather(tmp2[,2])-gather(tmp2[,1]))/2, (gather(tmp3[,2])-gather(tmp3[,1]))/2), 
                             method = c(rep('REHE Wald', ncol(tmp)),rep('REHE Quantile',ncol(tmpq)), rep('REML', nrow(tmp2)), rep('sREML', nrow(tmp3))),
                             n=dense_table_CI[i,1],
                             condition = dense_table_CI[i,4],
                             setting = dense_table_CI$setting[i])
      
      
      return(plot_data)
    }
    
    
    tmp = lapply(1:nrow(dense_table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('heri', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    # ggplot(tmp[tmp$method!='sREML',], aes(x=n, y=width, fill=method))+
    #   geom_boxplot()+facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    
    
    dense_CIw_plot1  = 
      ggplot(tmp[tmp$n==12000&tmp$setting==1&tmp$method!='sREML',], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      scale_color_manual(values=cbPalette[c(1, 1, 3)])+
      scale_fill_manual(values=c(cbPalette[2], cbPalette[3], 'transparent'))
    
    tmp2 = aggregate(data=tmp, width~method+n+condition+setting, median)
    
    
    dense_CIw_plot1_m = ggplot(tmp2[tmp2$setting==1&tmp2$method!='sREML', ], aes(x=n, y=width, group=method, color=method,linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label ='Heritability')+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1,  1, 3))+
      scale_color_manual(values=cbPalette[c(2,  3, 3)])
    
    
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('k2', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'REHE Wald', 'REHE Quantile'))
    
    dense_CIw_plot2   = 
      
      ggplot(tmp[tmp$n==12000&tmp$setting==4,], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption= bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1,  1, 3)])+
      scale_fill_manual(values=c(cbPalette[2], cbPalette[3], 'transparent'))
    
    tmp2 = aggregate(data=tmp, width~method+n+condition+setting, median)
    
    
    dense_CIw_plot2_m = ggplot(tmp2[tmp2$setting==4, ], aes(x=n, y=width, group=method, color=method,linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 3, 3)])
    
    
    supp_dense_CIw_con4 = ggplot(tmp[tmp$n==12000&tmp$setting==5,], aes(x=method, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption= bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent',cbPalette[3], 'transparent'))
    
    
    supp_dense_CIw_con4_m = ggplot(tmp2[tmp2$setting==5,], aes(x=n, y=width, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      # facet_wrap(.~condition+setting, nrow=2, scales='free',
      #          labeller = labeller(condition = label_condition, setting = label_setting))+    
      ylab('CI half width (median)')+
      ggtitle(label =bquote(sigma[1]^2))+
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 3)])
    
    
  }
  
  # dense_CIw_plot1_m
  # dense_CIw_plot2_m
  
  
}
  


#############################
## results of condition 1 2 3, used in appendix
################################
if(T){
  ## time for sparsification
  
    load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
    load('new_wrap_up_summary.RData')
    {
      tmp = data.frame(table_time)
      tmp$setting = rep(1:5, 4)
      # time under same sample size are very similar, so might just average over same sample size and show time for different methods
      time_agg = aggregate(tmp,by=list(tmp$n),function(x)min(x, na.rm=T))[-c(1, 3, 4, 5,7, 8, 12)]
      colnames(time_agg)<-c('n', 'sREHE CI', 'sREML',  'REML', 'sparsification')
      time_agg = reshape(time_agg, varying =list(names(time_agg)[2:ncol(time_agg)]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[2:ncol(time_agg)], timevar = 'method')
      time_agg$method = factor(time_agg$method, levels = c('REML', 'sREML', 'sREHE CI', 'sparsification'))
      
      sup_time_plot = ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
        geom_line(size=1.7, aes(linetype = method))+
        geom_point(size=2)+
        scale_linetype_manual(values=c(1, 3,  1, 3))+
        scale_color_manual(values=cbPalette[c(2,2, 3, 6 )])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
    }
  
  
  ##################MSE, one plot for sparse and nonsparse
    
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  {
  # for condition 1 and 2, sparse and nonsparse will be the same
  table_point_est = data.frame(table_point_est)
  table_point_est$setting = rep(rep(1:5, each=4), 4)
  
  # plot for: 7=e, 10=k, 13=heri
        start=7

        condition123 = table_point_est[table_point_est$condition%in%c(1,2,3),c(1:4, 59, seq(start,58, by=9))][,-c(2, 3)]
        colnames(condition123)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
        condition123 = reshape(condition123, varying =list(names(condition123)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                              v.names='MSE', times=names(condition123)[4:9], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML','sREML', 'REHE','REHE-sparse', 'HE', 'HE-sparse'))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        

        
        # png(paste0(start,'point_condition123.png'),width=1250, height=800, unit='px')
        # ggplot(condition123, aes(x=n, y=MSE, group=method, color=method, linetype = method))+
        #   geom_line(size=1.5)+
        #   scale_color_manual(values = c(2, 'red2', 3, 4))+
        #   scale_linetype_manual(values=c(1, 3, 2, 3))+
        #   facet_wrap(.~condition+setting, ncol=5, scales='free',
        #              labeller = labeller(condition = label_condition, setting = label_setting))
        # dev.off()
        # 
        
        
        
        start=7
        
        condition123 = table_point_est[table_point_est$condition%in%c(1,2,3),c(1:4, 59, seq(start,58, by=9))][,-c(2, 3)]
        colnames(condition123)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
        condition123 = reshape(condition123, varying =list(names(condition123)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[4:9], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML','sREML', 'REHE','REHE-sparse', 'HE', 'HE-sparse'))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {
        MSE11 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1' ))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        
        MSE12 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2' ))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE13 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 3' ))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        }
        
        # ggarrange(MSE11, MSE12, MSE13, labels=c('A', 'B', 'C'), nrow=1, common.legend=T, legend='bottom')
        
        
        start= 10
        
        condition123 = table_point_est[table_point_est$condition%in%c(1,2,3),c(1:4, 59, seq(start,58, by=9))][,-c(2, 3)]
        colnames(condition123)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
        condition123 = reshape(condition123, varying =list(names(condition123)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[4:9], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML','sREML', 'REHE','REHE-sparse', 'HE', 'HE-sparse'))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {
        MSE21 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[1]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE22 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[1]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE23 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote(sigma[1]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 3'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        }
        
        # ggarrange(MSE21, MSE22, MSE23, labels=c('D', 'E', 'F'), nrow=1, common.legend=T, legend='bottom')
        
        
        
        start= 13
        
        condition123 = table_point_est[table_point_est$condition%in%c(1,2,3),c(1:4, 59, seq(start,58, by=9))][,-c(2, 3)]
        colnames(condition123)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
        condition123 = reshape(condition123, varying =list(names(condition123)[4:9]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[4:9], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML','sREML', 'REHE','REHE-sparse', 'HE', 'HE-sparse'))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {        
          MSE31 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        
        MSE32 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE33 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=MSE, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('MSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4)])+
          scale_linetype_manual(values=c(1, 3, 2, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 3'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
}       
        
        # ggarrange(MSE31, MSE32, MSE33, common.legend = T, nrow=1, legend='bottom', labels = c('G', 'H', 'I'))
  }




  ##### comparison of confidence interval coverage, dense and sparse together
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  {

  table_CI = data.frame(table_CI)
  dense_table_CI = data.frame(dense_table_CI)
  table_CI$setting = rep(rep(1:5, each=4), 4)
  CI_NA = table_CI[, c(17,21)]
  
  # cbind(table_CI[(CI_NA[,1]>0),c(1,2,3,4,20)], CI_NA[CI_NA[,1]>0,])
  
  
  table_CI = table_CI[, -c(17,21)]
  table_CI[,c(1,4,20, seq(5, 19, by=3))]
  dense_table_CI[,c(1,4, seq(5, 19, by=3))]
  # start=5
  # 
  # 
  # 
  # start=6
  # start=7
  # 
  # 
  #       condition123= table_CI[table_CI$condition%in%c(1,2,3),c(1,4,20, seq(start, 19, by=3))]
  #       colnames(condition123)<-c( "n" ,"condition" ,"setting" , "REHE Wald" ,    "REHE-sparse Wald", "REHE Quantile" ,"REML", "sREML" )
  #       condition123 = reshape(condition123, varying =list(names(condition123)[4:8]), direction = 'long', idvar=c('n', 'condition','setting'), 
  #                             v.names='Coverage', times=names(condition123)[4:8], timevar = 'method')
  #       condition123$method = factor(condition123$method, levels= c("REML", "sREML" ,"REHE Wald" , "REHE Quantile",    "REHE-sparse Wald"))
  # 
  #       condition123 = condition123[!condition123$method%in% c('REHE-sparse Wald'),]
  #       
  #       png(paste0(start,'CI_condition123.png'),width=1250, height=800, unit='px')
  #       ggplot(condition123, aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
  #         geom_line(size=1.5)+
  #         geom_hline(yintercept = 0.95,linetype=10, color=1)+
  #         facet_wrap(.~condition+setting, ncol=5, scales='free',
  #                    labeller = labeller(condition = label_condition, setting = label_setting))
  #       dev.off()
  #       
  #       tmps= condition123
  #       ggplot(tmps[tmps$setting==5 & tmps$condition==3,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
  #         geom_line(size=1.5)+
  #         geom_hline(yintercept = 0.95,linetype=10, color=1)
        
        
        
        
        
        start = 5
        condition123= table_CI[table_CI$condition%in%c(1,2,3),c(1,4,20, seq(start, 19, by=3))]
        condition123 = cbind(condition123,   dense_table_CI[dense_table_CI$condition%in%c(1,2,3),-c(17,21)][, seq(start, 19, by=3)][,c(1, 3)])
        colnames(condition123)<-c( "n" ,"condition" ,"setting" , "sREHE Wald" ,    "REHE-spase Wald", "sREHE Quantile" ,"REML", "sREML" , 'REHE Wald', 'REHE Quantile')
        condition123 = condition123[,-5]
        condition123 = reshape(condition123, varying =list(names(condition123)[4:ncol(condition123)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                              v.names='Coverage', times=names(condition123)[4:ncol(condition123)], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
        

        {
        CI11 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('CI coverage')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
          scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          geom_hline(yintercept = 0.95,linetype=10, color=1)
        
        
        CI12 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('CI coverage')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
          scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          geom_hline(yintercept = 0.95,linetype=10, color=1)
        
        CI13 = ggplot(condition123[condition123$condition==3 & condition123$setting==5,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('CI coverage')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
          scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          geom_hline(yintercept = 0.95,linetype=10, color=1)
        }
        
   
        
        
        
        start = 6
        condition123= table_CI[table_CI$condition%in%c(1,2,3),c(1,4,20, seq(start, 19, by=3))]
        condition123 = cbind(condition123,   dense_table_CI[dense_table_CI$condition%in%c(1,2,3),-c(17,21)][, seq(start, 19, by=3)][,c(1, 3)])
        colnames(condition123)<-c( "n" ,"condition" ,"setting" , "sREHE Wald" ,    "REHE-spase Wald", "sREHE Quantile" ,"REML", "sREML" , 'REHE Wald', 'REHE Quantile')
        condition123 = condition123[,-5]
        condition123 = reshape(condition123, varying =list(names(condition123)[4:ncol(condition123)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='Coverage', times=names(condition123)[4:ncol(condition123)], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
        
        {
          CI21 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            ggtitle(label =bquote(sigma[1]^2))+
            labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
          
          
          CI22 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            ggtitle(label =bquote(sigma[1]^2))+
            labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
          
          CI23 = ggplot(condition123[condition123$condition==3 & condition123$setting==5,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            ggtitle(label =bquote(sigma[1]^2))+
            labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
        }
        
 
        
        
        start = 7
        condition123= table_CI[table_CI$condition%in%c(1,2,3),c(1,4,20, seq(start, 19, by=3))]
        condition123 = cbind(condition123,   dense_table_CI[dense_table_CI$condition%in%c(1,2,3),-c(17,21)][, seq(start, 19, by=3)][,c(1, 3)])
        colnames(condition123)<-c( "n" ,"condition" ,"setting" , "sREHE Wald" ,    "REHE-spase Wald", "sREHE Quantile" ,"REML", "sREML" , 'REHE Wald', 'REHE Quantile')
        condition123 = condition123[,-5]
        condition123 = reshape(condition123, varying =list(names(condition123)[4:ncol(condition123)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='Coverage', times=names(condition123)[4:ncol(condition123)], timevar = 'method')
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
        
        {
          CI31 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            
            ggtitle(label =bquote('Heritability'))+
            labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
          
          
          CI32 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            
            ggtitle(label =bquote('Heritability'))+
            labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
          
          CI33 = ggplot(condition123[condition123$condition==3 & condition123$setting==5,], aes(x=n, y=Coverage, group=method, color=method, linetype = method))+
            geom_line(size=1.5)+ylab('CI coverage')+
            geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 3, 4, 4)])+
            scale_linetype_manual(values=c(1, 3, 1, 3, 1, 3))+
            
            ggtitle(label =bquote('Heritability'))+
            labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
            theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                  legend.title = element_blank()) +
            geom_hline(yintercept = 0.95,linetype=10, color=1)
        }
        
}

  
        

  #### confidence interval width comparison
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  
  {
    table_CI = data.frame(table_CI)   
    table_CI$setting = rep(1:5, each=4)
    dense_table_CI = data.frame(dense_table_CI)
    dense_table_CI$setting = rep(1:5, each=4)

    get_plt_data =  function(target_name, i){   
      sig.e2 = table_CI[i,2]
      sig.k2 = table_CI[i,3]
      sig.heri = round(sig.k2/(sig.e2+sig.k2), 3)
      print(c(sig.e2, sig.k2, sig.heri))   
      
      if(is.null(rehe_full_est[[i]])){
        return(NULL)
      }
      # round(table_6000_CI[i,],2)
      # 
      # 
      if(target_name=='e2'){
        target = sig.e2
        tmp = rehe_full_est[[i]]$e_CI
        tmpq = rehe_full_est[[i]]$e_CI_q
        tmp2 = t(gather(reml_est[[i]]$e_CI))
        tmp3 = t(gather(reml_est[[i]]$e_CI_sparse))
      }else   if(target_name=='k2'){
        # 
        # 
        target=sig.k2
        tmp = rehe_full_est[[i]]$k_CI
        tmpq = rehe_full_est[[i]]$k_CI_q
        tmp2 = t(gather( reml_est[[i]]$k_CI))
        tmp3 = t(gather(reml_est[[i]]$k_CI_sparse))
      }else  if(target_name=='heri'){  
        target=sig.heri
        tmp = rehe_full_est[[i]]$heri_CI
        tmpq = rehe_full_est[[i]]$heri_CI_q
        tmp2 = t(gather(reml_est[[i]]$heri_CI))
        tmp3 = t(gather(reml_est[[i]]$heri_CI_sparse))
      }  else{
        print('error')
      }   
      
      plot_data = data.frame(width = c((tmp[2,]-tmp[1,])/2, (tmpq[2,]-tmpq[1,])/2, (gather(tmp2[,2])-gather(tmp2[,1]))/2, (gather(tmp3[,2])-gather(tmp3[,1]))/2), 
                             method = c(rep('sREHE Wald', ncol(tmp)),rep('sREHE Quantile',ncol(tmpq)), rep('REML', nrow(tmp2)), rep('sREML', nrow(tmp3))),
                             n=table_CI[i,1],
                             condition = table_CI[i,4],
                             setting = table_CI[i,'setting'])
      # coverage = c(mean(check_coverage(t(tmp), rep(target, ncol(tmp))), na.rm=T), 
      #              mean(check_coverage(t(tmpq), rep(target, ncol(tmpq))), na.rm=T),
      #              mean(check_coverage(tmp2, rep(target, nrow(tmp2))), na.rm=T), 
      #              mean(check_coverage(tmp3, rep(target, nrow(tmp3))), na.rm=T))
      # plt1 = ggplot(plot_data, aes(x=method, y=width, color=method))+
      #   geom_boxplot()+
      #   ggtitle(paste0('half CI width: ', target_name, '=',target))+
      #   geom_text(data=data.frame(), aes(x=1:length(coverage), y=rep( range(plot_data$width, na.rm=T)[1], length(coverage)),label=paste0('CI:', round(coverage, 3)), col='red', size=3))+
      #   theme(legend.position = 'none')
      
      return(plot_data)
    }
    get_plt_data_dense =  function(target_name, i){   
      sig.e2 = dense_table_CI[i,2]
      sig.k2 = dense_table_CI[i,3]
      sig.heri = round(sig.k2/(sig.e2+sig.k2), 3)
      print(c(sig.e2, sig.k2, sig.heri))   
      
      if(is.null(dense_rehe_full_est[[i]])){
        return(NULL)
      }
      
      
      if(target_name=='e2'){
        target = sig.e2
        tmp = dense_rehe_full_est[[i]]$e_CI
        tmpq = dense_rehe_full_est[[i]]$e_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$e_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$e_CI_sparse))
      }else   if(target_name=='k2'){
        # 
        # 
        target=sig.k2
        tmp = dense_rehe_full_est[[i]]$k_CI
        tmpq = dense_rehe_full_est[[i]]$k_CI_q
        tmp2 = t(gather( dense_reml_est[[i]]$k_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$k_CI_sparse))
      }else  if(target_name=='heri'){  
        target=sig.heri
        tmp = dense_rehe_full_est[[i]]$heri_CI
        tmpq = dense_rehe_full_est[[i]]$heri_CI_q
        tmp2 = t(gather(dense_reml_est[[i]]$heri_CI))
        tmp3 = t(gather(dense_reml_est[[i]]$heri_CI_sparse))
      }  else{
        print('error')
      }   
      
      plot_data = data.frame(width = c((tmp[2,]-tmp[1,])/2, (tmpq[2,]-tmpq[1,])/2, (gather(tmp2[,2])-gather(tmp2[,1]))/2, (gather(tmp3[,2])-gather(tmp3[,1]))/2), 
                             method = c(rep('REHE Wald', ncol(tmp)),rep('REHE Quantile',ncol(tmpq)), rep('REML', nrow(tmp2)), rep('sREML', nrow(tmp3))),
                             n=dense_table_CI[i,1],
                             condition = dense_table_CI[i,4],
                             setting = dense_table_CI$setting[i])
      
      
      return(plot_data)
    }
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('e2', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp=tmp[tmp$condition!=4,]
    tmp2 = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data_dense('e2', i)}, error = function(e) NULL))
    tmp2 = do.call(rbind, tmp2)
    tmp2$n = as.factor(tmp2$n)

    tmp = rbind(tmp, tmp2[tmp2$method%in%c( 'REHE Wald', 'REHE Quantile'),])
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML','REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    
    
    #tmp[tmp$method%in%c('REML', 'sREML')&is.na(tmp$width),]
    lookREMLNA = with(tmp, aggregate(is.na(tmp$width)~n+method+condition+setting,  FUN=mean))
    lookREMLNA[lookREMLNA$`is.na(tmp$width)`>0,]
    # png(paste0('e2CI_width.png'),width=1500, height=1250, unit='px')
    #   ggplot(tmp, aes(x=n, y=width, fill=method))+
    #     geom_boxplot()+
    #     facet_wrap(.~condition+setting,ncol=5, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()

  
    
    {
    CIw11 = ggplot(tmp[tmp$setting==2& tmp$condition==1,], aes(x=n, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[0]^2))+
      labs(caption=bquote( sigma[0]^2  == '0.04'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
    
    CIw12 = ggplot(tmp[tmp$setting==4& tmp$condition==2,], aes(x=n, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[0]^2))+
      labs(caption=bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
    
    CIw13 = ggplot(tmp[tmp$setting==5& tmp$condition==3,], aes(x=n, y=width, fill=method, color=method))+
      geom_boxplot()+
      ylab('CI half width')+
      ggtitle(label =bquote(sigma[0]^2))+
      labs(caption=bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
      scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
    
    }
    
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('k2', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp=tmp[tmp$condition!=4,]
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'sREHE Wald', 'sREHE Quantile'))
    tmp2 = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data_dense('k2', i)}, error = function(e) NULL))
    tmp2 = do.call(rbind, tmp2)
    tmp2$n = as.factor(tmp2$n)
    
    tmp = rbind(tmp, tmp2[tmp2$method%in%c( 'REHE Wald', 'REHE Quantile'),])
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML','REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    {
      CIw21 = ggplot(tmp[tmp$setting==2& tmp$condition==1,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote(sigma[1]^2))+
        labs(caption=bquote( sigma[0]^2  == '0.04'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
      CIw22 = ggplot(tmp[tmp$setting==4& tmp$condition==2,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote(sigma[1]^2))+
        labs(caption=bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
      CIw23 = ggplot(tmp[tmp$setting==5& tmp$condition==3,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote(sigma[1]^2))+
        labs(caption=bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
    }

    
    
    # png(paste0('k2CI_width.png'),width=1500, height=1250, unit='px')
    # ggplot(tmp[tmp$n==12000,], aes(x=method, y=width, fill=n))+
    #   geom_boxplot()+
    #   facet_wrap(.~condition+setting, ncol=5,scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    
    
    tmp = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data('heri', i)}, error = function(e) NULL))
    tmp = do.call(rbind, tmp)
    tmp$n = as.factor(tmp$n)
    tmp=tmp[tmp$condition!=4,]
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML', 'sREHE Wald', 'sREHE Quantile'))
    tmp2 = lapply(1:nrow(table_CI), function(i)tryCatch({
      print(i)
      get_plt_data_dense('heri', i)}, error = function(e) NULL))
    tmp2 = do.call(rbind, tmp2)
    tmp2$n = as.factor(tmp2$n)
    
    tmp = rbind(tmp, tmp2[tmp2$method%in%c( 'REHE Wald', 'REHE Quantile'),])
    tmp$method = factor(tmp$method, levels = c('REML', 'sREML','REHE Wald', 'sREHE Wald','REHE Quantile', 'sREHE Quantile'))
    
    {
      CIw31 = ggplot(tmp[tmp$setting==2& tmp$condition==1,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote('Heritability'))+
        labs(caption=bquote( sigma[0]^2  == '0.04'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
      CIw32 = ggplot(tmp[tmp$setting==4& tmp$condition==2,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote('Heritability'))+
        labs(caption=bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
      CIw33 = ggplot(tmp[tmp$setting==5& tmp$condition==3,], aes(x=n, y=width, fill=method, color=method))+
        geom_boxplot()+
        ylab('CI half width')+
        ggtitle(label =bquote('Heritability'))+
        labs(caption=bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01, setting 3'))+
        theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        scale_color_manual(values=cbPalette[c(1, 2, 1, 3, 1, 4)])+
        scale_fill_manual(values=c(cbPalette[2], 'transparent', cbPalette[3], 'transparent', cbPalette[4], 'transparent'))
      
    }
    
    
    
    
    
    # png(paste0('heriCI_width.png'),width=1500, height=1250, unit='px')
    # ggplot(tmp[tmp$n==12000,], aes(x=method, y=width, fill=n))+
    #   geom_boxplot()+
    #   facet_wrap(.~condition+setting,ncol=5, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    # dev.off()
    

    
  }
  
   
  }
  
  
##### get the HE estimation zero proportion

load("new_wrap_up_summary.RData")

zero_prop = cbind(table_time[,1:4],t(sapply(1:80, function(i)rowMeans(he_est[[i]][['e_k']]==0))))

zero_prop[zero_prop[,4]==4,]


#######
# wrap up data application results
#######

if(T){
  
  load("heritability.RData")
  
  
  CIvc = round(cbind(varCompCI(nullMM, prop=F)[c(4,1:3),],nullMMREHE[['varComp']][c(4, 1:3)], nullMMREHE[['varCompCI']]$CIvc$wald, nullMMREHE[['varCompCI']]$CIvc$quantile
  ),4)
  CIheri = round(cbind(varCompCI(nullMM, prop=T)[c(4,1:3),],nullMMREHE[['varComp']][c(4, 1:3)]/sum(nullMMREHE[['varComp']]), nullMMREHE[['varCompCI']]$CIheri$wald, nullMMREHE[['varCompCI']]$CIheri$quantile
  ),4) 
  
  colnames(CIvc) <- c('REML', ' ', ' ', 'REHE', ' ',' ', ' ', ' ')
  rownames(CIvc) <- c('noise','kinship', 'community', 'household')
  
  CIvc = CIvc[,c(1,4,2:3,5:8)]
  
  dataVC = data.frame(CIvc[, 1:2], 'REML CI'= paste('[', paste(CIvc[,3], CIvc[,4], sep=','),']'),
                      'REHE Wald' =  paste('[', paste(CIvc[,5], CIvc[,6], sep=','),']'),
                      'REHE Quantile' =  paste('[', paste(CIvc[,7], CIvc[,8], sep=','),']'))
  
  
  
  library(kableExtra)
  # kable(data.frame(CIvc[, 1:2], 'REML CI'= paste('[', paste(CIvc[,3], CIvc[,4], sep=','),']'),
  #                  'REHE Wald' =  paste('[', paste(CIvc[,5], CIvc[,6], sep=','),']'),
  #                  'REHE Quantile' =  paste('[', paste(CIvc[,7], CIvc[,8], sep=','),']')),
  #                  booktab=T, 'latex')%>%
  #   kable_styling()%>%
  #   add_header_above(c(' '=1, 'Point Est'=2, '95 CI' = 3))
  
  
  
  colnames(CIheri) <- c('REML', ' ', ' ', 'REHE', ' ',' ', ' ', ' ')
  rownames(CIheri) <- c('noise','kinship', 'community', 'house')
  
  CIheri = CIheri[,c(1,4,2:3,5:8)]
  
  
  
  # kable(data.frame(CIheri[, 1:2], 'REML CI'= paste('[', paste(CIheri[,3], CIheri[,4], sep=','),']'),
  #                    'REHE Wald' =  paste('[', paste(CIheri[,5], CIheri[,6], sep=','),']'),
  #                    'REHE Quantile' =  paste('[', paste(CIheri[,7], CIheri[,8], sep=','),']')),
  #         booktab=T, 'latex')%>%
  #     kable_styling()%>%
  #     add_header_above(c(' '=1, 'Point Est'=2, '95 CI' = 3))
  
  colnames(CIheri)<-c('REML', 'REHE', 'REML.l', 'REML.u', 'Wald.l', 'Wald.u', 'Q.l', 'Q.u')
  CIheri = data.frame(CIheri)
  CIheri$method = rownames(CIheri)
  
  CIheri = as.matrix(CIheri)
  ggCIheri = data.frame(est = as.numeric(as.vector(CIheri[,c(1,2,2)])),method =rep(c('REML', 'REHE Wald', 'REHE Quantile'), each=4), CI.l = as.numeric(as.vector(CIheri[,c(3,5,7)])), CI.u = as.numeric(as.vector(CIheri[,c(4,6,8)])), target = rep(CIheri[,9], 3))
  
  ggCIheri$target2 = rep(c(-0.5,0, 0.5), each=4)
  
  ggCIheri$method = factor(ggCIheri$method, levels = c('REML', 'REHE Wald', 'REHE Quantile'))
  ggCIheri$target = factor(ggCIheri$target, levels = c('noise','kinship', 'community', 'house'))
  
  
  dataCI =
    ggplot(ggCIheri, aes(x=target2, y=est, colour=method, group=method)) + 
    geom_point(aes(x=target2, y=est, color=method, group=method))+
    geom_errorbar(aes(ymin=CI.l, ymax=CI.u, group=method, linetype = method), width=.1) +
    theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+
    facet_wrap(~target, scales= 'free', ncol=4)+
    theme(axis.title.y = element_blank())+
    scale_color_manual(values = cbPalette[c(2, 3, 3)])+
    scale_linetype_manual(values = c(1, 1, 2))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), legend.position = 'bottom', legend.title=element_blank())
  
  
  
  load("chr1REHE.RData")
  load("chr1.RData")
  
  
  
  assoc_part = assoc[assoc$Score.pval<=5e-8,]
  assoc_part = assoc_part[assoc_part$Score.pval<=5e-8,]  
  
  assocREHE_part = assocREHE[assocREHE$Score.pval<=5e-8,]
  
  all((assoc_part$variant.id %in% assocREHE_part$variant.id)) # REHE contains all genes picked by REML
  
  assoc_part = assoc[assoc$variant.id %in% c(assocREHE_part$variant.id,assoc_part$variant.id ),]
  assocREHE_part =  assocREHE[assocREHE$variant.id %in% c(assocREHE_part$variant.id,assoc_part$variant.id ),]
  
  assoc_part[assoc_part$Score.pval>5e-8,]
  
  
  all(assoc_part$variant.id == assocREHE_part$variant.id)
  
  datapval = ggplot(data=data.frame(), aes(x=assocREHE_part$Score.pval, y=assoc_part$Score.pval))+
    geom_point()+
    xlab('REHE p values')+
    ylab('REML p values')+
    geom_hline(aes(yintercept = 5e-8), color=2, linetype=2)+
    geom_vline(aes(xintercept = 5e-8), color=2, linetype=2)+
    geom_abline(aes(intercept=0, slope=1), color=2)
  
  
  
  
  tmp = data.frame(method=c('REML', 'REHE'), time = c(1433,134) )
  tmp$method = factor(tmp$method, levels=c('REML', 'REHE'))
  datatime = ggplot(tmp)+
    geom_col(
      aes(x=method, y=time,fill=method), width=0.3 )+
    coord_flip()+
    ylab('Time')+
    scale_fill_manual(values=c(2, 3))+theme(legend.position='none',axis.title.y=element_blank() )
  
  
  data_vcpoint = ggtexttable(data.frame(method= c('REML', 'REHE'), t(as.matrix(dataVC[,1:2]))), rows=NULL, cols = c(' ', 'noise','kinship', 'census', 'household'))
  data_vcpoint = ggtexttable(data.frame(VC = c('noise', 'kinship', 'census', 'household'), dataVC[,1:2]), rows=NULL, cols = c('','REML', 'REHE'))
  
}


###########################
### generate plots
###########################

## main fig1
aligned_1 = align_plots(dense_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'bottom', legend.title=element_blank())
                        , datapval+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                        dense_MSE_plot_h+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), 
                        dense_MSE_plot_v+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), align="hv", axis="tblr")

main1 = ggarrange(ggarrange(plotlist = aligned_1, ncol=2, nrow=2, labels=c('a', 'b', 'c', 'd')), 
                  get_legend(dense_MSE_plot_v+theme(legend.position = 'bottom', legend.title = element_blank())), nrow=2, heights=c(8, 1))

ggsave(filename = 'fig1_main1.pdf', plot=main1, device='pdf', width=7, height=7)



## main fig2
aligned_2 = align_plots(dense_CIcov_plot1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none')+
                          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)), 
                        dense_CIcov_plot2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none')+
                          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                        dense_CIw_plot1_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'),
                        dense_CIw_plot2_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'),
                        align="hv", axis="tblr")

main2 = ggarrange(ggarrange(plotlist = aligned_2, ncol=2, nrow=2, labels=c('a', 'b', 'c', 'd')), 
                  get_legend(dense_CIw_plot2_m+theme(legend.position = 'bottom', legend.title = element_blank())), nrow=2, heights=c(8, 1))



ggsave(filename = 'fig2_main2.pdf', plot=main2, device='pdf', width=7, height=7)


## sup fig2
sup_data_CI = dataCI
ggsave(filename = 'figS2_sup_data_CI.pdf', plot=sup_data_CI, device='pdf', width=5*1.3, height = 2.5*1.3)


## sup fig3
ggsave(filename='figS3_sup_time_plot.pdf', plot=sup_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title = element_blank()), device='pdf', width=3.5*1.8, height = 2.5*1.8)

## sup fig4
MSE_main_sup = ggarrange(MSE_main_sup1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                         MSE_main_sup2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                         MSE_plot_v+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)), 
                         labels = c('a', 'b', 'c'),nrow=1, common.legend = T, legend='bottom')
ggsave(filename='figS4_MSE_main_sup.pdf', plot=MSE_main_sup, device='pdf', width=7*1.5, height=3*1.5)          

## sup fig5
aligned_1 = align_plots(CIcov_plot1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none')+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                        CIcov_plot2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none')+ scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                        align="hv", axis="tblr")

CIcov_main_sup = ggarrange(ggarrange(plotlist = aligned_1,nrow=1, labels=c('a', 'b')), 
                           get_legend(CIcov_plot1+theme(legend.title=element_blank())), nrow=1, widths=c(5, 1))

ggsave(plot=CIcov_main_sup, filename='figS5_CIcov_main_sup.pdf', device='pdf', width=7*1.5, height=3.5*1.5)

## sup fig6
aligned_1 = align_plots(CIw_plot1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank())+theme(legend.title=element_blank(), legend.position = 'none'),
                        CIw_plot2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank())+theme(legend.title=element_blank(), legend.position = 'none'),
                        supp_CIw_con4+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),
                                                        axis.text.x=element_blank(),
                                                        axis.ticks.x=element_blank())+theme(legend.title=element_blank(), legend.position = 'none'),
                        CIw_plot1_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                        CIw_plot2_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                        supp_CIw_con4_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                        align="hv", axis="tblr")

CIw_main_sup = ggarrange(
  ggarrange(plotlist=aligned_1[1:3],  nrow=1, labels=c('a', 'b', 'c')),
  get_legend(CIw_plot1+theme(legend.title=element_blank(), legend.position = 'bottom')),
  ggarrange(plotlist=aligned_1[4:6], nrow=1, labels=c('d', 'e', 'f')),
  get_legend(CIw_plot1_m+theme(legend.title=element_blank(), legend.position = 'bottom')),
  nrow=4, heights=c(6, 1, 6, 1))
ggsave(plot=CIw_main_sup, filename='figS6_CIw_main_sup.pdf', device='pdf', width=7*1.5, height=7*1.5)


## sup fig7
aligned =  align_plots(MSE11+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE12+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE13+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE21+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE22+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE23+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE31+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE32+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       MSE33+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                       align="hv", axis="tblr")

suppfig1 = ggarrange(ggarrange(plotlist=aligned,  labels=c('a', 'b', 'c','d', 'e', 'f','g', 'h', 'i'), nrow=3, ncol=3),
                     get_legend(MSE11+theme(legend.title=element_blank(), legend.position = 'bottom')),
                     nrow=2, heights=c(18, 1))
ggsave(filename ='figS7_suppfigure1.pdf', plot=suppfig1,  device='pdf',  width=7*1.8, height=7*1.8 )


## sup fig8
aligned = align_plots(CI11+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI12+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI13+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI21+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI22+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI23+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI31+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI32+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CI33+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title=element_blank(), legend.position = 'none'),
                      align="hv", axis="tblr"
)
suppfig2 =  ggarrange(
  ggarrange(plotlist=aligned,  labels=c('a', 'b', 'c','d', 'e', 'f','g', 'h', 'i'), nrow=3, ncol=3),
  get_legend(CI11+theme(legend.title=element_blank(), legend.position = 'bottom')),
  nrow=2, heights=c(18, 1))
ggsave(filename ='figS8_suppfigure2.pdf', plot=suppfig2,  device='pdf',  width=7*1.8, height=7*1.8 )



## sup fig9
aligned = align_plots(CIw11+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw12+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw13+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw21+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw22+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw23+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw31+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw32+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      CIw33+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title=element_blank(), legend.position = 'none'),
                      align="hv", axis="tblr")


suppfig3 = ggarrange(ggarrange(plotlist=aligned,  labels=c('a', 'b', 'c','d', 'e', 'f','g', 'h', 'i'), nrow=3, ncol=3),
                     get_legend(CIw11+theme(legend.title=element_blank(), legend.position = 'bottom')),
                     nrow=2, heights=c(18, 1))
ggsave(filename ='figS9_suppfigure3.pdf', plot=suppfig3,  device='pdf',  width=7*1.8, height=7*1.8 )

 
  
  
  
  
  