## based on stored data mattrices to obtain the plots

library(ggplot2)
library(ggpubr)
library(scales)
library(cowplot)
filepath = '...'
setwd(filepath)
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
  
  
  ##############################
  ##results for condition 4 only##
  ##############################
  
  #time
  load('new_reREHE_wrap_up_summary.RData')
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')

  
  {
  con4 = seq(4, 80, 4)
  table_time = table_time[con4,]
  dense_table_time = dense_table_time[con4,]
  reREHE_table_time = reREHE_table_time[seq(4, 160, 4),]
  
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
  
  tmp = aggregate(reREHE_table_time, by=list(reREHE_table_time$n, reREHE_table_time$ratio), function(x)mean(x, na.rm=T))
  
  time_agg$`reREHE 0.05` = tmp$reREHE.fit[tmp$ratio==0.05]
  time_agg$`reREHE 0.1` = tmp$reREHE.fit[tmp$ratio==0.1]
  
  
  time_agg = reshape(time_agg, varying =list(names(time_agg)[-1]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[-1], timevar = 'method')
  
  # png('time_agg_method4.png',width=600, height=300, unit='px')
  # ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
  #   geom_line(size=2)
  # dev.off()
  
  time_agg$method = factor(time_agg$method, levels = c('REML', 'sREML', 'sREML est', 'REHE est', 'sREHE est' , 'sREHE CI', 'Sparsification', 'reREHE 0.05', 'reREHE 0.1' ))
  
  time_plot = ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
    geom_line(size=2)
  
  
  time_plot = ggplot(time_agg[! time_agg$method %in% c('sREML', 'sREHE est'),],aes(x=n, y=log10(time), group=method, color=method))+
    geom_line(size=1.7, aes(linetype = method))+
    geom_point(size=2)+
    scale_linetype_manual(values=c(1, 3, 1, 3, 4, 5, 3))+
    scale_color_manual(values=cbPalette[c(2, 2, 3, 3, 6, 7, 7)])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
  
  }
 
  
  ## time for dense only
  {
  tmp = data.frame(dense_table_time)
  tmp$setting = rep(1:5, 4)
  # time under same sample size are very similar, so might just average over same sample size and show time for different methods
  time_agg = aggregate(tmp,by=list(tmp$n),function(x)mean(x, na.rm=T))[-c(1, 3, 4, 5, 7, 9, 11,12)]
  colnames(time_agg)<-c('n', 'REHE CI', 'REHE est',  'REML')
  
  tmp = aggregate(reREHE_table_time, by=list(reREHE_table_time$n, reREHE_table_time$ratio), function(x)mean(x, na.rm=T))
  
  time_agg$`reREHE 0.05` = tmp$reREHE.fit[tmp$ratio==0.05]
  time_agg$`reREHE 0.1` = tmp$reREHE.fit[tmp$ratio==0.1]
  
  
  time_agg = reshape(time_agg, varying =list(names(time_agg)[2:ncol(time_agg)]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[2:ncol(time_agg)], timevar = 'method')
  time_agg$method = factor(time_agg$method, levels = c('REML', 'REHE est', 'REHE CI', 'reREHE 0.05', 'reREHE 0.1'))
 
  dense_time_plot = ggplot(time_agg,aes(x=n, y=log10(time), group=method, color=method))+
    geom_line(size=1.7, aes(linetype = method))+
    geom_point(size=2)+
    scale_linetype_manual(values=c(1,  1, 3, 5, 3))+
    scale_color_manual(values=cbPalette[c(2, 3, 3, 7, 7 )])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
  }
  
  
  }
  
   time_plot
   dense_time_plot
  ##################
  #MSE (present with root MSE (RMSE))
  load('new_reREHE_wrap_up_summary.RData')
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
  
  tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
  
  
  colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
  
  condition34 = cbind(condition34,tmp)
  
  condition34 = reshape(condition34, varying =list(names(condition34)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                        v.names='MSE', times=names(condition34)[-(1:3)], timevar = 'method')
  
  
 
  # png(paste0(start,'point_condition4.png'),width=1250, height=800, unit='px')
  # ggplot(condition34, aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
  #   geom_line(size=1.5)+
  #   scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
  #   facet_wrap(.~condition+setting, nrow=2, scales='free',
  #              labeller = labeller(condition = label_condition, setting = label_setting))
  # dev.off()
  # 
  ggplot(condition34[!condition34$method%in%c('HE-sparse', 'REHE-sparse', 'sREML'),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    facet_wrap(.~condition+setting, nrow=2, scales='free',
               labeller = labeller(condition = label_condition, setting = label_setting))
    
  condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
  
  

  MSE_plot_h = 
  ggplot(condition34[condition34$setting==1& !(condition34$method %in% c('HE-sparse', 'REHE-sparse', 'reREHE 0.05 median','reREHE 0.1 median')),], 
         aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    ylab('RMSE')+ggtitle(label ='Heritability')+
    labs(caption = bquote(~ sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
    theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
    geom_point(size=2)+
    scale_linetype_manual(values=c(1, 3, 1, 4, 5, 3))+
    scale_color_manual(values=cbPalette[c(2, 2, 3, 4, 7, 7)])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    scale_y_continuous(labels=scientific)
  
  dense_MSE_plot_h = 
    ggplot(condition34[condition34$setting==1& !(condition34$method %in% c('HE-sparse', 'REHE-sparse', 'sREML', 'reREHE 0.05 median','reREHE 0.1 median')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    ylab('RMSE')+ggtitle(label ='Heritability')+
    labs(caption = bquote(~ sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1'))+
    theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
    geom_point(size=2)+
    scale_linetype_manual(values=c(1, 1, 4, 5, 3, 5, 3))+
    scale_color_manual(values=cbPalette[c(2,  3, 4, 7, 7, 8, 8)])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    scale_y_continuous(labels=scientific)
  
  start=10
  
  condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
  colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
  tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
              reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
  
  
  colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
  
  condition34 = cbind(condition34,tmp)
  
  condition34 = reshape(condition34, varying =list(names(condition34)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                        v.names='MSE', times=names(condition34)[-(1:3)], timevar = 'method')
  
  condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
  
  
  
  MSE_plot_v = 
  ggplot(condition34[condition34$setting==4& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    ylab('RMSE')+
    ggtitle(label =bquote(sigma[1]^2 ))+ 
    labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
    theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
    geom_point(size=2)+
    scale_linetype_manual(values=c(1, 3, 1, 4, 5, 3, 5, 3))+
    scale_color_manual(values=cbPalette[c(2, 2, 3, 4, 7, 7, 8, 8)])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    scale_y_continuous(labels=scientific)
  
  
  dense_MSE_plot_v = 
    ggplot(condition34[condition34$setting==4& !(condition34$method %in% c('HE-sparse', 'REHE-sparse', 'sREML', 'reREHE 0.05 median','reREHE 0.1 median')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    ylab('RMSE')+
    ggtitle(label =bquote(sigma[1]^2 ))+ 
    labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1'))+
    theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
    geom_point(size=2)+
    scale_linetype_manual(values=c(1,  1, 4, 5, 3, 5, 3))+
    scale_color_manual(values=cbPalette[c(2,  3, 4, 7, 7, 8, 8)])+
    scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
    scale_y_continuous(labels=scientific)
  
  
  
  # more for supplement
  {
    start=7
    
    condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
    colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
    tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
    
    
    colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
    
    condition34 = cbind(condition34,tmp)
    
    condition34 = reshape(condition34, varying =list(names(condition34)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='MSE', times=names(condition34)[-(1:3)], timevar = 'method')
    
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
    
    MSE_main_sup0 = ggplot(condition34[condition34$setting==2& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('RMSE')+
      ggtitle(label =bquote(sigma[0]^2 ))+ 
      labs(caption = bquote( sigma[0]^2  == '0.04'~',' ~ sigma[1]^2 == '0.1'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 4, 5, 3, 5, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 4, 7, 7, 8, 8)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    MSE_main_sup1 = ggplot(condition34[condition34$setting==5& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('RMSE')+
      ggtitle(label =bquote(sigma[0]^2 ))+ 
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 4, 5, 3, 5, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 4, 7, 7, 8, 8)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    
    start=10
    
    condition34 = table_point_est[table_point_est$condition%in%c(3,4),c(1:4, 59, seq(start,58, by=9))][, -c(2, 3)]
    colnames(condition34)<-c('n', 'condition', 'setting', 'REHE','REHE-sparse', 'HE', 'HE-sparse', 'REML', 'sREML')
    tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition==4&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
    
    
    colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
    
    condition34 = cbind(condition34,tmp)
    
    condition34 = reshape(condition34, varying =list(names(condition34)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                          v.names='MSE', times=names(condition34)[-(1:3)], timevar = 'method')
    
    condition34$method = factor(condition34$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
    
    MSE_main_sup2 = ggplot(condition34[condition34$setting==5& !(condition34$method %in% c('HE-sparse', 'REHE-sparse')),], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      ylab('RMSE')+
      ggtitle(label =bquote(sigma[1]^2 ))+ 
      labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.01'))+
      theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5) )+    
      geom_point(size=2)+
      scale_linetype_manual(values=c(1, 3, 1, 4, 5, 3, 5, 3))+
      scale_color_manual(values=cbPalette[c(2, 2, 3, 4, 7, 7, 8, 8)])+
      scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+
      scale_y_continuous(labels=scientific)
    
    }
  
  
  
  }
  #MSE_plot_h
  #MSE_plot_v
  #dense_MSE_plot_h
  #dense_MSE_plot_v
  
  MSE_main_sup = ggarrange(    
    MSE_main_sup0+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),

    MSE_main_sup1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                           MSE_main_sup2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                           MSE_plot_v+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)), 
                           
                           labels = c('a', 'b', 'c', 'd'),nrow=2, ncol=2, common.legend = T, legend='bottom')
  
  ggsave(filename='FigS3_MSE_main_sup.pdf', plot=MSE_main_sup, device='pdf', width=7.2, height=7.2)          
  
  ##### comparison of confidence interval coverage
  
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  # for sparsification
  {
 
  table_CI = data.frame(table_CI[con4,])
  table_CI$setting = rep(1:5, 4)
  CI_NA = table_CI[, c(17,21, 1:4)]
  
  cbind(table_CI[(CI_NA[,1]>0),c(1,2,3,4,20)], CI_NA[CI_NA[,1]>0,])
  
  
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

  #png(paste0(start,'CI_condition4.png'),width=1250, height=800, unit='px')
  ggplot(condition34, aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
    geom_line(size=1.5)+
    geom_hline(yintercept = 0.95,linetype=10, color=1)+
    facet_wrap(.~condition+setting, nrow=2, scales='free',
               labeller = labeller(condition = label_condition, setting = label_setting))
  #dev.off()
  

  
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

  aligned_1 = align_plots(CIcov_plot1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5), legend.position = 'none')+scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                          CIcov_plot2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5), legend.position = 'none')+ scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                         align="hv", axis="tblr")
  
  CIcov_main_sup = ggarrange(ggarrange(plotlist = aligned_1,nrow=1, labels=c('a', 'b')), 
                             get_legend(CIcov_plot1+theme(legend.title=element_blank(), legend.position = 'bottom')), nrow=2, heights =c(5, 1))
  

  ggsave(plot=CIcov_main_sup, filename='FigS4_CIcov_main_sup.pdf', device='pdf', width=7*7.2/7, height=4*7.2/7)
  
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
    
    #png(paste0(start,'CI_condition4.png'),width=1250, height=800, unit='px')
    ggplot(condition34, aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      facet_wrap(.~condition+setting, nrow=2, scales='free',
                 labeller = labeller(condition = label_condition, setting = label_setting))
    #dev.off()
    
    ggplot(condition34[condition34$method!='sREML',], aes(x=n, y=Coverage, group=method, color=method, linetype=method))+
      geom_line(size=1.5)+
      geom_hline(yintercept = 0.95,linetype=10, color=1)+
      facet_wrap(.~condition+setting, nrow=2, scales='free',
                 labeller = labeller(condition = label_condition, setting = label_setting))
    
    
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
  
  dense_CIcov_plot1
  dense_CIcov_plot2
  
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
    
    
    ggplot(tmp, aes(x=n, y=width, fill=method))+
      geom_boxplot()+facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    
    
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
  
  ggsave(plot=CIw_main_sup, filename='FigS5_CIw_main_sup.pdf', device='pdf', width=7.2, height=7.2)
  
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
    
    ggplot(tmp[tmp$method!='sREML',], aes(x=n, y=width, fill=method))+
      geom_boxplot()+facet_wrap(.~condition+setting, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    
    
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
  dense_CIw_plot1_m
  dense_CIw_plot2_m
 
  
}
  



#############################
## results of condition 1 2 3, used in appendix
################################
{
  ## time for sparsification
  load('new_reREHE_wrap_up_summary.RData')
    load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
    load('new_wrap_up_summary.RData')
    {
      tmp = data.frame(table_time)
      tmp$setting = rep(1:5, 4)
      # time under same sample size are very similar, so might just average over same sample size and show time for different methods
      time_agg = aggregate(tmp,by=list(tmp$n),function(x)min(x, na.rm=T))[-c(1, 3, 4, 5,7, 8, 12)]
      colnames(time_agg)<-c('n', 'sREHE CI', 'sREML',  'REML', 'sparsification')
      
      tmp = aggregate(reREHE_table_time, by=list(reREHE_table_time$n, reREHE_table_time$ratio), function(x)mean(x, na.rm=T))
      time_agg$`reREHE 0.05` = tmp$reREHE.fit[tmp$ratio==0.05]
      time_agg$`reREHE 0.1` = tmp$reREHE.fit[tmp$ratio==0.1]
      
      
      
      time_agg = reshape(time_agg, varying =list(names(time_agg)[2:ncol(time_agg)]), direction = 'long', idvar='n', v.names='time', times=names(time_agg)[2:ncol(time_agg)], timevar = 'method')
      time_agg$method = factor(time_agg$method, levels = c('REML', 'sREML', 'sREHE CI', 'sparsification', 'reREHE 0.05', 'reREHE 0.1'))
      
      sup_time_plot = ggplot(time_agg[time_agg$method%in%c('REML', 'sREML', 'sREHE CI', 'sparsification'),],aes(x=n, y=log10(time), group=method, color=method))+
        geom_line(size=1.7, aes(linetype = method))+
        geom_point(size=2)+
        scale_linetype_manual(values=c(1, 3,  1, 3, 5, 5))+
        scale_color_manual(values=cbPalette[c(2,2, 3, 6 , 7, 8)])+
        scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
    }
    ggsave(filename='FigS2_sup_time_plot.pdf', plot=sup_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title = element_blank()), 
           device='pdf', width=3.5*1.8, height = 2.5*1.8)
  
  
  ##################MSE, one plot for sparse and nonsparse
  load('new_reREHE_wrap_up_summary.RData')
    
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
        tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
        
        
        colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
        
        condition123 = cbind(condition123,tmp)
        
        condition123 = reshape(condition123, varying =list(names(condition123)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                              v.names='MSE', times=names(condition123)[-(1:3)], timevar = 'method')
        
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        

        
        # png(paste0(start,'point_condition123.png'),width=1250, height=800, unit='px')
        # ggplot(condition123, aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
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
        tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
        
        
        colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
        
        condition123 = cbind(condition123,tmp)
        
        condition123 = reshape(condition123, varying =list(names(condition123)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[-(1:3)], timevar = 'method')
        
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {
        MSE11 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote( sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1' ))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        
        MSE12 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote(sigma[0]^2))+
          labs(caption = bquote( sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2' ))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE13 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+
          ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
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
        tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
        
        
        colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
        
        condition123 = cbind(condition123,tmp)
        
        condition123 = reshape(condition123, varying =list(names(condition123)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[-(1:3)], timevar = 'method')
        
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {
        MSE21 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
            scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote(sigma[1]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE22 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote(sigma[1]^2))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE23 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
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
        tmp = cbind(reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] ,
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==1 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.05, start+2, drop=F] , 
                    reREHE_table_point_est[reREHE_table_point_est$mean1median2==2 & reREHE_table_point_est$condition%in% c(1:3)&reREHE_table_point_est$ratio==0.1, start+2, drop=F] )
        
        
        colnames(tmp) = c('reREHE 0.05','reREHE 0.1','reREHE 0.05 median','reREHE 0.1 median')
        
        condition123 = cbind(condition123,tmp)
        
        condition123 = reshape(condition123, varying =list(names(condition123)[-(1:3)]), direction = 'long', idvar=c('n', 'condition','setting'), 
                               v.names='MSE', times=names(condition123)[-(1:3)], timevar = 'method')
        
        condition123$method = factor(condition123$method, levels = c('REML', 'sREML', 'REHE', 'HE', 'REHE-sparse', 'HE-sparse', colnames(tmp)))
        condition123 = condition123[!(condition123$method%in%c('HE-sparse', 'REHE-sparse')),]
        
        {        
          MSE31 = ggplot(condition123[condition123$condition==1 & condition123$setting==1,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
            scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
            scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.1'~',' ~ sigma[1]^2 == '0.1, setting 1'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
            scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        
        MSE32 = ggplot(condition123[condition123$condition==2 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 2'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
        MSE33 = ggplot(condition123[condition123$condition==3 & condition123$setting==4,], aes(x=n, y=sqrt(MSE), group=method, color=method, linetype = method))+
          geom_line(size=1.5)+ylab('RMSE')+
          geom_point(size=2)+
          scale_color_manual(values = cbPalette[c(2, 5, 3, 4, 7, 7, 8, 8)])+
          scale_linetype_manual(values=c(1, 3, 2, 3, 5, 3, 5, 3))+
          ggtitle(label =bquote('Heritability'))+
          labs(caption = bquote(sigma[0]^2  == '0.01'~',' ~ sigma[1]^2 == '0.1, setting 3'))+
          theme(plot.title = element_text(hjust=0.5), plot.subtitle =element_text(hjust=0.5),
                legend.title = element_blank()) +
          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))
        
}       
        
        # ggarrange(MSE31, MSE32, MSE33, common.legend = T, nrow=1, legend='bottom', labels = c('G', 'H', 'I'))
  }

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
            get_legend(MSE11+theme(legend.title=element_blank(), legend.position = 'bottom', legend.background = element_rect(fill = "transparent", colour = NA))),
            nrow=2, heights=c(18, 1))
  ggsave(filename ='FigS6_suppfigure1.pdf', plot=suppfig1,  device='pdf',  width=8.5, height=9 )



  ##### comparison of confidence interval coverage, dense and sparse together
  load('new_dense_wrap_up_summary.RData') # mark with dense_ for three results tables and collections
  load('new_wrap_up_summary.RData')
  
  {

  table_CI = data.frame(table_CI)
  dense_table_CI = data.frame(dense_table_CI)
  table_CI$setting = rep(rep(1:5, each=4), 4)
  CI_NA = table_CI[, c(17,21, 1:5)]
  
  cbind(table_CI[(CI_NA[,1]>0),c(1,2,3,4,20)], CI_NA[CI_NA[,1]>0,])
  
  
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
    get_legend(CI11+theme(legend.title=element_blank(), legend.position = 'bottom', legend.background = element_rect(fill = "transparent", colour = NA))),
    nrow=2, heights=c(18, 1))
  ggsave(filename ='FigS7_suppfigure2.pdf', plot=suppfig2,  device='pdf',  width=8.5, height=9)
  
        

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
    png(paste0('e2CI_width.png'),width=1500, height=1250, unit='px')
      ggplot(tmp, aes(x=n, y=width, fill=method))+
        geom_boxplot()+
        facet_wrap(.~condition+setting,ncol=5, scales = 'free', labeller = labeller(condition = label_condition, setting = label_setting))
    dev.off()

  
    
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
  get_legend(CIw11+theme(legend.title=element_blank(), legend.position = 'bottom', legend.background = element_rect(fill = "transparent", colour = NA))),
  nrow=2, heights=c(18, 1))
  ggsave(filename ='FigS8_suppfigure3.pdf', plot=suppfig3,  device='pdf',  width=8.5, height=9 )
  
   
  }
  
  


##### get the HE estimation zero proportion


zero_prop = cbind(table_time[,1:4],t(sapply(1:80, function(i)rowMeans(he_est[[i]][['e_k']]==0))))

zero_prop[zero_prop[,4]==4,]
zero_prop[zero_prop[,5]>0 | zero_prop[,6]>0,]



zero_prop_sub = cbind(reREHE_table_point_est[,1:6],
                  rbind(
                    t(sapply(1:160, function(i)rowMeans(rehe_sub_est[[i]][['e_k_median']]==0))), 
                  t(sapply(1:160, function(i)rowMeans(rehe_sub_est[[i]][['e_k_mean']]==0)))))

zero_prop_sub[zero_prop_sub[,7]>0 | zero_prop_sub[,8]>0,]


cbind(zero_prop, zero_prop_sub[1:160,5:8])[81:160,]


##### get the real data application presentation (files are too large to share)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
source('../R/GENESIS/varCompCI.R')

load("../nullModelREHE_LABA2.RData")
timeNMM_REHE_noCI = timeNMM_REHE

load("../nullModelREHE_LABA2_with_CI.RData")
load("../nullModel_LABA2.RData")

name = c('01_mean', '01_median', '005_mean','005_median')
timeNMM_reREHE_list=list()
nullMMreREHE_list=list()
for (i in 1:4){
  load(paste0('../nullModelreREHE',name[i], '_LABA2.RData'))
  timeNMM_reREHE_list[[i]] = timeNMM_reREHE
  nullMMreREHE_list[[i]] = nullMMreREHE
}

# look at time
timeNMM_REHE_noCI
timeNMM_REHE
timeNMM_REML
timeNMM_reREHE_list

# for point estimation
{
  nullMM$varComp
  nullMMREHE$varComp

  name
  nullMMreREHE_list[[1]]$varComp
  nullMMreREHE_list[[2]]$varComp
  nullMMreREHE_list[[3]]$varComp
  nullMMreREHE_list[[4]]$varComp
  
  # library(kableExtra)
  # kable(data.frame(CIvc[, 1:2], 'REML CI'= paste('[', paste(CIvc[,3], CIvc[,4], sep=','),']'),
  #                  'REHE Wald' =  paste('[', paste(CIvc[,5], CIvc[,6], sep=','),']'),
  #                  'REHE Quantile' =  paste('[', paste(CIvc[,7], CIvc[,8], sep=','),']')),
  #       booktab=T, 'latex')%>%
  #   kable_styling()%>%
  #   add_header_above(c(' '=1, 'Point Est'=2, '95 CI' = 3))
  
}

# for confidence interval
{
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



# library(kableExtra)
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

}


## for association testing (files are too large to share)

{
load("../chr1REHE.RData")
load("../chr1.RData")

assocreREHE_list = list()
assocreREHE_part = list()
assoc_part = assoc[assoc$Score.pval<=5e-6,]
assoc_part = assoc_part[assoc_part$Score.pval<=5e-8,] 

for (i in 1:4){
  load(paste0('../chr1reREHE',name[i], '.RData'))
  assocreREHE_list[[i]] = assocreREHE
  assocreREHE_part[[i]] = assocreREHE[assocreREHE$Score.pval<=5e-8,]
  all((assoc_part$variant.id %in% assocreREHE_part$variant.id)) # if reREHE contains all genes picked by REML
}

 

assocREHE_part = assocREHE[assocREHE$Score.pval<=5e-8,]

all((assoc_part$variant.id %in% assocREHE_part$variant.id)) # REHE contains all genes picked by REML

include_index = unique(c(assocREHE_part$variant.id,assoc_part$variant.id, 
                         do.call(c, sapply(1:2, function(w) assocreREHE_part[[w]]$variant.id))))

assoc_part = assoc[assoc$variant.id %in% include_index,]
assocREHE_part =  assocREHE[assocREHE$variant.id %in% include_index,]

for (i in 1:4){
  assocreREHE_part[[i]] = assocreREHE_list[[i]][assocreREHE_list[[i]]$variant.id %in% include_index, ]
}

all(assoc_part$variant.id == assocREHE_part$variant.id)
all(assoc_part$variant.id == assocreREHE_part[[2]]$variant.id)


datapvalreREHE1 = ggplot(data=data.frame(), 
                         aes(x=-log10(assocreREHE_part[[1]]$Score.pval), 
                             y=-log10(assoc_part$Score.pval)))+
  geom_point()+
  xlab('reREHE p values (-log10)')+
  ylab('REML p values (-log10)')+
  geom_hline(aes(yintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_vline(aes(xintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_abline(aes(intercept=0, slope=1), color=2)

datapvalreREHE2 = ggplot(data=data.frame(), aes(x=-log10(assocreREHE_part[[2]]$Score.pval),
                                                y=-log10(assoc_part$Score.pval) ))+
  geom_point()+
  xlab('reREHE p values (-log10)')+
  ylab('REML p values (-log10)')+
  geom_hline(aes(yintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_vline(aes(xintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_abline(aes(intercept=0, slope=1), color=2)


datapval = ggplot(data=data.frame(), aes(x=-log10(assocREHE_part$Score.pval), 
                                         
                                         y=-log10(assoc_part$Score.pval)))+
  geom_point()+
  xlab('REHE p values (-log10)')+
  ylab('REML p values (-log10)')+
  geom_hline(aes(yintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_vline(aes(xintercept = -log10(5e-8)), color=2, linetype=2)+
  geom_abline(aes(intercept=0, slope=1), color=2)
  

datapvalREHEreREHE = ggplot(data=data.frame(), aes(x=assocreREHE_part[[1]]$Score.pval, y=assocREHE_part$Score.pval))+
  geom_point()+
  xlab('reREHE p values')+
  ylab('REHE p values')+
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

save(file='real_data_plots.RData',  list=c('datapval','datapvalreREHE1', 'datapvalreREHE2', 'dataCI', 'datapvalREHEreREHE'))


## data objects are saved here
load('real_data_plots.RData')

# main1 = ggarrange(
#   
#   ggarrange(dense_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'bottom', legend.title=element_blank()), datapval, heights = c(4, 4), nrow=1, labels=c('a', 'b'), align='hv', axis='tblr'),
#           
#   ggarrange(dense_MSE_plot_h+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.title = element_blank()), dense_MSE_plot_v+scale_y_continuous(labels = scientific), nrow=1, common.legend=T, legend='bottom', labels = c('c', 'd'), heights=c(0.8, 1)),
#           nrow=2,heights=c(0.8, 1), align = "hv")

aligned_1 = align_plots(dense_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'bottom', legend.title=element_blank())
                        , datapval+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5)),
                        dense_MSE_plot_h+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), 
                        dense_MSE_plot_v+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), align="hv", axis="tblr")

main1 = ggarrange(ggarrange(plotlist = aligned_1, ncol=2, nrow=2, labels=c('a', 'b', 'c', 'd')), 
                  get_legend(dense_MSE_plot_v+theme(legend.position = 'bottom', legend.title = element_blank(), legend.background = element_rect(fill = "transparent", colour = NA))), nrow=2, heights=c(8, 1))

# ggsave(filename = 'main1.pdf', plot=main1, device='pdf', width=7.2, height=7.2)



aligned_1 = align_plots(ggplot() + theme_void(), dense_time_plot+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'bottom', legend.title=element_blank())
                        , ggplot() + theme_void(),
                        dense_MSE_plot_h+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), 
                        dense_MSE_plot_v+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'), align="hv", axis="tblr")



main1 = ggarrange(ggarrange(plotlist=aligned_1[1:3], ncol=3, nrow=1, labels=c('', 'a', ''), widths=c(0.25, 0.5, 0.25)),
                  ggarrange(plotlist = aligned_1[4:5], ncol=2, nrow=1, labels=c('b', 'c')), 
                  get_legend(dense_MSE_plot_v+theme(legend.position = 'bottom', legend.title = element_blank())), 
                  nrow=3, heights=c(8, 8, 1))

ggsave(filename = 'Fig3_biost_main1.pdf', plot=main1, device='pdf', width=7.2, height=7.2)


name[1]
aligned_1.1 = align_plots(datapval+ theme_bw()+theme(axis.text=element_text(size=7)), 
                           datapvalreREHE1+ theme_bw()+theme(axis.text=element_text(size=7)),
                           #datapvalREHEreREHE+xlim(c(NA, 8e-8))+ theme_bw()+theme(axis.text=element_text(size=7)), 
                          
                          align="hv", axis="tblr")

main0 = ggarrange(ggarrange(plotlist=aligned_1.1, ncol=2, labels=c('a', 'b')),
dataCI, nrow=2, labels=c('', 'c'), heights = c(1, 1))
ggsave(filename = 'Fig1_biost_main0.pdf', plot=main0, device='pdf', width=7.2, height=7.2)



main2= ggarrange(
  
  ggarrange(dense_CIcov_plot1+
              scale_x_continuous(breaks=c(3000, 6000, 9000, 12000))+theme(legend.title = element_blank()), 
            dense_CIcov_plot2+
              scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)), 
          
  dense_CIw_plot1_m+theme(legend.title = element_blank()), dense_CIw_plot2_m, 
  
  nrow=2,ncol=2,  common.legend=T, legend='bottom', labels = c('a', 'b','c', 'd'))
  
          
)

aligned_2 = align_plots(dense_CIcov_plot1+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none')+
                          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)), 
                        dense_CIcov_plot2+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none')+
                          scale_x_continuous(breaks=c(3000, 6000, 9000, 12000)),
                        dense_CIw_plot1_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'),
                        dense_CIw_plot2_m+ theme_bw()+theme(legend.title=element_blank(),plot.title = element_text(hjust = 0.5))+theme(legend.position = 'none'),
                         align="hv", axis="tblr")

main2 = ggarrange(ggarrange(plotlist = aligned_2, ncol=2, nrow=2, labels=c('a', 'b', 'c', 'd')), 
                  get_legend(dense_CIw_plot2_m+theme(legend.position = 'bottom', legend.title = element_blank())), nrow=2, heights=c(8, 1))



ggsave(filename = 'Fig4_main2.pdf', plot=main2, device='pdf', width=7.2, height=7.2)












 
  
  
  
  
  