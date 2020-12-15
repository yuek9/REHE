# to summarize the output for netgsa simulation
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

filepath = '...'
setwd(filepath)
muvals=  0.1
s2e = 0.5
s2g = 0.1
iter = 1
today = '20200827'

out_path = paste0('Output/semi_synthetic/e',s2e, 'g', s2g, 'mu',muvals)


load(file = paste0(out_path,'/', today, '_iter_', iter, '.rda'))

process = function(res){
  res[order(res$pathway),]
}

npath = nrow(out$REHE$results)

# the true power order follows the pathways_mat matrix, should make sure it match the order
load('../data/breastcancer2012.rda')
pathways_mat_name = rownames(pathways_mat)

s_vals = matrix(c(0.1, 1, 1, 0.1, 0.5, 0.1, 0.1, 0.5, 0.1, 0.1), ncol=2, byrow=T)
all_res = lapply(1:nrow(s_vals), function(s_row){
  s2e = s_vals[s_row, 1]
  s2g = s_vals[s_row, 2]
  all_res = lapply(c(0.1, 0.2, 0.3), function(muvals){
    
    out_path = paste0('Output/semi_synthetic/e',s2e, 'g', s2g, 'mu',muvals)
    
    tmp =  matrix(NA, nrow=npath, ncol = 200, dimnames = list(process(out$REHE$results)$pathway, NULL))
    
    pval_sum = list(REHE = tmp, reREHE = tmp, REML =tmp)
    time_sum = matrix(NA, nrow=200, ncol=3, dimnames = list(NULL, c('REHE',  'reREHE', 'REML')))
    vc_e <- vc_g <- matrix(NA, nrow=200, ncol=3, dimnames = list(NULL, c('REHE', 'reREHE', 'REML')))
    
    
    for (iter in 1:200){
      print(c(muvals, iter))
      tryCatch(
        {
          load(file = paste0(out_path,'/', today, '_iter_', iter, '.rda'))
          tmp = process(out$REHE$results)
          pval_sum$REHE[tmp$pathway,iter] = tmp$pFdr
          
          tmp = process(out$reREHE$results)
          pval_sum$reREHE[tmp$pathway,iter] = tmp$pFdr
          
          tmp = process(out$REML$results)
          pval_sum$REML[tmp$pathway,iter] = tmp$pFdr
          
          time_sum[iter, colnames(time)] = time
          
          c(out$REHE$s2.epsilon, out$REHE$s2.gamma)
          
          vc_e[iter, ] = c(out$REHE$s2.epsilon, out$reREHE$s2.epsilon, out$REML$s2.epsilon)
          vc_g[iter, ] = c(out$REHE$s2.gamma, out$reREHE$s2.gamma, out$REML$s2.gamma)
          
        }
        , error = function(e)e
      )
    }
    
    load(file = paste0(out_path,'/truepower.rda'))
    
    time_sum = data.frame(time = colMeans(time_sum, na.rm = T), method = colnames(time_sum))
    alpha = 0.05
    power_sum = data.frame(power = c(rowMeans(pval_sum$REHE<alpha, na.rm = T), 
                                     rowMeans(pval_sum$reREHE<alpha, na.rm = T), 
                                     rowMeans(pval_sum$REML<alpha, na.rm = T),
                                     truepower),
                           pathway = c(rep(rownames(pval_sum$REHE), 3), pathways_mat_name),
                           method = rep(c('REHE', 'reREHE', 'REML', 'truepower'), each = npath))
    
    time_sum$mu <- power_sum$mu <- muvals
    
    vc_e = cbind(vc_e, mu = muvals)
    vc_g = cbind(vc_g, mu = muvals)
    
    return(list(time = time_sum, power = power_sum, vc_e = vc_e, vc_g = vc_g))
  })
  
  time = do.call(rbind, sapply(all_res, `[`, 1))
  time = aggregate(time$time, by = list(method = time$method), FUN=mean)
  power =  do.call(rbind, sapply(all_res, `[`, 2))
  
  
  vc_e =  do.call(rbind, sapply(all_res, `[`, 3))
  vc_g =  do.call(rbind, sapply(all_res, `[`, 4))
  
  return(list(time = time, power=power, vc_e = vc_e, vc_g = vc_g))
})

power = do.call(rbind, lapply(1:nrow(s_vals), function(i){
  tmp = data.frame(all_res[[i]][[2]])
  tmp$s2e = s_vals[i,1]
  tmp$s2g = s_vals[i,2]
  tmp
  }))

time = do.call(rbind,lapply(1:nrow(s_vals), function(i){
  tmp = all_res[[i]][[1]]
  tmp$s2e = s_vals[i,1]
  tmp$s2g = s_vals[i,2]
  tmp
}))

head(power)
head(time)

library(ggplot2)
library(cowplot)
theme_set(theme_bw())

ggplot(data = power[power$mu==0.1 & power$s2e==0.1 & power$s2g == 0.1, ], aes(x=pathway, y=power, group = method, color=method))+
  geom_point()

power_wide = reshape(power,v.names = c( 'power'), idvar = c('pathway', 'mu', 's2e', 's2g'), timevar = 'method', direction='wide')
head(power_wide)

plt=ggplot(power_wide, aes(x=power.REHE, y=power.REML))+
  geom_point(size=1.1)+
  xlab('REHE empirical power')+
  ylab('REML empirical power')+
  geom_abline(slope=1, intercept = 0, color=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_bquote(cols = {mu == .(mu)}*","*{sigma[0]==.(s2e)}*", "*{sigma[1]==.(s2g)}), nrow = 3)
plt
# ggsave('reml_rehe_netgsa_sim_power.pdf', plt, width = 6*1.2, height = 4*1.2, units = 'in', device = 'pdf')

tmp1 = power_wide[,c(1,2,3,4,7,8)]
colnames(tmp1)[5]<-'power'
tmp2 = power_wide[,c(1,2,3,4,5,8)]
colnames(tmp2)[5]<-'power'
tmp3 = power_wide[,c(1,2,3,4,6,8)]
colnames(tmp3)[5]<-'power'
power_wide_tmp = rbind(data.frame(tmp1, method='REML'), 
                       data.frame(tmp2, method='REHE'),
                       data.frame(tmp3, method='reREHE'))

power1 = ggplot(power_wide_tmp, aes(x=power.truepower, y=power, group=method, color = method, shape=method))+
  geom_point(size=1.1)+
  xlab('true power')+
  ylab('empirical power')+
  theme(legend.title = element_blank())+
  scale_color_manual(values = cbPalette[c(2,3, 7)])+
  scale_shape_manual(values = c(6, 1, 3))+
  geom_abline(slope=1, intercept = 0, color=1, linetype=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_bquote(cols = {mu == .(mu)}*", "*{sigma[0]==.(s2e)}*", "*{sigma[1]==.(s2g)}), nrow = 3)

ggsave(power1, file='FigS10_netgsa_power_reml_rehe_rerehe.pdf', width = 8*1.2, height = 6*1.2, units = 'in', device = 'pdf')


tmp3 = power_wide[,c(1,2,3,4,6,8)]
colnames(tmp3)[5]<-'power'
power_wide_tmp = rbind(data.frame(tmp1, method='REML'), 
                       data.frame(tmp3, method='reREHE'))

power2 = ggplot(power_wide_tmp, aes(x=power.truepower, y=power, group=method, color = method, shape=method))+
  geom_point(size=1.1)+
  xlab('true power')+
  ylab('empirical power')+
  theme(legend.title = element_blank())+
  scale_color_manual(values = cbPalette[c(2,7)])+
  geom_abline(slope=1, intercept = 0, color=1, linetype=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_bquote(cols = {mu == .(mu)}*", "*{sigma[0]==.(s2e)}*", "*{sigma[1]==.(s2g)}), nrow = 3)

# ggsave(power2, file='plot/netgsa_power_reml_rerehe.pdf', width = 8*1.2, height = 6*1.2, units = 'in', device = 'pdf')

plt=ggplot(power_wide, aes(x=power.reREHE, y=power.REML))+
  geom_point(size=0.8)+
  xlab('reREHE empirical power')+
  ylab('REML empirical power')+
  geom_abline(slope=1, intercept = 0, color=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_both, nrow = 3)
plt
# ggsave('reml_rerehe_netgsa_sim_power.pdf', plt, width = 6*1.2, height = 4*1.2, units = 'in', device = 'pdf')


plt=ggplot(power_wide, aes(x=power.REML, y=power.truepower))+
  geom_point(size=0.8)+
  xlab('REML empirical power by pathway')+
  ylab('True power by pathway')+
  geom_abline(slope=1, intercept = 0, color=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_both, nrow = 3)
plt

plt=ggplot(power_wide, aes(x=power.REHE, y=power.truepower))+
  geom_point(size=0.8)+
  xlab('REHE empirical power by pathway')+
  ylab('True power by pathway')+
  geom_abline(slope=1, intercept = 0, color=2)+
  facet_wrap(~mu+s2e+s2g, labeller = label_both, nrow = 3)
plt


## look at variance component values
vc_e =  do.call(rbind,lapply(1:nrow(s_vals), function(i){
  tmp = data.frame(all_res[[i]][[3]])
  tmp$s2e = s_vals[i,1]
  tmp$s2g = s_vals[i,2]
  tmp
}))
vc_g =  do.call(rbind,lapply(1:nrow(s_vals), function(i){
  tmp = data.frame(all_res[[i]][[4]])
  tmp$s2e = s_vals[i,1]
  tmp$s2g = s_vals[i,2]
  tmp
}))
head(vc_e)

vc_e_long = reshape(data.frame(vc_e), direction='long', varying = c('REML', 'REHE', 'reREHE'),  v.names = 'value', times =c('REML', 'REHE', 'reREHE') )

vc_e_long$time = factor(vc_e_long$time, levels = c('REML', 'REHE', 'reREHE' ))

ggplot(data = vc_e_long, aes(x=time, y=value, group=time, color=time))+
  geom_boxplot()+
  ylab('noise variance component')+
  xlab('method')+
  scale_color_manual(name=NULL,values=c(1, 2,3 ))+
  facet_wrap(~mu+s2e+s2g, labeller = label_both, nrow=3)





vc_g_long = reshape(data.frame(vc_g), direction='long', varying = c('REML', 'REHE', 'reREHE'),  v.names = 'value', times =c('REML', 'REHE', 'reREHE') )

vc_g_long$time = factor(vc_g_long$time, levels = c('REML', 'REHE', 'reREHE' ))

ggplot(data = vc_g_long, aes(x=time, y=value, group=time, color=time))+
  geom_boxplot()+
  ylab('network information variance component')+
  xlab('method')+
  scale_color_manual(name=NULL,values=c(1, 2,3 ))+
  facet_wrap(~mu+s2e+s2g, labeller = label_both, nrow=3)


# the results not affected by mu, so just look at mu=0.1 results
plte = ggplot(data = vc_e_long, aes(x=time, y=value-s2e, group=time, color=time))+
  geom_boxplot()+
  ylab(bquote(hat(sigma[0]^2)-sigma[0]^2))+
  xlab(NULL)+
  theme(axis.text.x = element_blank())+
  scale_color_manual(name=NULL,values=cbPalette[c(2, 3, 7)])+
  facet_wrap(~s2e+s2g, labeller = label_bquote(cols = {sigma[0]==.(s2e)}*", "*{sigma[1]==.(s2g)}), nrow=1)+
  geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=c(-0.30, -0.15, -0, 0.15, 0.30, 0.45, 0.60, 0.75))

pltg = ggplot(data = vc_g_long, aes(x=time, y=value-s2g, group=time, color=time))+
  geom_boxplot()+
  ylab(bquote(hat(sigma[1]^2)-sigma[1]^2))+
  xlab(NULL)+
  theme(axis.text.x = element_blank())+
  scale_color_manual(name=NULL,values=cbPalette[c(2, 3, 7)])+
  facet_wrap(~s2e+s2g, labeller = label_bquote(cols = {sigma[0]==.(s2e)}*", "*{sigma[1]==.(s2g)}), nrow=1)+
  geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=c(-0.75, -0.60, -0.45, -0.30, -0.15, -0, 0.15, 0.30))



library(ggpubr)
out1 = ggarrange(plte, pltg, nrow=2, common.legend = T, legend = 'right')
ggsave(out1, file='FigS9_netgsa_e_g_estimation.pdf', width = 6*1.2, height = 4*1.2, units = 'in', device = 'pdf')


aggregate(vc_e_long$value, FUN=mean, by = list(vc_e_long$mu, vc_e_long$time))
aggregate(vc_g_long$value, FUN=mean, by = list(vc_g_long$mu, vc_g_long$time))
