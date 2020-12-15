setwd('...')
load('20200827_data_application_res.rda')

time

result = data.frame(pathway = out$REML$results$pathway[order(out$REML$results$pathway)], 
                    REMLpval = out$REML$results$pval[order(out$REML$results$pathway)], 
                    REHEpval = out$REHE$results$pval[order(out$REHE$results$pathway)], 
                    reREHEpval = out$reREHE$results$pval[order(out$reREHE$results$pathway)])
library(ggplot2)
library(ggpubr)
library(cowplot)
theme_set(theme_bw())

net1 = ggplot(result, aes(x=-log10(REHEpval), y=-log10(REMLpval)))+
  geom_point()+
  xlab('REHE p values (-log10)')+
  ylab('REML p values (-log10) ')+
  geom_abline(aes(intercept=0, slope=1), color=2)+
  coord_cartesian(xlim=c(0, 60), ylim=c(0, 60))
  #geom_hline(aes(yintercept = log10(5e-2)), color=2, linetype=2)+
  #geom_vline(aes(xintercept = log10(5e-2)), color=2, linetype=2)
  
net2 = ggplot(result, aes(x=-log10(reREHEpval), y=-log10(REMLpval)))+
  geom_point()+
  xlab('reREHE p values (-log10)')+
  ylab('REML p values (-log10)')+
  geom_abline(aes(intercept=0, slope=1), color=2)+
  #geom_hline(aes(yintercept = 5e-2), color=2, linetype=2)+
  #geom_vline(aes(xintercept = 5e-2), color=2, linetype=2)+
  theme_bw()+
  coord_cartesian(xlim=c(0, 60), ylim=c(0, 60))

# net3 = ggplot(result, aes(x=reREHEpval, y=REHEpval))+
#   geom_point()+
#   xlab('reREHE p values')+
#   ylab('REHE p values')+
#   geom_abline(aes(intercept=0, slope=1), color=2)+
#   geom_hline(aes(yintercept = 5e-2), color=2, linetype=2)+
#   geom_vline(aes(xintercept = 5e-2), color=2, linetype=2)+
#   theme_bw()

aligned_1 = align_plots(net1,net2,align="hv", axis="tblr")

main1 = ggarrange(plotlist = aligned_1, ncol=length(aligned_1), nrow=1, labels=c('a', 'b')) 

ggsave(filename = 'Fig2_biost_main_netgsa.pdf', plot=main1, device='pdf', width=7.2, height=3.5)

result[order(result$REMLpval)[1:2],]

out$REML$s2.epsilon
out$REML$s2.gamma
