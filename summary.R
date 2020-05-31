VAE_part = read.table('\\\\fs2-vip-nfs.nfs.biost.priv\\students\\yuek\\Desktop\\VAE/VAE_res_out.txt')
colnames(VAE_part)<-c('method','n','p','distribution','theta','z_scale','time','frob','l1','shannon', 'simpson', 'KL', 'mean_sparsity')
VAE_part

library(SeqVarTools)
path = '\\\\fs2-vip-nfs.nfs.biost.priv\\students\\yuek\\Desktop\\VAE\\composition-estimate-master\\Rres'
setwd(path)
files = list.files('\\\\fs2-vip-nfs.nfs.biost.priv\\students\\yuek\\Desktop\\VAE\\composition-estimate-master\\Rres')

for (ii in 1:length(files)){
  load(files[ii])
  print(files[ii])
  tmp = rowMeans(eval)
  VAE_part= rbind(VAE_part,  data.frame('method'='Cao', 'n'=n, 'p'=p, 'distribution'=option, 'theta'=theta, 'z_scale'=z_scale, 
                        'time'=mean(do.call(c, Cao_time)), 'frob'=tmp[1], 'l1'=tmp[2],'shannon'=tmp[3], 'simpson'=tmp[4], 'KL' = tmp[5], 'mean_sparsity' = mean(sapply(1:nreps, function(w)mean(count_data[[w]]==0)))))
}

(with(VAE_part, VAE_part[n==300 & p==100 & distribution=='logisticnormal' & theta==8,])) # higher sparsity makes VAE stands out
(with(VAE_part, VAE_part[n==300 & p==100 & distribution=='logisticnormal' & theta==4,]))

(with(VAE_part, VAE_part[n==300 & p==100 & distribution=='dirichlet' & z_scale==10,]))
(with(VAE_part, VAE_part[n==300 & p==100 & distribution=='dirichlet' & z_scale==15,]))
(with(VAE_part, VAE_part[n==300 & p==100 & distribution=='dirichlet' & z_scale==25,])) # VAE always much better in terms of KL distance and Shannon index, even with lower sparsity, 
