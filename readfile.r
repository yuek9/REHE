library(reticulate)

# load the virtual environment, Jing you may need to run: use_virtualenv('xx')
use_condaenv('C:\\Users\\yuek\\AppData\\Local\\r-miniconda\\envs\\r-reticulate', required=T) #BOX

# load VAE simulation results
setwd('\\\\fs2-vip-nfs.nfs.biost.priv\\students\\yuek\\Desktop/VAE')
source_python('readfile.py')

setwd('VAEres')

for (input in 1:9){
  input_list = read.table('input_arg.txt', header=T)
  
  n = as.numeric(input_list[input,1])
  p = as.numeric(input_list[input,2])
  library_scale = as.numeric(input_list[input,3])
  theta = as.numeric(input_list[input,4])
  option  = input_list[input,5]
  z_scale = as.numeric(input_list[input,6])
  print(input_list[input,])
  
  true_comp=readfile('', n, p, library_scale, theta, option, z_scale, '_true_comp.npy')
  true_count = readfile('', n, p, library_scale, theta, option, z_scale, '_true_count.npy')
  VAE_comp = readfile('', n, p, library_scale, theta, option, z_scale, '_VAE_comp.npy')
  naive_comp = readfile('', n, p, library_scale, theta, option, z_scale, '_naive_comp.npy')
  VAE_time = readfile('', n, p, library_scale, theta, option, z_scale, '_VAE_time.npy')
  
  # the loaded results will be arrays 
}



