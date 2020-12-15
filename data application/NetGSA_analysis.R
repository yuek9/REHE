# data application: analyze the breast cancer data set provided by the netgsa package (from mike's github)
library(netgsa)
library(graphite)
library(data.table)
library(tidyr)
library(glassoFast)
ls()

today <- '20200827'

filepath = '..'

setwd(filepath)


load('breastcancer2012.rda')

set.seed(2020)

n = as.vector(table(group))
p = nrow(x)

x = list(x[,group==1], x[,group==2])
group = group[order(group)]

x = do.call(cbind, x)

load(paste0(today, 'edgeList_database.rda'))


network_info <- prepareAdjMat(x, group, database_search,
                              cluster = TRUE, 
                              file_e = "edgelist.txt", 
                              file_ne = "nonedgelist.txt")


time <- matrix(0, nrow=1, ncol=3, dimnames = list(NULL, c('REHE', 'reREHE', 'REML')))


cat('run REHE \n')

time[1,'REHE'] = sum(system.time({
  out_REHE <- NetGSA(network_info[["Adj"]], 
                     x, group, pathways_mat,
                     lklMethod = "REHE", sampling = F)})[1:2])

cat('run reREHE \n')

time[1,'reREHE'] = sum(system.time({
  out_reREHE <- NetGSA(network_info[["Adj"]], 
                       x, group, pathways_mat,
                       lklMethod = "REHE", sampling = T, 
                       sample_n = 0.1, sample_p=0.5)})[1:2])


cat('run REML \n')

time[1,'REML'] = sum(system.time({
  out_REML <- NetGSA(network_info[["Adj"]], 
                     x, group, pathways_mat,
                     lklMethod = "REML")})[1:2])


out = list(REHE = out_REHE, reREHE = out_reREHE, REML = out_REML)
save(list= c('time', 'out'), file = paste0(today, '_data_application_res.rda'))


# some inspection of the data features:

dim(x)
table(group)
dim(pathways_mat)