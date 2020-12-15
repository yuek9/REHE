##
## Kun Yue
## yuek@uw.edu
## 12/14/2020

# using Mike's github code and data set to estimate the full network
# 
# devtools::install_github('mikehellstern/netgsa')
 
library(netgsa)
library(graphite)
library(data.table)
library(tidyr)
library(glassoFast)
ls()

today <- '20200827'

filepath = '...'

set.seed(100)

setwd(filepath)
load(paste0(filepath, '../data/breastcancer2012.rda'))
write.table(edgelist, file = 'edgelist.txt',row.names = F)
write.table(nonedgelist, file= 'nonedgelist.txt', row.names = F)


# note that the edge list is not symmetric/is directed
database_search <- obtainEdgeList(rownames(x), c("kegg", "reactome", "biocarta"))
save(list='database_search', file = paste0(today, 'edgeList_database.rda'))

## only use a subsample of the edge lists
# load('20200827edgeList_database.rda')
set.seed(10)
tmp = database_search$edgelist
# database_search$edgelist =  tmp[sample(nrow(tmp),size=500, replace = F),]
save(list='database_search', file = paste0(today, 'edgeList_database.rda'))



group = group[sample(1:length(group), replace=F)]

network_info <- prepareAdjMat(x, group, database_search,
                              cluster = TRUE, 
                              file_e = "edgelist.txt", 
                              file_ne = "nonedgelist.txt")

save(list=c('network_info', 'group'), file = paste0(today, 'permute_breastcancer2012_network.rda'))

# load(paste0(today, 'full_breastcancer2012_network.rda'))





