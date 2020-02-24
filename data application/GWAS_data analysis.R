## Kun Yue, yuek@uw.edu
## 02/24/2020

# file for data application 

filepath = '...'
setwd(paste0(filepath, '/HCHS_SOL'))

library(gdsfmt)
library(GWASTools) 
library(GENESIS)
sessionInfo()


if(T){

# Load scan annot
# obtain the outcome variable: LABA2
scanAnnot<- getobj(paste0(filepath, "/scanAnnot_316987.RData"))

outcome <- "LABA2"
covars <- c("factor.sex",
            "factor.CENTER",
            "AGE",
            "factor.CIGARETTE_USE",
            "factor.gengrp6",
            "log.WEIGHT_FINAL_NORM_OVERALL",
            "EV1", "EV2","EV3","EV4","EV5")


head(pData(scanAnnot))
scan.include <- scanAnnot$scanID[!is.na(scanAnnot$LABA2)]
length(scan.include) # [1] 12704, no missing data for LABA2

# Load kinship matrix
phi <- getobj(paste0(filepath, "/pcreap_phi_matrix_without_asian_outliers.RData"))
dim(phi) # [1] 12784 12784
phi[1:5,1:5]
phi.ids <- colnames(phi)
all(scan.include%in% phi.ids) # FALSE; some subjects not available in kinship matrix, remove those
scan.include = scan.include[(scan.include%in% phi.ids)]
length(scan.include) #12686
phi = phi[as.character(scan.include), as.character(scan.include)]

# Load census_block matrix
census <- getobj(paste0(filepath, "/covar_matrix_PSU_ID.RData"))
dim(census) 
all(scan.include %in% colnames(census)) #True
census = census[as.character(scan.include), as.character(scan.include)]

# Load household matrix
house <- getobj(paste0(filepath,"/covar_matrix_HH_ID.RData"))
dim(house) # [1] 12803 12803
all(scan.include %in% colnames(house)) #True
house = house[as.character(scan.include), as.character(scan.include)]


n=length(scan.include)
scan.include_tmp = scan.include[1:n]
scanAnnot_tmp = scanAnnot
pData(scanAnnot_tmp)<-pData(scanAnnot)[scanAnnot$scanID%in%scan.include_tmp,]
phi_tmp = phi[as.character(scan.include_tmp), as.character(scan.include_tmp)]
census_tmp = census[as.character(scan.include_tmp), as.character(scan.include_tmp)]
house_tmp = house[as.character(scan.include_tmp), as.character(scan.include_tmp)]





# Perform Association Testing

# Run Null Model 
if(T){
  timeNMM_REML = system.time({
    nullMM <- fitNullModel(scanAnnot,      			
                           outcome=outcome, 
                           covars = covars,
                           cov.mat=list(phi=phi_tmp,census=census_tmp,house=house_tmp),
                           sample.id = scan.include_tmp,						 
                           family = gaussian)
    
  })
  save(list=c('timeNMM_REML','nullMM', 'scanAnnot_tmp', 'scan.include_tmp'), file="nullModel_LABA2.RData")

}


cat('NullModel fit completed by GENESIS \n')


library(quadprog)

timeNMM_REHE = system.time({
nullMMREHE <- fitNullModelREHE(scanAnnot,      			
                       outcome=outcome, 
                       covars = covars,
                       cov.mat=list(phi=phi_tmp,census=census_tmp,house=house_tmp),
                       sample.id = scan.include_tmp,						 
                       family = gaussian)
})
cat('fit REHE success, now saving file \n')
save(list=c('timeNMM_REHE','nullMMREHE','scanAnnot_tmp', 'scan.include_tmp'), file="nullModelREHE_LABA2.RData")

cat('NullModel fit completed by REHE \n')

rm(list=c('phi', 'house', 'census'))

save.image(file = 'assotestchr1res.RData')

}

#########################################################

load('nullModelREHE_LABA2.RData')
load('nullModel_LABA2.RData')

cat('work on gds data \n')
# load genotype data
gdata <- openfn.gds(paste0(filepath, "/SOL_freeze3xup_chr-1.gds"))

# make scanAnnot match sample.ids in GDS file
sample.id <- read.gdsn(index.gdsn(gdata, "sample.id"))
genotype.include = (1:length(sample.id))
dat <- data.frame(scanID=sample.id[genotype.include], stringsAsFactors=FALSE)
dat <- dplyr::left_join(dat, pData(scanAnnot_tmp), by='scanID') 
pData(scanAnnot_tmp) <- dat



# gdsSubset(parent.gds='SOL_freeze3xup_chr-1.gds', sub.gds='subset_data/SOL_freeze3xup_chr-1.gds',
#           sample.include=sample.id[genotype.include],compress = 'ZIP',
#           verbose=TRUE)
# 
# 
# gdata <- openfn.gds(paste0(filepath, "/thornton/HCHS_SOL/subset_data/SOL_freeze3xup_chr-1.gds"))



gds <- GdsGenotypeReader(filename = gdata)


cat('start association testing \n')

if(T){
time_assotest_prep = system.time({
genoData <- GenotypeData(gds, scanAnnot=scanAnnot_tmp)
iterator <- GenotypeBlockIterator(genoData)
})


time_assotest = system.time({
assoc <- assocTestSingle(iterator, nullMM, test = "Score")
})

save(list=c('assoc','time_assotest', 'time_assotest_prep'), file="chr1.RData")
}



nullMMREHE = nullMMREHE[-c(20, 22,27)]


genoData <- GenotypeData(gds, scanAnnot=scanAnnot_tmp)
iterator <- GenotypeBlockIterator(genoData)

time_assotest2 = system.time({
assocREHE <- assocTestSingle(iterator, nullMMREHE, test = "Score")
})
cat('association testing done \n')


save(list=c('assocREHE', 'time_assotest2'), file="chr1REHE.RData")



#write.table(file="/results/chr1.txt", sep = " ", quote = F, row.names = F, assoc)

#write.table(file="/results/chr1REHE.txt", sep = " ", quote = F, row.names = F, assocREHE)


save.image(file='assotestchr1resREHE.RData')

close(gds)







