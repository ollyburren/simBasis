## Generate Scenario

library(data.table)
library(yaml)
library(magrittr)

## change this
YAML.DIR <- '/home/ob219/rds/rds-cew54-wallace-share/Projects/simBasis/scenarios'

source <- 'jp_ld'
snps <- readRDS("/home/ob219/rds/rds-cew54-wallace-share/Projects/simBasis/support/chr1_maf_0.01_pCV.RDS")
## check to see if there are any ld blocks that have no basis snps in as we don't
## need to simulate these
use <- snps[,list(has_asb=any(asb)),by=`source`][has_asb==TRUE,][[source]]
snps <- snps[get(`source`) %in% use,]
ld.blocks <- unique(snps[[`source`]])

# Scenario 1 basis contains no overlaps but CV's are in 0.8 with at least one basis SNP
scenario.name <- 'no_share'
effect.size <- log(2)
N1 <- 5000
N0 <- 5000
## number of basis GWAS
n.gwas <- 10
## number of CV
cv.length <- 12


gwas <- lapply(split(sample(ld.blocks,size=n.gwas*cv.length),paste0('GWAS',1:n.gwas)),function(x){
  pCV <- paste('possCV',source,sep='_')
  idx <- which(snps[[source]] %in% x & snps[[pCV]]==TRUE)
  CV <- sapply(split(snps[idx,]$id,snps[[source]][idx]),sample,size=1)
  names(CV) <- NULL
  gammaCV <- rep(effect.size,length(CV))
  list(CV=CV,gammaCV=gammaCV,N1=N1,N0=N0,basis=TRUE)
})

cases <- c(500,2000,5000)
controls <- c(3000,2000,5000)

## add a simulation of exactly same GWAS

tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
snamestub <- 'GWAS10_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

## share 50:50 with two of the diseases



tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
tgwas$CV <- c(sample(gwas[['GWAS1']]$CV,cv.length/2),sample(gwas[['GWAS10']]$CV,cv.length/2))
snamestub <- 'share_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

## random CV

tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
pCV <- paste('possCV',source,sep='_')
tgwas$CV <- sample(snps[get(`pCV`),list(CV=sample(id,1)),by=`source`]$CV,cv.length)
snamestub <- 'random_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

write_yaml(gwas,file.path(YAML.DIR,paste(scenario.name,'yml',sep='.')))

## scenario 2 tiled overlap

#########
    #########
        #########

scenario.name <- 'tile_share'
N1 <- 5000
N0 <- 5000
effect.size <- log(2)
n.gwas <- 10

## number of blocks per trait

cv.length <- 12/2

s2.ldblocks <- ld.blocks[1:(cv.length * (n.gwas+1))]

#s2.ldblocks <- head(ld.blocks,length(ld.blocks) %% cv.length * -1)

idx.b <- matrix(s2.ldblocks,ncol=cv.length) %>% t %>% split(.,1:(n.gwas+1))

gwas <- lapply(idx.b,function(x){
  pCV <- paste('possCV',source,sep='_')
  idx <- which(snps[[source]] %in% x & snps[[pCV]]==TRUE)
  CV <- sapply(split(snps[idx,]$id,snps[[source]][idx]),sample,size=1)
  names(CV) <- NULL
  gammaCV <- rep(effect.size,length(CV))
  list(CV=CV,gammaCV=gammaCV,N1=N1,N0=N0,basis=TRUE)
})


gwas2 <- list()
for(i in 1:(length(gwas)-1)){
  gw1 <- gwas[[i]]
  gw2 <- gwas[[i+1]]
  gwas2[[i]] <- list(CV=c(gw1$CV,gw2$CV),gammaCV=c(gw1$gammaCV,gw2$gammaCV),N1=N1,N0=N0,basis=TRUE)
}

gwas <- gwas2
names(gwas) <- paste('GWAS',1:10,sep="")



cases <- c(500,2000,5000)
controls <- c(3000,2000,5000)

## add a simulation of exactly same GWAS

tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
snamestub <- 'GWAS10_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

## share 50:50 with two of the diseases

tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
tgwas$CV <- c(sample(gwas[['GWAS1']]$CV,cv.length),sample(gwas[['GWAS10']]$CV,cv.length))
snamestub <- 'share_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

## random CV

tgwas <- gwas[['GWAS10']]
tgwas$basis <- FALSE
pCV <- paste('possCV',source,sep='_')
tgwas$CV <- sample(snps[get(`pCV`),list(CV=sample(id,1)),by=`source`]$CV,cv.length*2)
snamestub <- 'random_%d_%d'

for(i in 1:length(cases)){
  message(i)
  tgwas$N1 <- cases[i]
  tgwas$N0 <- controls[i]
  sname <- sprintf(snamestub,cases[i],controls[i])
  gwas[[sname]] <- tgwas
}

write_yaml(gwas,file.path(YAML.DIR,paste(scenario.name,'yml',sep='.')))
