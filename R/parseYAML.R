library(yaml)
yaml.load_file("yaml/scenario1.yml")

## can we sim a gwas with a list of CV


#!/usr/bin/env Rscript
library(data.table)
#library(randomFunctions)
library(devtools)
install_github("chr1swallace/simGWAS")
library(simGWAS)
library(mvtnorm)
library(corpcor)

BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='/home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/chr1.bcf.gz'


args <- getArgs(defaults=list(N=1000,NEFF=1,NSIM=100,
                              file.ldd="/home/cew54/newscratch/Data/reference/lddetect/EUR/fourier_ls-chr22.bed",
                              file.vcf="/home/cew54/newscratch/Data/reference/UK10K/chr22.bcf.gz"),
                numeric=c("N","NEFF","NSIM"))

## get everything from annotation file (including ld blocks)

snps <- readRDS("/home/ob219/rds/hpc-work/simBasis/support/chr1.RDS")
snps[,id:=paste(chr,start,sep=':')]

getLD <- function(source=c('jp_ld','hm_ld')){
  tmp<-snps[,list(start1=min(start),stop1=max(start)),by=`source`]
  tmp[,chr:=unique(snps$chr)]
  setnames(tmp,`source`,'blocknum')
  tmp[,comm:=sprintf(
    "%s view %s --min-af 0.01:minor --max-alleles 2 --min-alleles 2 -r %s:%s-%s -Ov",
    BCF_TOOLS,BCF_FILE,chr,start1,stop1)]
  setnames(tmp,c('start1','stop1'),c('start','stop'))
  tmp
}

ldd <- getLD('jp_ld')
## working directory
d <- "/rds/user/cew54/hpc-work/simgwas"


gethap <- function(i,chr) {
    y=fread(ldd$comm[i])
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste(y[['#CHROM']],y$POS,sep=':')
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}

setCV <- 'chr1:1866893'
setGammaCV <- 1.2
N0 <- 1000
N1 <- 1000
NSIM <- 10

simv <- simz <- simbeta <- vector("list",nrow(ldd))
for(i in 1:nrow(ldd)) {
    h <- gethap(i)
    use <- apply(h,2,var)>0
    h <- h[,use,drop=FALSE]
    freq <- as.data.frame(h+1)
    freq$Probability <- 1/nrow(freq)
    LD <- cor2(h)
    variants <- colnames(h)
    CV <- variants[which(variants %in% setCV)]
    gammaCV <- log(setGammaCV[setCV %in% CV])
    if(length(CV)==0){
      CV <- sample(variants,1)
      gammaCV <- 0
    }
    vbeta(N0,N1,snps[id %in% variants]$MAF)^2
    ## method 1 - simulate Z scores and adjust by simulated variance to get beta
    EZ <- est_statistic(N0=N0, # number of controls
                        N1=N1, # number of cases
                        snps=variants, # column names in freq of SNPs for which Z scores should be generated
                        W=CV, # causal variants, subset of snps
                        gamma1=gammaCV, # odds ratios
                        freq=freq, # reference haplotypes
                        GenoProbList=FP) # FP above
    simz[[i]] <- t(rmvnorm(n = args$NSIM, mean = EZ, sigma = LD))
    tmpsimv <- sim_vbeta(N0=args$N, # number of controls
                         N1=args$N, # number of cases
                         snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                         W=CV, # causal variants, subset of snps
                         gamma1=log(g1), # odds ratios
                         freq=freq, # reference haplotypes
                         GenoProbList=FP,
                         nsim=args$NSIM)
    simv[[i]] <- 1/do.call("cbind",tmpsimv)
    simbeta[[i]] <- simz[[i]] * sqrt(simv[[i]])
}


vbeta <- function(N0,N1,f){
  sqrt(1/2) * sqrt((N0+N1)/(N0*N1)) * sqrt(1/f + 1/1-f)
}
## altmethod

simvalt <- simzalt <- simbetaalt <- vector("list",nrow(ldd))

for(i in 1:nrow(ldd)) {
    h <- gethap(i)
    use <- apply(h,2,var)>0
    h <- h[,use,drop=FALSE]
    freq <- as.data.frame(h+1)
    freq$Probability <- 1/nrow(freq)
    LD <- cor2(h)
    variants <- colnames(h)
    CV <- variants[which(variants %in% setCV)]
    gammaCV <- log(setGammaCV[setCV %in% CV])
    if(length(CV)==0){
      simzalt[[i]] <- t(rmvnorm(n = NSIM, mean = rep(0,length(variants)), sigma = LD))
      tmpsimvy <- vbeta(N0,N1,snps[id %in% variants]$MAF)^2
      ## chris seems to generate the reciprocal of the variance ?
      simvalt[[i]] <- tmpsimvy
      mybeta <- simzalt[[i]] * sqrt(simvalt[[i]])
      #simbetaalt[[i]] <- simzalt[[i]] * sqrt(simvalt[[i]])
    }

    nulltime <- function(){
      simzalt[[i]] <- t(rmvnorm(n = NSIM, mean = rep(0,length(variants)), sigma = LD))
      tmpsimvy <- vbeta(N0,N1,snps[id %in% variants]$MAF)^2
      ## chris seems to generate the reciprocal of the variance ?
      simvalt[[i]] <- tmpsimvy
      mybeta <- simzalt[[i]] * sqrt(simvalt[[i]])
    }

    christime <- function(){
      FP <- make_GenoProbList(snps=colnames(h),W=CV,freq=freq)
      ## method 1 - simulate Z scores and adjust by simulated variance to get beta
      EZ <- est_statistic(N0=N0, # number of controls
                          N1=N1, # number of cases
                          snps=variants, # column names in freq of SNPs for which Z scores should be generated
                          W=CV, # causal variants, subset of snps
                          gamma1=gammaCV, # odds ratios
                          freq=freq, # reference haplotypes
                          GenoProbList=FP) # FP above
      simzalt[[i]] <- t(rmvnorm(n = NSIM, mean = EZ, sigma = LD))
      tmpsimv <- sim_vbeta(N0=N0, # number of controls
                           N1=N1, # number of cases
                           snps=variants, # column names in freq of SNPs for which Z scores should be generated
                           W=CV, # causal variants, subset of snps
                           gamma1=gammaCV, # odds ratios
                           freq=freq, # reference haplotypes
                           GenoProbList=FP,
                           nsim=NSIM)
      simvalt[[i]] <- 1/do.call("cbind",tmpsimv)
      chrisbeta <- simzalt[[i]] * sqrt(simvalt[[i]])
    }

    system.time(christime())
    system.time(nulltime())



    FP <- make_GenoProbList(snps=colnames(h),W=CV,freq=freq)
    ## method 1 - simulate Z scores and adjust by simulated variance to get beta
    EZ <- est_statistic(N0=N0, # number of controls
                        N1=N1, # number of cases
                        snps=variants, # column names in freq of SNPs for which Z scores should be generated
                        W=CV, # causal variants, subset of snps
                        gamma1=gammaCV, # odds ratios
                        freq=freq, # reference haplotypes
                        GenoProbList=FP) # FP above
    simzalt[[i]] <- t(rmvnorm(n = NSIM, mean = EZ, sigma = LD))
    tmpsimv <- sim_vbeta(N0=N0, # number of controls
                         N1=N1, # number of cases
                         snps=variants, # column names in freq of SNPs for which Z scores should be generated
                         W=CV, # causal variants, subset of snps
                         gamma1=gammaCV, # odds ratios
                         freq=freq, # reference haplotypes
                         GenoProbList=FP,
                         nsim=NSIM)
    simvalt[[i]] <- 1/do.call("cbind",tmpsimv)
    chrisbeta <- simzalt[[i]] * sqrt(simvalt[[i]])
    simbetaalt[[i]] <- simzalt[[i]] * sqrt(simvalt[[i]])
}



## REFACTORED CODE BELOW HERE

## simGWAS helper functions


library(data.table)
library(devtools)
install_github("chr1swallace/simGWAS")
library(simGWAS)
library(mvtnorm)


BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='/home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/chr1.bcf.gz'

setCV <- 'chr1:1866893'
setGammaCV <- 1.2
N0 <- 1000
N1 <- 1000
NSIM <- 10

## note minor allele freq is hard coded

getLD <- function(snps.DT,source=c('jp_ld','hm_ld'),bcf_tools=BCF_TOOLS,bcf_file=BCF_FILE,chr='chr1'){
  tmp<-snps.DT[,list(start1=min(start),stop1=max(start)),by=`source`]
  setnames(tmp,`source`,'blocknum')
  tmp[,comm:=sprintf(
    "%s view %s --min-af 0.01:minor --max-alleles 2 --min-alleles 2 -r %s:%s-%s -Ov",
    bcf_tools,bcf_file,chr,start1,stop1)]
  setnames(tmp,c('start1','stop1'),c('start','stop'))
  tmp
}

snps <- readRDS("/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS")

ldd <- getLD(snps,'jp_ld')

nullSim <- function(variants,sigma.LD,snps.DT,nsim=NSIM){
  z <- t(rmvnorm(n = nsim, mean = rep(0,length(variants)), sigma = LD))
  vbeta <- vbeta(N0,N1,snps.DT[id %in% variants]$MAF)^2
  beta <- z * sqrt(vbeta)
  list(z=z,beta=beta,snps=variants)
}

cvSim <- function(cv,variants,freq,LD,g1,nsim=NSIM){
  FP <- make_GenoProbList(snps=variants,W=CV,freq=freq)
  EZ <- est_statistic(N0=N0, # number of controls
                      N1=N1, # number of cases
                      snps=variants, # column names in freq of SNPs for which Z scores should be generated
                      W=CV, # causal variants, subset of snps
                      gamma1=gammaCV, # odds ratios
                      freq=freq, # reference haplotypes
                      GenoProbList=FP) # FP above
  z <- t(rmvnorm(n = NSIM, mean = EZ, sigma = LD))
  tmpsimv <- sim_vbeta(N0=N0, # number of controls
                       N1=N1, # number of cases
                       snps=variants, # column names in freq of SNPs for which Z scores should be generated
                       W=CV, # causal variants, subset of snps
                       gamma1=gammaCV, # odds ratios
                       freq=freq, # reference haplotypes
                       GenoProbList=FP,
                       nsim=NSIM)
  tmpsimv <- 1/do.call("cbind",tmpsimv)
  beta <- z * sqrt(tmpsimv)
  list(z=z,beta=beta,snps=variants)
}
