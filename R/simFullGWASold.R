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

getLD <- function(snps.DT,source=c('jp_ld','hm_ld'),bcf_tools=BCF_TOOLS,bcf_file=BCF_FILE,chr='chr1',bcf_maf=0.01){
  tmp<-snps.DT[,list(start1=min(start),stop1=max(start)),by=`source`]
  setnames(tmp,`source`,'blocknum')
  tmp[,comm:=sprintf(
    "%s view %s --min-af %f:minor --max-alleles 2 --min-alleles 2 -r %s:%s-%s -Ov",
    bcf_tools,bcf_file,bcf_maf,chr,start1,stop1)]
  setnames(tmp,c('start1','stop1'),c('start','stop'))
  tmp
}


gethap <- function(cmd) {
    y=fread(cmd)
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste(y[['#CHROM']],y$POS,sep=':')
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}


vbeta <- function(N0,N1,f){
  sqrt(1/2) * sqrt((N0+N1)/(N0*N1)) * sqrt(1/f + 1/1-f)
}

nullSim <- function(ld.cmd,snps.DT,nsim=NSIM){
  h <- gethap(ld.cmd)
  variants <- colnames(h)
  use <- apply(h,2,var)>0
  h <- h[,use,drop=FALSE]
  LD <- cor2(h)
  z <- t(rmvnorm(n = nsim, mean = rep(0,length(variants)), sigma = LD))
  vbeta <- vbeta(N0,N1,snps.DT[id %in% variants]$MAF)^2
  beta <- z * sqrt(vbeta)
  list(z=z,beta=beta,snps=variants,cv=character(),gammaCV=numeric())
}

simCVLDBlock <- function(ld.cmd,CV,gammaCV,nsim=NSIM){
  h <- gethap(ld.cmd)
  variants <- colnames(h)
  use <- apply(h,2,var)>0
  h <- h[,use,drop=FALSE]
  freq <- as.data.frame(h+1)
  freq$Probability <- 1/nrow(freq)
  LD <- cor2(h)
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
  list(z=z,beta=beta,snps=variants,cv=CV,gammaCV=gammaCV)
}

# where we get blocks from
source <- 'jp_ld'
snps <- readRDS("/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS")
## check to see if there are any ld blocks that have no basis snps in as we don't
## need to simulate these
use <- snps[,list(has_asb=any(asb)),by=`source`][has_asb==TRUE,][[source]]
snps <- snps[get(`source`) %in% use,]
ldd <- getLD(snps,source)
ldd <- split(ldd$comm,ldd$blocknum)


library(magrittr)
nCV <- 13
CV <- sapply(split(snps$id,snps[[source]]),sample,size=1) %>% sample(.,nCV)
effect <- seq(1.2,1.6,by=0.1) %>% sample(.,nCV,replace=TRUE)
## create a lookup mapping CV and effect to LD block
idx <- which(snps$id %in% CV)
ld.lu.cv <- split(snps[idx,]$id,snps[[source]][idx])
ld.lu.gamma <- split(effect,snps[[source]][idx])

## code for simulating a block under alternative

## benchmarking stuff
ld.block <- "35"
library(microbenchmark)
mbm <- microbenchmark(
  alt=simCVLDBlock(ldd[[ld.block]],ld.lu.cv[[ld.block]],ld.lu.gamma[[ld.block]]),
  null=nullSim(ldd[[ld.block]],snps),
  times=50
)
library(ggplot2)
autoplot(mbm)

## build null first
alt.blocks <- names(ldd)[!names(ldd) %in% names(ld.lu.cv)]

alt <- lapply(alt.blocks,function(ld.block){
  message(ld.block)
  nullSim(ldd[[ld.block]],snps)
})

## very slow - the slow bit is getting GT and computing LD matrix ?
## read in all scenarios and loop through.

#scenario is a list of CV and their associated effect sizes

library(yaml)


simCVLDBlock <- function(ld.cmd,snps.DT,scenarios,nsim=NSIM,prune=TRUE){
  h <- gethap(ld.cmd)
  variants <- colnames(h)
  use <- apply(h,2,var)>0
  h <- h[,use,drop=FALSE]
  freq <- as.data.frame(h+1)
  freq$Probability <- 1/nrow(freq)
  LD <- cor2(h)
  ## loop through scenarios
  res <- lapply(names(scenarios),function(gw){
    message(paste("Simulation of",gw))
    S <- scenarios[[gw]]
    if(length(S$nsim)!=0)
      nsim <- S$nsim
    scen.name <- S$scenario
    CV <- variants[variants %in% S$CV]
    gammaCV <- S$gammaCV[S$CV==CV]
    if(length(CV)==0){
      # under the null
      z <- t(rmvnorm(n = nsim, mean = rep(0,length(variants)), sigma = LD))
      #vbeta <- vbeta(S$N0,S$N1,snps.DT[id %in% variants]$MAF)^2
      #beta <- z * sqrt(vbeta)
      vbeta <- vbeta(S$N0,S$N1,snps.DT[id %in% variants]$MAF)
      beta <- z * sqrt(vbeta)
      if(prune){
        keep <- which(variants %in% snps.DT[asb==TRUE,]$id)
        return(list(z=z[keep,],beta=beta[keep,],snps=variants[keep],N0=S$N0,N1=S$N1,cv=character(),gammaCV=numeric()))
      }
      return(list(z=z,beta=beta,snps=variants,N0=S$N0,N1=S$N1,cv=character(),gammaCV=numeric()))
    }else{
      # under the alternative
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
      if(prune){
        keep <- which(variants %in% snps.DT[asb==TRUE,]$id)
        return(list(z=z[keep,],beta=beta[keep,],snps=variants[keep],N0=S$N0,N1=S$N1,cv=CV,gammaCV=gammaCV))
      }
      return(list(z=z,beta=beta,snps=variants,N0=S$N0,N1=S$N1,cv=CV,gammaCV=gammaCV))
    }
  })
}

ld.block <- "42"
scenarios <- read_yaml("/home/ob219/rds/hpc-work/simBasis/support/scenarios/Scenario1.yml")
foo<-simCVLDBlock(ldd[[ld.block]],snps,scenarios,10)
## possibly don't prune as this gives us the opportunity to see what happens when we
## upsample the number of variants 
