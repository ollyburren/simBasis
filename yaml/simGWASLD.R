## simGWAS helper functions

library(optparse)
library(data.table)
library(magrittr)

TEST<-FALSE
option_list = list(
        make_option(c("-l", "--ldblock"), type="character", default=NULL,
              help="ld block to process", metavar="character"),
        make_option(c("-s", "--scenario_file"), type="character", default=NULL,
              help="scenario file to use", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="output directory", metavar="character"),
        make_option(c("-b", "--ld_source"), type="character", default=NULL,
              help="ld block source to use one of jp_ld or hm_ld",metavar="character"),
        make_option(c("-n", "--nsim"), type="numeric", default=10,
              help="Number of simulations to run",metavar="numeric")
        )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}else{
  args <- list(
      ldblock='77',
      scenario_file="/home/ob219/rds/hpc-work/simBasis/support/scenarios/Scenario1.yml",
      out_dir="/home/ob219/rds/hpc-work/simBasis/simulations/",
      ld_source='jp_ld',
      nsim=10
  )
}

print(args)

library(data.table)
library(devtools)
install_github("chr1swallace/simGWAS")
library(simGWAS)
library(mvtnorm)
library(yaml)


BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='/home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/chr1.bcf.gz'
SNP_SUPPORT='/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS'
# THIS DICTATES WHETHER WE DOWNSAMPLE SO THAT ONLY HAVE SNPS IN THE BASIS
#AS_BASIS_FILTER=TRUE
#OUT.DIR <- '/home/ob219/rds/hpc-work/simBasis/simulations/'
#args$scenario_file <- "/home/ob219/rds/hpc-work/simBasis/support/scenarios/Scenario1.yml"
## one of jp_ld or hm_ld - first is j pickrells ld blocks and second is mine using hapmap
#args$ld_source <- 'jp_ld'
## default number of simulations per GWAS - can override in function call
#NSIM <- 10
# need to do this otherwise not nice with data tables
ld_source <- args$ld_source

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

#scenario is a list of CV and their associated effect sizes

simCVLDBlock <- function(ld.cmd,snps.DT,scenarios,nsim,prune=TRUE){
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
      FP <- simGWAS::make_GenoProbList(snps=variants,W=CV,freq=freq)
      #EZ <- est_statistic(N0=S$N0, # number of controls
      #                    N1=S$N1, # number of cases
      #                    snps=variants, # column names in freq of SNPs for which Z scores should be generated
      #                    W=CV, # causal variants, subset of snps
      #                    gamma1=gammaCV, # odds ratios
      #                    freq=freq, # reference haplotypes
      #                    GenoProbList=FP) # FP above

      #expected_z_score<-function(N0,N1,snps,W,gamma.CV,freq)
      #EZ <- expected_z_score(N0=S$N0, # number of controls,
      #  N1=S$N1, # number of cases
      #  snps=variants, # column names in freq of SNPs for which Z scores should be generated
      #  W=CV, # causal variants, subset of snps
      #  gamma.W=gammaCV, # odds ratios
      #  freq=freq,
      #  GenoProbList=FP) # reference haplotype) # FP above
      #z <- t(rmvnorm(n = nsim, mean = EZ, sigma = LD))

      z <- simulated_z_score(N0=S$N0,
                          N1=S$N1,
                          snps=variants,
                          W=CV,
                          gamma.W=gammaCV,
                          freq=freq,
                          GenoProbList=FP,
                          nrep=nsim)
      tmpsimv <- simulated_vbeta(N0=S$N0,
                          N1=S$N1,
                          snps=variants,
                          W=CV,
                          gamma.W=gammaCV,
                          freq=freq,
                          GenoProbList=FP,
                          nrep=nsim)


      #tmpsimv <- sim_vbeta(N0=S$N0, # number of controls
      #                     N1=S$N1, # number of cases
      #                     snps=variants, # column names in freq of SNPs for which Z scores should be generated
      #                     W=CV, # causal variants, subset of snps
      #                     gamma1=gammaCV, # odds ratios
      #                     freq=freq, # reference haplotypes
      #                     GenoProbList=FP,
      #                     nsim=nsim)
      #tmpsimv <- 1/do.call("cbind",tmpsimv)
      tmpsimv <- 1/tmpsimv
      beta <- z * sqrt(tmpsimv)
      if(prune){
        keep <- which(variants %in% snps.DT[asb==TRUE,]$id)
        return(list(z=z[keep,],beta=beta[keep,],snps=variants[keep],N0=S$N0,N1=S$N1,cv=CV,gammaCV=gammaCV))
      }
      return(list(z=z,beta=beta,snps=variants,N0=S$N0,N1=S$N1,cv=CV,gammaCV=gammaCV))
    }
  })
}


snps <- readRDS(SNP_SUPPORT)[get(`ld_source`)==args$ldblock,]
scenarios <- read_yaml(args$scenario_file)

if(nrow(snps)==0)
  stop("No snps ib block !")

ld <- getLD(snps,ld_source,bcf_maf=0.01)
tmp <- simCVLDBlock(ld$comm,snps,scenarios,prune=FALSE)
names(tmp) <- names(scenarios)

fname <- sprintf("%s:%s.RDS",sub("yml$","RDS",
  basename(args$scenario_file)),
  args$ldblock,
  args$nsim)

ofile <- file.path(args$out_dir,fname)
saveRDS(tmp,file=ofile)
message("Sucess")

if(FALSE){
  library(data.table)
  SNP_SUPPORT='/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS'
  snps <- readRDS(SNP_SUPPORT)
  cmd <- "Rscript /home/ob219/git/simBasis/R/simGWASLD.R -l %s  -s /home/ob219/rds/hpc-work/simBasis/support/scenarios/Scenario1.yml -o /home/ob219/rds/hpc-work/simBasis/simulations/ -b jp_ld"
  write(sprintf(cmd,unique(snps$jp_ld)),"/home/ob219/git/simBasis/sh/rs.txt")
}
