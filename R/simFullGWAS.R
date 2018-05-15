## simGWAS helper functions

library(optparse)
library(data.table)
library(magrittr)

## atm set to true untill fully validated
TEST<-TRUE
option_list = list(
        #make_option(c("-l", "--ldblock"), type="character", default=NULL,
        #      help="ld block to process", metavar="character"),
        make_option(c("-s", "--scenario_file"), type="character", default=NULL,
              help="scenario file to use", metavar="character"),
        #make_option(c("-o", "--out_dir"), type="character", default=NULL,
        #      help="output directory", metavar="character"),
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
      #ldblock='77',
      #scenario_file="/home/ob219/rds/hpc-work/simBasis/support/scenarios/scenario1.yml",
      scenario_file="/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios/scenario1_detec.yml",
      #out_dir="/home/ob219/rds/hpc-work/simBasis/simulations/",
      ld_source='jp_ld',
      nsim=10
  )
}

print(args)

library(data.table)
library(devtools)
install_github("chr1swallace/simGWAS")
install_github("ollyburren/cupcake")
library(simGWAS)
library(cupcake)
library(mvtnorm)
library(yaml)
library(ggplot2)
library(cowplot)



BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='/rds/project/cew54/rds-cew54-wallace-share/Data/reference/UK10K/chr1.bcf.gz'
SNP_SUPPORT='/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/support/chr1_maf_0.01_pCV.RDS'
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

gethapPos <- function(DT,scenarios,bcf_tools=BCF_TOOLS,bcf_file=BCF_FILE,chr='chr1',bcf_maf=0.01) {
  all.cvs <- lapply(scenarios,'[[','CV') %>% do.call('c',.)
  cv <- DT[id %in% all.cvs,]$id
  basis_variants <- c(DT[asb==TRUE,]$id,cv) %>% unique
  rfile <- tempfile()
  write(gsub(":","\t",basis_variants),file=rfile)
  cmd <- sprintf(
    "%s view %s --min-af %f:minor --max-alleles 2 --min-alleles 2 -R %s -Ov",
    bcf_tools,bcf_file,bcf_maf,rfile)
  message(cmd)
  y=fread(cmd)
  unlink(rfile)
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

simCVLDBlock <- function(snps.DT,scenarios,nsim){
  h <- gethapPos(snps.DT,scenarios)
  variants <- colnames(h)
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
    ## we supply OR but software needs log(or)
    gammaCV <- S$gammaCV[S$CV==CV]
    if(length(CV)==0){
      # under the null
      z <- t(rmvnorm(n = nsim, mean = rep(0,length(variants)), sigma = LD))
      #vbeta <- vbeta(S$N0,S$N1,snps.DT[id %in% variants]$MAF)^2
      #beta <- z * sqrt(vbeta)
      vbeta <- vbeta(S$N0,S$N1,snps.DT[id %in% variants]$MAF)
      beta <- z * sqrt(vbeta)
      keep <- which(variants %in% snps.DT[asb==TRUE,]$id)
      return(list(z=z[keep,],beta=beta[keep,],snps=variants[keep],N0=S$N0,N1=S$N1,cv=character(),gammaCV=numeric()))
    }else{
      # under the alternative
      FP <- simGWAS::make_GenoProbList(snps=variants,W=CV,freq=freq)
      z <- simulated_z_score(N0=S$N0,
                          N1=S$N1,
                          snps=variants,
                          W=CV,
                          gamma.W=gammaCV,
                          freq=freq,
                          GenoProbList=FP,
                          nrep=nsim) %>% t
      tmpsimv <- simulated_vbeta(N0=S$N0,
                          N1=S$N1,
                          snps=variants,
                          W=CV,
                          gamma.W=gammaCV,
                          freq=freq,
                          GenoProbList=FP,
                          nrep=nsim) %>% t
      tmpsimv <- 1/tmpsimv
      beta <- z * sqrt(tmpsimv)
      keep <- which(variants %in% snps.DT[asb==TRUE,]$id)
      list(z=z[keep,],beta=beta[keep,],snps=variants[keep],N0=S$N0,N1=S$N1,cv=CV,gammaCV=gammaCV)
    }
  })
  return(res)
}

# code to unpack simulations and create DT for basis creation

buildDT <- function(simlist,snp.DT){
  tmp.DT <- lapply(names(simlist),function(bname){
    b <- simlist[[bname]]
    #message(names(b))
    lapply(names(b),function(gw){
      s <- b[[gw]]
      ns <- ncol(s[['beta']])
      or <- as.vector(exp(s[['beta']]))
      ## if z is above a certain value then we get p=0
      ## this messes up downstream analysis
      p <- as.vector(pnorm(abs(s[['z']]),lower.tail=FALSE) * 2)
      #p <- (pnorm(abs(s[['z']]),lower.tail=FALSE,log=TRUE) + log(2)) %>% as.vector %>% exp
      pid <- rep(s[['snps']],ns)
      n <- rep(s$N1 + s$N0,ns)
      sim <- rep(1:ns,each=length(s[['snps']]))
      data.table(pid=pid,or=or,p.value=p,trait=gw,sim=sim,n=n,n1=s$N1,ld.block=as.numeric(bname))
    }) %>% rbindlist
  }) %>% rbindlist
  setkey(tmp.DT,pid)
  maf <- snp.DT[,.(pid=id,maf=MAF)]
  setkey(maf,pid)
  tmp.DT <- tmp.DT[maf]
  split(tmp.DT[,.(pid,or,p.value,trait,sim,n,n1,ld.block,maf)],tmp.DT$sim)
}



snps <- readRDS(SNP_SUPPORT)
## check to see if there are any ld blocks that have no basis snps in as we don't
## need to simulate these
use <- snps[,list(has_asb=any(asb)),by=`ld_source`][has_asb==TRUE,][[ld_source]]
snps <- snps[get(`ld_source`) %in% use,]

## read in scenarios from YML
scen <- read_yaml(args$scenario_file)
## do simulations
all.sims <- lapply(split(snps,snps[[ld_source]]),function(DT){
  tmp <- simCVLDBlock(DT,scen,args$nsim)
  names(tmp) <- names(scen)
  tmp
})

## build data structure for basis generation
DT.sims<-buildDT(all.sims,snps)
basis.sims <- names(scen)[(sapply(scen,'[[','basis'))]
proj.sims <- names(scen)[!(sapply(scen,'[[','basis'))]

## for a given set of basis sims we can do projection at the same time

basis <- basis.sims
proj <- proj.sims



## test plotting code to check simulations make sense
if(FALSE
  DT <- DT.sims[[2]]
  ## compute z values for qqplot
  DT[,z:=qnorm(p.value/2,lower.tail=FALSE) * sign(log(or))]
  DT[,pos:=as.numeric(sub("chr1:","",pid,fixed=TRUE))]
  ## annotate causal variants for each study

  DT[,CV:=FALSE]

  for(na in names(scen)){
    message(na)
    cvs <- scen[[na]]$CV
    lds <- unique(snps[id %in% cvs,]$jp_ld)
    DT[(ld.block %in%  lds) & trait==na,]$CV <- TRUE
  }

  ## Manhattan
  ggplot(DT[grep("GWAS10",trait),],aes(x=pos,y=z,col=CV)) +
  geom_point() + facet_wrap(~trait,ncol=1) + geom_hline(yintercept=c(-log10(5e-8),log10(5e-8)),col='firebrick1',lty=2)

  ## qqplot
  colmat <- ifelse(DT[trait=="GWAS10",]$CV,'red','black')
  qqnorm(DT[trait=="GWAS10",]$z,col=colmat)
  #dcast(DT,pid~or+trait)
}


## create a basis from a simulation and project on a set of linked GWAS scenarios

createBasisAndProj<-function(DT,basis,proj){
  basis.idx <- which(DT$trait %in% basis)
  proj.idx <- which(DT$trait %in% proj)


  ## helper func that build matrix for PC input where no shrinkage e.g. Z and beta

  build_matrix_raw <- function(DT){
    B <- dcast(DT,pid ~ trait,value.var='metric')
    snames <- B[,1]$pid
    tmp.mat <- as.matrix(B[,-1]) %>% t()
    colnames(tmp.mat) <- snames
    tmp.mat
    #tmp.mat <- rbind(tmp.mat,control=rep(0,ncol(tmp.mat)))
    #prcomp(tmp.mat,center=TRUE,scale=FALSE)
  }


  ## helper func that build matrix for PC input where some weighting shrinkage is involved

  build_matrix_shrink <- function(bDT,sDT,vmethod){
    message(sprintf("Using %s",vmethod))
    stmp<-sDT[,c('pid',vmethod),with=FALSE]
    tmp<-bDT[stmp]
    tmp$metric <- tmp[[vmethod]] * log(tmp$or)
    B <- dcast(tmp,pid ~ trait,value.var='metric')
    snames <- B[,1]$pid
    tmp.mat <- as.matrix(B[,-1]) %>% t()
    colnames(tmp.mat) <- snames
    tmp.mat
    #tmp.mat <- rbind(tmp.mat,control=rep(0,ncol(tmp.mat)))
    #prcomp(tmp.mat,center=TRUE,scale=FALSE)
  }

  ## beta
  RES <- list()
  doRAW <- function(){
    b <- build_matrix_raw(DT[basis.idx,])
    ## add control
    b <- rbind(b,control=rep(0,ncol(b)))
    ## build basis
    pc <- prcomp(b,center=TRUE,scale=FALSE)
    list(basis=pc,proj=predict(pc,newdata=build_matrix_raw(DT[proj.idx,])))
  }
  DT[,metric:=log(or)]
  RES[['beta']] <- doRAW()
  DT[,metric:=sign(log(or)) * qnorm(p.value/2,lower.tail=FALSE)]
  RES[['z']] <- doRAW()
  ## compute the shrinkage
  shrink.DT <- cupcake::compute_shrinkage_metrics(DT[basis.idx,])
  shrink.DT[,c('r_emp_maf_se','r_est_maf_se'):=list(1/emp_maf_se,1/est_maf_se)]
  setkey(shrink.DT,'pid')
  doSHRINK <- function(metric){
    b <- build_matrix_shrink(DT[basis.idx,],shrink.DT,metric)
    ## add control
    b <- rbind(b,control=rep(0,ncol(b)))
    ## build basis
    pc <- prcomp(b,center=TRUE,scale=FALSE)
    list(basis=pc,proj=predict(pc,newdata=build_matrix_shrink(DT[proj.idx,],shrink.DT,metric)))
  }
  metrics <- c('r_est_maf_se','r_emp_maf_se','est_shrinkage','emp_shrinkage','ws_est_shrinkage','ws_emp_shrinkage')
  for(m in metrics){
    RES[[m]] <- doSHRINK(m)
  }
  return(RES)
}

stop()

# code to plot bi plots for an example simulation to check things

if(FALSE){
  library(data.table)
  library(magrittr)

  sim1 <- createBasisAndProj(DT.sims[[1]],basis.sims,proj.sims)

  biplot.DT <- lapply(names(sim1),function(n){
    S<-sim1[[n]]
    tDT <- data.table(rbind(S$basis$x,S$proj))
    tDT[,trait:=c(rownames(S$basis$x),rownames(S$proj))]
    tDT[,metric:=n]
    tDT
  }) %>% rbindlist

  library(ggrepel)

  ggplot(biplot.DT,aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel() + facet_wrap(~metric)
}

## create distance metrics

getDist <- function(S,ref='control'){
  R <- S$basis$x[ref,]
  tDT <- data.table(reference=ref,t(apply(S$proj,1,function(x) sqrt(sum((x - R)^2)))))
}

## compute dist matrix for all simulations considered
dist.res <- lapply(DT.sims,function(db){
  ts <- createBasisAndProj(db,basis.sims,proj.sims)
  ctres <- lapply(ts,getDist) %>% rbindlist
  ctres[,metric:=names(ts)]
  ares <- lapply(ts,getDist,'GWAS10') %>% rbindlist
  ares[,metric:=names(ts)]
  rbind(ctres,ares)
})


all.distances <- rbindlist(dist.res)

mall <- melt(all.distances,id.vars=c('metric','reference'))

mall[,scale.value:=scale(value),by=c('reference','metric')]

mall[,variable:=factor(variable,levels=c('GWAS10_500_3000','GWAS10_2000_2000','GWAS10_5000_5000','share_500_3000','share_2000_2000','share_5000_5000','random_500_3000','random_2000_2000','random_5000_5000'))]


pp<-ggplot(mall,aes(x=variable,y=value,col=reference)) + geom_boxplot() +
facet_wrap(~metric,scale="free",ncol=4) + theme(axis.text.x=element_text(angle = -90, hjust = 0))

pb<-ggplot(mall,aes(x=variable,y=scale.value,col=reference)) + geom_boxplot() +
facet_wrap(~metric,scale="free",ncol=4) + theme(axis.text.x=element_text(angle = -90, hjust = 0))




stop()


## BELOW HERE IS EXPERIMENTAL CODE FOR malhalanobis distance stuff.


## I think that this is equivalent to the malhalanobis distance
## here we multiply the loading for a given axis by the variance explained by that

mahalanobis_pairwise<-function(pc,proj){
  ## add in zero for mahalanobis distance
  load <- pc$x
  den <- pc$sdev^2
  if(exists("proj"))
    load<-rbind(load,proj)
  ## pairwise distance function
  #dist<-apply(load,1,function(p){
  #  apply(pc$x,1,function(x) (x-p)^2/pc$sdev^2) %>% colSums %>% sqrt
  #})
  ## equivalent and faster
  ## note we can't use scale if we add in the projections - we assume centred
  ## to standardise divide through by the sdev
  stand <- apply(load,1,function(x) x/pc$sdev) %>% t
  as.matrix(dist(stand))
  #mahalanobis
  #man <- apply(pc$x,1,function(x) x^2/den) %>% colSums
}


## not sure that this is what we want as it massively inflates distances for
## components that explain a tiny amount of variance but I need to read more about it
## see code in mahalanobis.R for an explanation
mahalanobis_pairwise(sim1$emp_shrinkage$basis,sim1$emp_shrinkage$proj)


## check to see what is going on using a biplot

pd <- lapply(names(sim1),function(n){
  D <- sim1[[n]]
  m <- rbind(D$basis$x,D$proj)
  m.DT <- data.table(trait=rownames(m),m)
  m.DT[,stat:=n]
}) %>% rbindlist

library(ggrepel)
pd.f <- pd[grep("GWAS",trait),]
ggplot(pd.f[grep("emp_shrinkage",stat),],aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel() + facet_wrap(~stat,scales="free")


ggplot(pd[stat=='ws_emp_shrinkage',],aes(x=PC1,y=PC2,label=trait)) + geom_point() + geom_text_repel()

## beta

test[,metric:=log(or)]
pc.beta <- build_pca_special(test)

## z

test[,metric:=sign(log(or)) * qnorm(p.value/2,lower.tail=FALSE)]
pc.z <- build_pca_special(test)

## these need the shrinkage

shrink.DT <- cupcake::compute_shrinkage_metrics(test)
shrink.DT[,c('r_emp_maf_se','r_est_maf_se'):=list(1/emp_maf_se,1/est_maf_se)]
setkey(shrink.DT,'pid')

# here e is for estimated a is for actual - empirical version

## gamma hat using analytical estimate MAF SE

pc.e.ghat <- build_pca_shrink(test,shrink.DT,'r_est_maf_se')

## gamma hat using empirical MAF SE

pc.a.ghat <- build_pca_shrink(test,shrink.DT,'r_emp_maf_se')

## shrinkage 1  using analytical estimate MAF SE

pc.e.sh1 <- build_pca_shrink(test,shrink.DT,'est_shrinkage')

## shrinkage 1 using empirical MAF SE

pc.a.sh1 <- build_pca_shrink(test,shrink.DT,'emp_shrinkage')

## shrinkage 2 using analytical estimate MAF SE

pc.e.sh2 <- build_pca_shrink(test,shrink.DT,'ws_est_shrinkage')

## shrinkage 2 using empirical MAF SE

pc.a.sh2 <- build_pca_shrink(test,shrink.DT,'ws_emp_shrinkage')


## see how this looks



## beta

create_ds_matrix <- function(bDT,sDT,method=c('emp','est','memp','mest')){
  if(missing(method)){
    method='emp'
  }
  message(sprintf("Using %s",method))
  vmethod = sprintf("%s_shrinkage",method)
  stmp<-sDT[,c('pid',vmethod),with=FALSE]
  tmp<-bDT[stmp]
  tmp$metric <- tmp[[vmethod]] * log(tmp$or)
  B <- dcast(tmp,pid ~ trait,value.var='metric')
  snames <- B[,1]$pid
  tmp.mat <- as.matrix(B[,-1]) %>% t()
  colnames(tmp.mat) <- snames
  return(tmp.mat)
}


buildBasis <- function(pno){
  tmp.DT <- lapply(names(all.sims),function(bname){
    b <- all.sims[[bname]]
    lapply(names(b),function(gw){

      s <- b[[gw]]
      or <- exp(s[['beta']][,pno])
      p <- pnorm(abs(s[['z']][,pno]),lower.tail=FALSE) * 2
      pid <- s[['snps']]
      n <- s$N1 + s$N0
      data.table(pid=pid,or=or,p.value=p,trait=gw,n=n,n1=s$N1,ld.block=as.numeric(bname))
    }) %>% rbindlist
  }) %>% rbindlist
  ## need to add MAF
  setkey(tmp.DT,pid)
  basis.DT <- tmp.DT[maf.DT]
  shrink.DT<-cupcake::compute_shrinkage_metrics(basis.DT)
  basis.mat.emp <- cupcake::create_ds_matrix(basis.DT,shrink.DT)
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
}


## this unpacks as a long list
pd <- lapply(all.sims,function(b){
  lapply(names(b),function(gw){
    s <- b[[gw]]
    tmp <- data.table(gwas=gw,variant=s[['snps']],mlp=-log10(pnorm(abs(s[['z']][,1]),lower.tail=FALSE) * 2))
    tmp[,CV:=variant %in% s$cv]
  }) %>% rbindlist
}) %>% rbindlist

pd[,pos:=as.numeric(sub("chr1:","",variant,fixed=TRUE))]

library(ggplot2)
ggplot(pd,aes(x=pos,y=mlp)) + geom_point()  +
geom_vline(data=cv,aes(xintercept=as.numeric(pos))) + facet_wrap(~gwas,ncol=1)

lapply(all.sims)
