## simGWAS helper functions

library(optparse)
library(magrittr)


## atm set to true untill fully validated
TEST<-FALSE

option_list = list(
        make_option(c("-s", "--scenario_file"), type="character", default=NULL,
            help="scenario file to use", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="where to put plots and source data generated", metavar="character"),
        make_option(c("-b", "--ld_source"), type="character", default='jp_ld',
              help="ld block source to use one of jp_ld or hm_ld",metavar="character"),
        make_option(c("-n", "--nsim"), type="numeric", default=100,
              help="Number of simulations to run",metavar="numeric"),
        make_option(c("-p", "--prefix"), action="store_true", default=FALSE,
              help="prefix to add to output files, automatically turns off plotting",metavar="numeric")
)


if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}else{
  args <- list(
      out_dir="/home/ob219/tmp/test_scen",
      #ldblock='77',
      #scenario_file="/home/ob219/rds/hpc-work/simBasis/support/scenarios/scenario1.yml",
      scenario_file="/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios//backbone_share.yml",
      #scenario_file="/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios/scenario2.yml",
      #out_dir="/home/ob219/rds/hpc-work/simBasis/simulations/",
      ld_source='jp_ld',
      nsim=100,
      prefix=FALSE
  )
}

print(args)

PLOT <- TRUE
ofile <- gsub(".yml","",basename(args$scenario_file))
if(args$prefix){
  PLOT <- FALSE
  PREFIX <- sample(letters,6) %>% paste(.,collapse="")
  ofile <- paste(PREFIX,ofile,sep='_')
}

if(file.exists(args$out_dir)){
  out.file.stub <- file.path(args$out_dir,ofile)
}else{
  stop(sprintf("out_dir %s does not exist",args$out_dir))
}


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
library(latex2exp)



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
  sqrt(1/2) * sqrt((N0+N1)/(N0*N1)) * sqrt(1/f + 1/(1-f))
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
  message(sprintf("Processing %d",unique(DT[[ld_source]])))
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
if(PLOT){
  DT <- DT.sims[[1]]
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
  DT[,trait:=factor(trait,levels=paste('GWAS',1:10,sep=''))]
  tp <- ggplot(DT[grep("^GWAS[0-9]+$",trait),],aes(x=pos,y=z,col=CV)) +
  geom_point() + facet_wrap(~trait,ncol=1) +
  geom_hline(yintercept=c(-log10(5e-8),log10(5e-8)),col='firebrick1',lty=2)
  save_plot(paste(out.file.stub,'design','pdf',sep="."),tp,base_height=10,base_width=10)

  #ggplot(DT[-grep("^GWAS[2-9]$",trait),],aes(x=pos,y=z,col=CV)) +
  #geom_point() + facet_wrap(~trait,ncol=1) + geom_hline(yintercept=c(-log10(5e-8),log10(5e-8)),col='firebrick1',lty=2)
  ## qqplot
  #colmat <- ifelse(DT[trait=="GWAS10",]$CV,'red','black')
  #qqnorm(DT[trait=="GWAS10",]$z,col=colmat)
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

  build_matrix_shrink_beta <- function(bDT,sDT,vmethod){
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

  build_matrix_shrink_z <- function(bDT,sDT,vmethod){
    message(sprintf("Using %s",vmethod))
    stmp<-sDT[,c('pid',vmethod),with=FALSE]
    tmp<-bDT[stmp]
    tmp$metric <- tmp[[vmethod]] * tmp$z
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
  #copy z score as we need it later !
  DT[,z:=metric]
  RES[['z']] <- doRAW()
  ## compute the shrinkage
  shrink.DT <- cupcake::compute_shrinkage_metrics(DT[basis.idx,])
  shrink.DT[,c('r_emp_maf_se','r_est_maf_se'):=list(1/emp_maf_se,1/est_maf_se)]
  shrink.DT[,'shrinkage_nog':=bshrink,by=pid]
  shrink.DT[,'ws_shrinkage_nog':=ws_ppi,by=pid]
  setkey(shrink.DT,'pid')
  doSHRINKBeta <- function(metric){
    b <- build_matrix_shrink_beta(DT[basis.idx,],shrink.DT,metric)
    ## add control
    b <- rbind(b,control=rep(0,ncol(b)))
    ## build basis
    pc <- prcomp(b,center=TRUE,scale=FALSE)
    list(basis=pc,proj=predict(pc,newdata=build_matrix_shrink_beta(DT[proj.idx,],shrink.DT,metric)))
  }
  #metrics <- c('r_est_maf_se','r_emp_maf_se','est_shrinkage','emp_shrinkage','ws_est_shrinkage','ws_emp_shrinkage')
  metrics <- c('r_emp_maf_se','emp_shrinkage','ws_emp_shrinkage','shrinkage_nog','ws_shrinkage_nog')
  for(m in metrics){
    n <- 'beta'
    RES[[paste(n,m,sep='_')]] <- doSHRINKBeta(m)
  }
  doSHRINKZ <- function(metric){
    b <- build_matrix_shrink_z(DT[basis.idx,],shrink.DT,metric)
    ## add control
    b <- rbind(b,control=rep(0,ncol(b)))
    ## build basis
    pc <- prcomp(b,center=TRUE,scale=FALSE)
    list(basis=pc,proj=predict(pc,newdata=build_matrix_shrink_z(DT[proj.idx,],shrink.DT,metric)))
  }
  #metrics <- c('r_est_maf_se','r_emp_maf_se','est_shrinkage','emp_shrinkage','ws_est_shrinkage','ws_emp_shrinkage')
  for(m in metrics){
    n <- 'z'
    RES[[paste(n,m,sep='_')]] <- doSHRINKZ(m)
  }
  return(RES)
}


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



getPairDist <- function(S){
  tmp.DT <- dist(S$proj) %>% as.matrix %>% melt %>% data.table
  tmp.DT[,c('qtype','qcase','qctrl'):=tstrsplit(Var1,'_')]
  tmp.DT[,c('mtype','mcase','mctrl'):=tstrsplit(Var2,'_')]
  tmp.DT[,.(qtype,qcase,qctrl,mtype,mcase,mctrl,d=value)]
  #tmp.DT <- tmp.DT[qtype!=mtype & qcase==mcase & qctrl==mctrl,][,.(qtype,mtype,qcase,qctrl,d=value)]
  #tmp.DT[!duplicated(d),]
}

all.pd <- lapply(DT.sims,function(db){
  ts <- createBasisAndProj(db,basis.sims,proj.sims)
  lapply(names(ts),function(n){
    tmp <- getPairDist(ts[[n]])
    tmp[,metric:=n]
    tmp
  }) %>% rbindlist

  #ctres <- lapply(ts,getPairDist) %>% rbindlist
  #ctres[,metric:=names(ts)]
})


all.pd <- rbindlist(all.pd)
all.pd[,label1:=paste(qtype,qcase,qctrl,sep='_')]
all.pd[,label2:=paste(mtype,mcase,mctrl,sep='_')]
if(PLOT){
sum.all.pd <- all.pd[,list(mean.d=median(d)),by=c('label1','label2','metric')]

lab.lev <- c('GWAS10_500_3000','GWAS10_2000_2000',
'GWAS10_5000_5000','share_500_3000','share_2000_2000',
'share_5000_5000','random_500_3000','random_2000_2000','random_5000_5000')

lab.lev2 <- c('ID1','ID2',
'ID3','S1','S2',
'S3','R1','R2','R3')

sum.all.pd[,label1:=factor(label1,levels=lab.lev)]
sum.all.pd[,label2:=factor(label2,levels=lab.lev)]
levels(sum.all.pd$label1)<-lab.lev2
levels(sum.all.pd$label2)<-lab.lev2


all.plots <- lapply(split(sum.all.pd,sum.all.pd$metric),function(dat){
  ggplot(data = dat, aes(x=label1, y=label2, fill=mean.d)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
    geom_tile() + ggtitle(unique(dat$metric)) +
    xlab("Scenario") + ylab("Scenario")
})

ppd <- plot_grid(plotlist=all.plots[c('z','beta','beta_r_emp_maf_se','beta_emp_shrinkage','beta_shrinkage_nog','beta_ws_emp_shrinkage','beta_ws_shrinkage_nog')])
save_plot(paste(out.file.stub,'pw_dist_beta','pdf',sep="."),ppd,base_height=10,base_width=15,base_aspect=1)

ppd <- plot_grid(plotlist=all.plots[c('z','beta','beta_r_emp_maf_se','z_emp_shrinkage','z_shrinkage_nog','z_ws_emp_shrinkage','z_ws_shrinkage_nog')])
save_plot(paste(out.file.stub,'pw_dist_z','pdf',sep="."),ppd,base_height=10,base_width=15,base_aspect=1)
}
saveRDS(all.pd,file=paste(out.file.stub,'pw_dist','RDS',sep="."))

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

#mall[,scale.value:=scale(value),by=c('reference','metric')]

mall[,variable:=factor(variable,levels=c('GWAS10_500_3000','GWAS10_2000_2000',
'GWAS10_5000_5000','share_500_3000','share_2000_2000',
'share_5000_5000','random_500_3000','random_2000_2000','random_5000_5000'))]



#mall[metric=='beta',metric:=latex2exp("$\beta")]

metric_names <- c(
TeX('$\\hat{\\beta}$'),
  TeX('$\\hat{\\gamma}$'),
  TeX('$Z$'),
  #TeX('$\\hat{\\gamma}_{SS}$'),
  #TeX('$\\hat{\\gamma}_{MAF}$'),
  TeX('Method 1'),
  #TeX('Method1$_{SS}$'),
  #TeX('Method1$_{MAF}$'),
  TeX('Method 2')
  #TeX('Method2$_{SS}$'),
  #TeX('Method2$_{MAF}$')
)

mall[,metric:=factor(metric,levels=c('beta','r_emp_maf_se','z','emp_shrinkage','ws_emp_shrinkage'),labels=metric_names)]
mall[,variable:=factor(gsub("_"," ",as.character(variable)),levels=gsub("_"," ",levels(variable)))]

if(PLOT){
  pp <- ggplot(mall,aes(x=variable,y=value,col=reference)) + geom_boxplot() +
  facet_wrap( ~ metric,scale="free",ncol=3,labeller = label_parsed) +
  xlab("Scenario") + ylab("Distance") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5),
  strip.background =element_rect(fill="grey95")) +
  scale_colour_manual("Origin",values=c("control"="steelblue1","GWAS10"="firebrick1"),labels=c('Control','GWAS10')) +
  background_grid(major = "xy", minor = "none")
  save_plot(paste(out.file.stub,'pdf',sep="."),pp,base_height=10,base_width=10)
  #save_plot(paste(out.file.stub,'scale','pdf',sep="."),pb,base_height=10,base_width=10)
}
saveRDS(mall,file=paste(out.file.stub,'RDS',sep="."))

## code for running on the Q
if(FALSE){
  SCEN.DIR <- '/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios'
  TOTAL <- 1000
  BATCH <- 100
  OUT_DIR <- '/home/ob219/tmp/test_scen/'

  n <- TOTAL/BATCH

  cmd <- sapply(list.files(path=SCEN.DIR,pattern="*.yml",full.names=TRUE),function(f){
    sapply(1:n,function(i){
      sprintf("Rscript /home/ob219/git/simBasis/R/simFullGWAS.R -s %s -o %s -n %d -p",f,OUT_DIR,BATCH)
    })
  })

  write(cmd,file="/home/ob219/git/simBasis/sh/scen.txt")
}



## BELOW HERE IS EXPERIMENTAL CODE FOR malhalanobis distance stuff.


## I think that this is equivalent to the malhalanobis distance
## here we multiply the loading for a given axis by the variance explained by that

if(FALSE){

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

}
