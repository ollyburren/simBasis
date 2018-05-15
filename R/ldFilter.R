## simGWAS helper functions

library(optparse)
library(data.table)
library(magrittr)

TEST<-FALSE
option_list = list(
        make_option(c("-l", "--ldblock"), type="character", default=NULL,
              help="ld block to process", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
              help="output directory", metavar="character"),
        make_option(c("-b", "--ld_source"), type="character", default=NULL,
              help="ld block source to use one of jp_ld or hm_ld",metavar="character")
        )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
}else{
  args <- list(
      ldblock='77',
      out_dir="/home/ob219/rds/hpc-work/simBasis/support/ld_support",
      ld_source='jp_ld'
  )
}

print(args)

library(data.table)


BCF_TOOLS='/usr/local/Cluster-Apps/bcftools/1.2/bin/bcftools'
BCF_FILE='/home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/chr1.bcf.gz'
SNP_SUPPORT='/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS'

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


ld_source <- args$ld_source
snps <- readRDS(SNP_SUPPORT)[get(`ld_source`)==args$ldblock,]
ld.cmd <- getLD(snps,ld_source,bcf_maf=0.01)$comm


ldFilter <- function(ld.cmd,snps.DT,r_thresh=0.8){
  h <- gethap(ld.cmd)
  LD <- cor2(h)
  keep <- which(rownames(LD) %in% snps.DT[asb==TRUE,]$id)
  LDF <- LD[keep,]
  tmp <- apply(LDF,2,function(x) any(x>r_thresh))
  snps.DT[,poss.CV:=id %in% names(tmp[tmp==TRUE])]
}

ofile <- sprintf("%s_%s.RDS",args$ld_source,args$ldblock)

saveRDS(ldFilter(ld.cmd,snps),file=file.path(args$out_dir,ofile))
message("Sucess")
if(FALSE){
  library(data.table)
  SNP_SUPPORT='/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS'
  snps <- readRDS(SNP_SUPPORT)
  cmd_jp <- "Rscript /home/ob219/git/simBasis/R/ldFilter.R -l %s  -o /home/ob219/rds/hpc-work/simBasis/support/ld_support -b jp_ld"
  jp <- sprintf(cmd_jp,unique(snps$jp_ld))
  #write(sprintf(cmd,unique(snps$jp_ld)),"/home/ob219/git/simBasis/sh/ld_jp.txt")
  cmd_hm <- "Rscript /home/ob219/git/simBasis/R/ldFilter.R -l %s  -o /home/ob219/rds/hpc-work/simBasis/support/ld_support -b hm_ld"
  hm <- sprintf(cmd_hm,unique(snps$hm_ld))
  write(c(jp,hm),"/home/ob219/git/simBasis/sh/ld.txt")
  # ~/git/slurmer/qlines_csd3.rb -t 01:00:00 ld.txt
  ## recombine to make support files
  recomb <- function(ld_source){
    f <- list.files(path="/home/ob219/rds/hpc-work/simBasis/support/ld_support",pattern="*.RDS",full.names=TRUE)
    f<-f[grep(ld_source,f)]
    DT <- lapply(f,readRDS) %>% rbindlist
    DT <- DT[order(start),]
    setnames(DT,'poss.CV',paste('possCV',ld_source,sep='_'))
    DT
  }
  jp_ld.DT <- recomb('jp_ld')
  setkey(jp_ld.DT,'id')
  hm_ld.DT <- recomb('hm_ld')
  setkey(hm_ld.DT,'id')
  final <- jp_ld.DT[hm_ld.DT[,.(id,possCV_hm_ld)]][order(start),]
  saveRDS(final,file="/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01_pCV.RDS")
}
