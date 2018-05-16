## code to create simplified chr1 files for parsing

module load bcftools/1.2
export BCFTOOLS_PLUGINS=/usr/local/Cluster-Apps/bcftools/1.2/plugins/

## create a master list of allele fre
cd /home/ob219/rds/rds-cew54-wallace-share/Data/reference/UK10K/
bcftools +fill-AN-AC chr1.bcf.gz | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\n' | gzip -c > /home/ob219/rds/hpc-work/simBasis/support/chr1.AFs.tab.gz

## create a list of SNP identifiers which are found above 5% in the dataset
library(data.table)

snps <- fread('zcat /home/ob219/rds/hpc-work/simBasis/support/chr1.AFs.tab.gz')
## an - total allele counts
## ac - alternate allele counts
setnames(snps,c('chr','start','ref','alt','an','ac'))
snps[,MAF:=ac/an]
snps[MAF>0.5,MAF:=1-MAF]




#snps[,id:=paste(gsub("chr","",chr),pos,sep=":")][,.(id,pos)]
#snps[,c('jp_ld','hm_ld'):=list(-1,-1)]

## annotate with LD blocks from different methods
library(IRanges)
snps.gr <- with(snps,IRanges(start=start,width=1L))
jp_ld <- fread("/home/ob219/rds/rds-cew54-wallace-share/Data/reference/lddetect/EUR/fourier_ls-chr1.bed")
jp_ld.gr <- with(jp_ld,IRanges(start=start,end=stop))
old_ld <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/0.1cM_regions.b37.bed")[V1=="1",]
old_ld.gr <-  with(old_ld,IRanges(start=V2+1,end=V3))

ol<-as.matrix(findOverlaps(snps.gr,jp_ld.gr))
snps[ol[,1],jp_ld:=ol[,2]]
ol<-as.matrix(findOverlaps(snps.gr,old_ld.gr))
snps[ol[,1],hm_ld:=ol[,2]]

## next add in those that are in the basis

as_basis <- fread("/home/ob219/rds/hpc-work/as_basis/support_tab/as_basis_snps.tab")[chr=="1",]
snps[,asb:=FALSE]
snps[start %in% as_basis$position,asb:=TRUE]
snps[,id:=paste(chr,start,sep=':')]
snps <- snps[,.(start,MAF,jp_ld,hm_ld,asb,id)]
## some are missing why ?
## as_basis[!position %in% snps$start,]

saveRDS(snps,"/home/ob219/rds/hpc-work/simBasis/support/chr1.RDS")
saveRDS(snps[MAF>0.01,],"/home/ob219/rds/hpc-work/simBasis/support/chr1_maf_0.01.RDS")
