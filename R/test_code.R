## messing around code to build basis and check it a bit

foo<-lapply(names(all.sims),function(ld.block){
  block <- all.sims[[ld.block]]
  test <- lapply(names(block),function(sname){
    SL<-block[[sname]]
    DT <- data.table(pid=rep(SL$snps,NSIMS),z=as.vector(SL$z),
                    beta=as.vector(SL$beta),
                    pos=as.numeric(sub("chr1:","",SL$snps,fixed=TRUE)),
                    sim=paste0('sim',rep(1:NSIMS,times=rep(length(SL$snps),NSIMS))),
                    n1=SL$N1,
                    n=SL$N1+SL$N0,
                    trait=sname,
                    ld.block=ld.block)
    DT.f <- DT[pid %in% snps[snps$asb==TRUE,]$id, ]
    DT.f[,c('or','p.value'):=list(exp(beta),2 * pnorm(abs(z),lower.tail=FALSE))]
    ## add in maf
    maf.DT <- snps[id %in% DT.f$pid,.(id,MAF)]
    setkey(DT.f,pid)
    setkey(maf.DT,id)
    DT.f <- DT.f[maf.DT]
    setnames(DT.f,'MAF','maf')
    DT.f[,.(pid,or,p.value,trait,n,n1,maf,ld.block,sim)]
  })
  rbindlist(test)
})

## test code for creating basis (to see if works)
test <- rbindlist(foo)
test <- split(test,test$sim)


foobar <- lapply(test,function(basis.DT){
  setkey(basis.DT,'pid')
  shrink.DT<-cupcake::compute_shrinkage_metrics(basis.DT)
  basis.mat.emp <- cupcake::create_ds_matrix(basis.DT,shrink.DT)
  ## need to add control where beta is zero
  basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
  pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
  pc.DT<-data.table(pc.emp$x) %>% cbind(trait=rownames(pc.emp$x),.)

  pc.DT <- melt(pc.DT,id.var='trait')
  pc.DT <- dcast(pc.DT[variable %in% c('PC1','PC2'),],trait~variable)
})

foobar <- rbindlist(foobar)
library(ggplot2)
library(ggrepel)
ggplot(foobar,aes(x=PC1,y=PC2,color=trait,label=trait)) + geom_point() + geom_text_repel()



basis.DT <- test[[1]]
setkey(basis.DT,'pid')
shrink.DT<-cupcake::compute_shrinkage_metrics(basis.DT)
basis.mat.emp <- cupcake::create_ds_matrix(basis.DT,shrink.DT)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)
pc.DT<-data.table(pc.emp$x) %>% cbind(trait=rownames(pc.emp$x),.)

pc.DT <- melt(pc.DT,id.var='trait')
pc.DT <- dcast(pc.DT[variable %in% c('PC1','PC2'),],trait~variable)
ggplot(pc.DT,aes(x=PC1,y=PC2)) + geom_point()

## function to process simulations so have a complete set in long format.

SL <- all.sims[[1]]$GWAS1
NSIMS<-10

DT <- data.table(pid=rep(SL$snps,NSIMS),z=as.vector(SL$z),
                beta=as.vector(SL$beta),
                pos=as.numeric(sub("chr1:","",SL$snps,fixed=TRUE)),
                sim=paste0('sim',rep(1:NSIMS,times=rep(length(SL$snps),NSIMS))),
              n1=SL$N1,
              n=SL$N1+SL$N0)
## should call snps$id snps$pid for consistency !
DT.f <- DT[pid %in% snps[snps$asb==TRUE,]$id, ]
## compute p.vals
DT.f[,c('or','p.value'):=list(exp(beta),2 * pnorm(abs(z),lower.tail=FALSE))]
## add in maf
maf.DT <- snps[id %in% DT.f$pid,.(id,MAF)]
setkey(DT.f,pid)
setkey(maf.DT,id)
DT.f <- DT.f[maf.DT]
setnames(DT.f,'MAF','maf')




library(cupcake)

DT.f[,ld.block:=1]
shrink.DT<-cupcake::compute_shrinkage_metrics(DT.f)

basis.mat.emp <- create_ds_matrix(basis.DT,shrink.DT,SHRINKAGE_METHOD)
## need to add control where beta is zero
basis.mat.emp<-rbind(basis.mat.emp,control=rep(0,ncol(basis.mat.emp)))
pc.emp <- prcomp(basis.mat.emp,center=TRUE,scale=FALSE)

## wakefield abf
N <- SL$N0 + SL$N1
prop <- SL$N1/N
DT.f[,pp:=cupcake::wakefield_pp(pval,MAF,N,prop),by=sim]



## possibly don't prune as this gives us the opportunity to see what happens when we
## upsample the number of variants

## build a basis based on the SNPs

#cupcake is setup to read things from disk that but wants things in long format
# also whether to filter - that is down sample so that we have things the same
# as the
if(AS_BASIS_FILTER){

}
