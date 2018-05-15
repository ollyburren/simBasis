## code to see if we can recapitulate manahbolis distance metric
library(MASS)
library(magrittr)

## here simulate some data from the mvn
N=1000
cr = runif(1,min=-1,max=1)
A = matrix(c(1,cr,cr,1),2)
e<-mvrnorm(n = N,rep(0,2),A)
## compute pca
pc<-prcomp(e,center=TRUE,scale=FALSE)
## compute malhalanobis distance from mean of to multivariate variables
man <- mahalanobis(e,colMeans(e),cov(e))
## compute the same from pca info
man.from.pc <- apply(pc$x,1,function(x) x^2/pc$sdev^2) %>% colSums
## plot one against the other to check that they are the same
plot(cbind(man.from.pc,man))
abline(a=0,b=1,col='red')


## accepts a prcomp object and an additional set of projections
## computes the mahalanobis distance between each loading
mahalanobis_pairwise<-function(pc,proj){
  ## add in zero for mahalanobis distance
  load <- pc$x
  den <- pc$sdev^2
  if(exists(proj))
    load<-rbind(load,proj$x)
  ## pairwise distance function
  #dist<-apply(load,1,function(p){
  #  apply(pc$x,1,function(x) (x-p)^2/pc$sdev^2) %>% colSums %>% sqrt
  #})
  ## equivalent and faster
  ## note we can't use scale if we add in the projections - we assume centred
  ## to standardise divide through by the sdev
  stand <- apply(load,1,function(x) x/pc$sdev) %>% t
  dist2 <- as.matrix(dist(stand))
  #mahalanobis
  #man <- apply(pc$x,1,function(x) x^2/den) %>% colSums
}
