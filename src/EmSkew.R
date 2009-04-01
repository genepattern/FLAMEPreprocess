#*****************************************************
#  EM algorithm for Multivariate Skew Normal/T Mixture Models
#  Version: 1.0-4
#  
#
#
#  Code By Kui Wang (kwang@maths.uq.edu.au)
#  Updated on 28 June, 2008
#  
#*****************************************************

#----------------------------------------------
#
#  Interface functions to library emskew.dll
#
#----------------------------------------------



EmSkew<-function(dat,ng,dist,ncov=3,clust=NULL,init=NULL, method=1,itmax=1000,epsilon=0.0001,nkmean=10,nrandom=0,nhclust=0,seed=123456,OS)
{
if(dist>4|dist<1|ncov<1|ncov>5) 
stop("model specified not vailable yet")
dat<-as.matrix(dat)
n<-nrow(dat);kk<-9

set.seed(seed)

# start from initial values
if(!missing(init))
{
obj<-emskewfit2(dat,init$pro,init$mu,init$sigma,init$dof,init$delta,ng,dist,ncov,itmax,epsilon,method)
}
else
{
if(missing(clust)) 
clust<-init.mix(dat,ng,dist,ncov,nkmean,nrandom,nhclust,method)

# start from initial partition
obj<-emskewfit(dat,ng,clust,dist,ncov,itmax,epsilon,method)
}



error <- obj$error

#   error code

#   error = 0, converge!
#   error = 1, not reach convergence after itmax iterations;  
#   error = 2, initial values not found;
#   error = 3, estep fails (singular matrix);


msg<-switch(error,
'1' = paste("did not converge within ", itmax, " iterations, try to increase itmax or epsilon"),
'2' = paste("failed to find a initial values! try more initial partitions"),
'3' = paste("estep fails (singular covariance matrix),try more initial partitions"))

cat('\n-----------------------\n\n')
cat(msg,"\n")
cat('\n-----------------------\n\n')

ret<-NULL


if(!error) {

# 2. check the loglikelihood
lk<-obj$lk
#x11();plot(1:length(lk),lk,xlab="itereations",ylab="loglikelihood")
print(lk)
# 3. print/summarize the results: 

msg<-switch(dist, 
'1'=paste(ng,"- Component Multivariate Normal Mixture Model"),
'2'=paste(ng,"- Component Multivariate t      Mixture Model"),
'3'=paste(ng,"- Component Multivariate Skew Normal Mixture Model"),
'4'=paste(ng,"- Component Multivariate Skew t      Mixture Model"))


cat('\n-----------------------\n\n')
cat(msg,"\n")
cat('\n-----------------------\n\n')
 
print(obj[1:kk])
cat('\n-----------------------\n')
ret<-obj
}


ret
}




EmSkewDA<-function(test, ng, dist, ncov, 
training=list(dat=dat,clust=clust), popPAR = NULL,itmax=500,epsilon=0.0001)
# Arguments:
# 1. test:   data to be classified, matrix;
# 2. ng  :   number of groups;
# 3. dist:   1 -- mvn; 2-- mvt; 3 -- msn; 4 -- mst;
# 4. ncov:   1 -- common covaraince; 2 --- common diagnal covariance; ...
# 5. training:  list of training data (matrix) and cluster label (vector);
# 6. popPAR: list of parameters of each component; 
# 7. itmax: maximum number of iterations.
# note that popPAR is a list including: pro, mu, sigma, dof,  delta 

{
if(dist>4|dist<1|ncov<1|ncov>5) 
stop("model specified not vailable yet")


if(missing(training)) 
{
  if(missing(popPAR)) stop("either 'training' or 'popPAR' shoud be given")
}
else
  popPAR<-emskewda(training$dat,ng,dist,ncov,training$clust,itmax,epsilon)


test<-as.matrix(test)

emskewpred(test,dist,ng,popPAR$pro,popPAR$mu,popPAR$sigma,popPAR$dof,popPAR$delta)
}



EmSkewMOD<-function(dist,mu,sigma,step,delta,dof)
{
if(dist>4|dist<1) 
stop("model specified not available yet")

mu<-as.matrix(mu)
p<-nrow(mu)
g<-ncol(mu)

if(dist==3|dist==1) dof<-rep(0,g)
emskewmod(dist,mu,sigma,dof,delta=delta, p,g,step)
}
# calculate the Scale free Weighted (mahalonobis distance) Ratio (SWR)and UWR



getSWR<-function(dat,ng,sigma, clust, tau)
{

intra <- intradist(dat,ng,sigma, clust, tau) 
inter <- interdist(dat,ng,sigma, clust, tau) 

list(SWR=sqrt(intra$OverallIntraDist1/inter$OverallInterDist1),
UWR=sqrt(intra$OverallIntraDist2/inter$OverallInterDist2))

}


#calculate the ICL criteria

getICL<-function(x,ng,dist,ncov,pro,mu,sigma,dof,delta,clust)
{
x<-as.matrix(x)
p=ncol(x)
n=nrow(x)


den<-(as.matrix((ddmix(x, ng, dist, pro, mu, sigma, dof, delta))))

loglik <- 0

for(h in 1:ng)
  loglik=loglik+sum(ifelse(clust==h,log(den[,h]) ,0) )

# number of parameters

nc<-switch(ncov,
    '1' = (ng-1)+ng*p+p*(1+p)/2,  #common covariance
    '2' = (ng-1)+ng*p+p,          #common diagonal covariance
    '3' = (ng-1)+ng*(p+p*(1+p)/2),#general covariance
    '4' = (ng-1)+ng*(p+p),        #general diagonal covariance
    '5' = (ng-1)+ng*(p+1)  )      #sigma(h)*I_p

nc <- switch(dist,
"mvn" = nc,
"mvt" = nc + ng,
"msn" = nc + ng*p,
"mst" = nc + ng*p + ng)

ICL = loglik - nc/2*log(n)

list(ICL=ICL)
}

# mahalonobis distance

mahalonobis<-function(p, ng, mu, sigma) 
{
load.emskew()

#obj<-.C('mahalonobis_',
#as.integer(p),as.integer(ng),as.double(mu),as.double(sigma), dist = double(ng*ng), error = integer(1)) 

#if(obj$error !=0) 
obj<-.C('mahalonobis2_',
as.integer(p),as.integer(ng),as.double(mu),as.double(sigma), dist = double(ng*ng), error = integer(1)) 

if(obj$error ==0) 
matrix(obj$dist, ncol=ng)
}


# end


#---------------------------------------------------------------------------
# fit of Skew Mixture models

# start from initial partition

emskewfit<-function(dat,ng,clust,dist,ncov,itmax,epsilon,method)
{
load.emskew()

if(dist>4|dist<1|ncov<1|ncov>5) 
stop("model specified not vailable yet")

dat<-as.matrix(dat)


obj<-.Fortran('emskew',
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(ng),as.integer(dist),as.integer(ncov),as.integer(itmax),
pro   = double(ng),mu  = double(ng*m),sigm  = double(ng*m*m),
mdof   = double(ng),mdelta= double(ng*m),clust = as.integer(clust),method=as.integer(method),
tao   = double(ng*n),loglik= double(1),aic= double(1),bic= double(1),lk= double(itmax),
double(ng*n),double(ng*n),double(ng*n),double(ng*n),error = integer(1),as.double(epsilon))
lk<-obj$lk;lk<-lk[lk!=0]



list(error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(m,ng)),
sigma=array(obj$sigm,c(m,m,ng)),dof=obj$mdof,delta=array(obj$mdelta,c(m,ng)),
clust=obj$clust,tao=array(obj$tao,c(n,ng)),lk=lk)

}



# start from initial values

emskewfit2<-function(dat,pro,mu,sigma,dof,delta,ng,dist,ncov,itmax,epsilon,method)
{
load.emskew()

if(dist>4|dist<1|ncov<1|ncov>5) 
stop("model specified not vailable yet")

dat<-as.matrix(dat)
obj<-.Fortran('emskew2',
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(ng),as.integer(dist),as.integer(ncov),as.integer(itmax),
pro   = as.double(pro),mu  = as.double(mu),sigm  = as.double(sigma),
mdof   = as.double(dof),mdelta= as.double(delta),clust = integer(n),method=as.integer(method),
tao   = double(ng*n),loglik= double(1),aic= double(1),bic= double(1),lk= double(itmax),
double(ng*n),double(ng*n),double(ng*n),double(ng*n),error = integer(1),as.double(epsilon))

lk<-obj$lk;lk<-lk[lk!=0]
list(error=obj$error,loglik=obj$loglik,bic=obj$bic,aic=obj$aic,
pro=obj$pro,mu= array(obj$mu,c(m,ng)),
sigma=array(obj$sigm,c(m,m,ng)),dof=obj$mdof,delta=array(obj$mdelta,c(m,ng)),
clust=obj$clust,tao=array(obj$tao,c(n,ng)),lk=lk)

}


# Discrimination Analysis (DA)

emskewda<-function(dat,ng,dist,ncov,clust,itmax,epsilon)
{
load.emskew()

dat<-as.matrix(dat)

obj<-.Fortran('emskewda',
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(ng),as.integer(dist),as.integer(ncov),
clust = as.integer(clust),as.integer(itmax),
pro  = double(ng),mu  = double(ng*m),sigm  = double(ng*m*m),
dof   = double(ng),delta = double(ng*m),
loglik= double(1),lk=double(itmax),
tao   = double(ng*n),double(ng*n),double(ng*n),double(ng*n),double(ng*n),
error = integer(1),as.double(epsilon)) 

lk<-obj$lk;lk<-lk[lk!=0]

list(error=obj$error,loglik=obj$loglik,
pro=obj$pro,mu= array(obj$mu,c(m,ng)),
sigma=array(obj$sigm,c(m,m,ng)),dof=obj$dof,
delta=array(obj$delta,c(m,ng)),lk=lk)
}

# mode of skew mixture distribution

emskewmod<-function(dist,mu,sigma,dof,delta,p,ng,step)
{
load.emskew()

obj<-.Fortran('emskewmod',as.integer(dist),
as.double(mu),as.double(sigma),as.double(dof),
as.double(delta),as.integer(p),as.integer(ng),as.double(step),
emu = double(ng*p),esigm = double(ng*p*p),modpt = double(ng*p),
error = integer(1)) 

list(error=obj$error,emu= array(obj$emu,c(p,ng)),
esigma=array(obj$esigm,c(p,p,ng)),modpts=array(obj$modpt,c(p,ng)))
}


# prediction of skew mixture model

emskewpred<-function(dat,dist,ng,pro,mu,sigma,dof,delta)
{
load.emskew()

dat<-as.matrix(dat)
n<-nrow(dat)
m<-ncol(dat)

obj<-.Fortran('emskewpred',
as.double(dat),as.integer(dist),pro=as.double(pro),mu  =as.double(mu),
sigm  =as.double(sigma),dof=as.double(dof),delta =as.double(delta),
as.integer(n),as.integer(m),as.integer(ng),
tao=double(ng*n),error= integer(1)) 


if(obj$error !=0) stop(paste("error:",obj$error))

tao2clust(array(obj$tao,c(n,ng)))

}


#*****************************************
#initial values
#*****************************************


initEmSkew<-function(dat,ng,clust,dist,ncov,method)
{
load.emskew()

if(dist>4|dist<1|ncov<1|ncov>5)
stop("model specified not vailable yet")

dat<-as.matrix(dat)

obj<-.Fortran('initfit',
as.double(dat),
as.integer((n=nrow(dat))),
as.integer((m=ncol(dat))),
as.integer(ng),as.integer(dist),as.integer(ncov),
as.integer(clust),pro   =double(ng),mu  =double(ng*m),
sigm  =double(ng*m*m),mdof =double(ng),delta =double(ng*m),
loglik= double(1),as.double(array(0,c(n,ng))),
as.double(array(0,c(n,ng))),as.double(array(0,c(n,ng))),
as.double(array(0,c(n,ng))),as.double(array(0,c(n,ng))),
error = integer(1),method=as.integer(method)) 


list(error=obj$error,loglik=obj$loglik,
pro=obj$pro,mu= array(obj$mu,c(m,ng)),
sigma=array(obj$sigm,c(m,m,ng)),
dof=obj$mdof,delta=array(obj$delta,c(m,ng)))

}

init.mix<-function(dat,ng,dist,ncov,nkmean,nrandom,nhclust,method)
{
found<-list()
found$loglik<--Inf

n<-nrow(dat)
clust<-rep(1,n)
mclust<-NULL


if(ng>1){
#-----------------------------------------------
#1. keams starts
if(nkmean>0)
for(i in 1:nkmean)
{
clust<-kmeans(dat,ng,nstart=5)$cluster
if(min(table(clust))<10) next
obj<-initEmSkew(dat,ng,clust,dist,ncov,method)
if(length(obj)!=7 | obj$error!=0) next
if(obj$loglik>found$loglik)
{
found<-obj
mclust<-clust
}
} #end of i, nkmean

#2. random starts
if(nrandom>0)
for(i in 1:nrandom) {
clust<-sample(1:ng,n,replace=T)
if(min(table(clust))<10) next
obj<-initEmSkew(dat,ng,clust,dist,ncov,method)
if(length(obj)!=7 | obj$error!=0) next
if(obj$loglik>found$loglik)
{
found<-obj 
mclust<-clust
}

} #end of nrandom 

#3. nature start
if(0) {
clust<-rep(1:ng,rep(n/ng,ng))
obj<-initEmSkew(dat,ng,clust,dist,ncov,method)
if(length(obj)!=7 | obj$error!=0) next
if(obj$loglik>found$loglik)
found<-obj ;mclust<-clust }


#4. herarchical clustering starts
#methods<-c( "ward", "single", "complete", "average", "mcquitty", "median","cen")
methods<-c("complete")
#
if(nhclust>0)
{
dd <- as.dist((1 - cor(t(dat)))/2)  
#Use correlations between variables ``as distance''

for(j in methods){	
clust<- cutree(hclust(dd, j), k = ng)
if(min(table(clust))<10) next

obj<-try(initEmSkew(dat,ng,clust,dist,ncov,method))
if(length(obj)!=7 | obj$error!=0) next

if(obj$loglik>found$loglik)
{
found<-obj;
mclust<-clust
}

} #end of for j
} #end if
#------------------------------------------------
}
else
{
obj<-try(initEmSkew(dat,ng,clust,dist,ncov,method))
if(length(obj)!=7 | obj$error!=0)  
stop("failed to find a initial values")
found<-obj;mclust<-clust
}

if(is.null(mclust)) stop("failed to find a initial values!")
mclust
}




intradist<-function(dat,ng,sigma, clust, tau) 
{
load.emskew()

dat<-as.matrix(dat)

intraobj<-.C('intradist2_',
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(ng),as.integer(clust),sigm  =as.double(sigma),tau=as.double(tau),
dist1=double(ng+1),dist2 = double(ng+1), error = integer(1)) 

if(intraobj$error==0)
list(dist1 = intraobj$dist1[1:ng],dist2 = intraobj$dist2[1:ng],
OverallIntraDist1=intraobj$dist1[ng+1],OverallIntraDist2=intraobj$dist2[ng+1])
}


interdist<-function(dat,ng,sigma, clust, tau) 
{
load.emskew()

dat<-as.matrix(dat)

interobj<-.C('interdist2_',
as.double(dat),as.integer((n=nrow(dat))),as.integer((m=ncol(dat))),
as.integer(ng),as.integer(clust),sigm  =as.double(sigma),tau=as.double(tau),
dist1=double(ng*ng+1),dist2 = double(ng*ng+1), error = integer(1)) 

if(interobj$error ==0)
list(dist1 = matrix(interobj$dist1[1:(ng*ng)],ncol=ng),dist2 = matrix(interobj$dist2[1:(ng*ng)],ncol=ng), 
OverallInterDist1 = interobj$dist1[(ng*ng)+1],OverallInterDist2 = interobj$dist2[(ng*ng)+1]   )
}




ddmix2 <- function(x, ng, dist, props, mu, sigma, dof, delta)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")


  # multiple component mixture

  dens <- matrix(0,ncol=ng,nrow=nrow(x))

  # sum of each normal density value from each component at x  
if(ng>1)
  for (h in 1:ng)
      dens[,h] <- (props[h]*switch(dist,
      msn=  ddmsn(x,mu=mu[,h],sigma=sigma[,,h],           delta=delta[,h]),
      mst=  ddmst(x,mu=mu[,h],sigma=sigma[,,h],dof=dof[h],delta=delta[,h]),
      mvn=  ddmvn(x,mu=mu[,h],sigma=sigma[,,h]                           ),
      mvt=  ddmvt(x,mu=mu[,h],sigma=sigma[,,h],dof=dof[h]               )))
else
      dens[,h] <- (switch(dist,
      msn=  ddmsn(x,mu=mu,sigma=sigma,        delta=delta),
      mst=  ddmst(x,mu=mu,sigma=sigma,dof=dof,delta=delta),
      mvn=  ddmvn(x,mu=mu,sigma=sigma                           ),
      mvt=  ddmvt(x,mu=mu,sigma=sigma,dof=dof               )))


  return(dens)
}   

#end



logddmix <- function(dat, g, dist, mu, sigma, dof, delta)
{
load.emskew()

dat<-as.matrix(dat)
n<- nrow(dat)
p<- ncol(dat)


obj<-.Fortran('ddskew',as.double(dat),as.integer(n),
as.integer(p),as.integer(g),as.integer(dist),
as.double(mu),as.double(sigma),as.double(dof),as.double(delta),
den = double(n*g),error = integer(1))[10:11] 


den<- NULL

if(obj$error) 
 warning("error in calculate the density functin!")
else
den <- matrix(obj$den,ncol=g)


den
}



ddmix <- function(dat, g, dist, prop, mu, sigma, dof, delta)
{

ndist<-switch(dist,
"mvn"=1,
"mvt"=2,
"msn"=3,
"mst"=4,5)

den <- logddmix(dat, g, ndist, mu, sigma, dof, delta)

if(!is.null(den))
 ret<- (t(t(exp(den))*prop))
else
  ret<-NULL

  ret

}


dmvskew.mix <- function(x, g, dist, prop, mu, sigma, dof, delta)
{
  if (!(identical(all.equal(sum(prop), 1), TRUE)))
    stop("Proportions don't sum to one\n")


  # multiple component mixture

ndist<-switch(dist,
"mvn"=1,
"mvt"=2,
"msn"=3,
"mst"=4,5)

den <- logddmix(x, g, ndist, mu, sigma, dof, delta)


ret<-NULL

if(!is.null(den))
ret<- (t(t(exp(den))*prop))
rowSums(ret)

}   







#--------
#-----------------------------------------------------
#  1.  density of multivariate normal distribution
#-----------------------------------------------------




denmvn<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)))
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)


if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)


if((det=det(sigma))==0)
{
warning("singular covariance matrix")
return
}


inv<-solve(sigma)
fdist<-function(x) c(t(x-mu)%*%inv%*%(x-mu))
dist<-apply(y,FUN=fdist,1)

lk<--0.5*(dist+log(2*pi)*p+log(det))
list(den=exp(lk),dist=dist,loglik=sum(lk),lk=lk)
}


ddmvn<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)))
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)


if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)


if((det=det(sigma))==0)
{
warning("singular covariance matrix")
return(NULL)
}


inv<-solve(sigma)
fdist<-function(x) c(t(x-mu)%*%inv%*%(x-mu))
dist<-apply(y,FUN=fdist,1)

lk<--0.5*(dist+log(2*pi)*p+log(det))
exp(lk)
}


#-----------------------------------------------------

#  2.  density of multivariate t distribution

#-----------------------------------------------------


denmvt<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),dof=1)
{

y<-as.matrix(y)
sigma<-as.matrix(sigma)
mu=c(mu)

if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)
#
inv<-as.matrix(solve(sigma))
fdist<-function(x) c(t(x-mu)%*%inv%*%(x-mu))

dist<-apply(y,FUN=fdist,1)

lk<-( (lgamma((dof+p)/2)-log(dof*pi)*(p/2)-lgamma(dof/2)-0.5*log(det))-log(1+dist/dof)*((dof+p)/2))

list(den=exp(lk),dist=dist,loglik=sum(lk),lk=lk)
}


ddmvt<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),dof=1)
{

y<-as.matrix(y)
sigma<-as.matrix(sigma)
mu=c(mu)

if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)
#
inv<-as.matrix(solve(sigma))
fdist<-function(x) c(t(x-mu)%*%inv%*%(x-mu))

dist<-apply(y,FUN=fdist,1)

lk<-( (lgamma((dof+p)/2)-log(dof*pi)*(p/2)-lgamma(dof/2)-0.5*log(det))-log(1+dist/dof)*((dof+p)/2))

exp(lk)
}



#-----------------------------------------------------

# 3. density of Multivariate Skew Normal/t Distribution

#-----------------------------------------------------


denmsn<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),
delta=rep(0,length(mu)))
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)
sigmd<-sigma+delta%*%t(delta)
#

if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

obj<-denmvn(y,mu,sigma=sigmd)

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)

inv<-solve(sigmd)
#  mu and variance of V|Y=y (truncated normal)
mu.v<-c(  t(delta)%*%inv%*%(t(y)-mu))
sg.v<-c(1-t(delta)%*%inv%*%delta)

#
# E(I(V>0)|Y) = pvalue
pvalue<-2*pnorm(mu.v/sqrt(sg.v))

#
lk<-obj$lk+log(pvalue)
list(den=exp(lk),loglik=sum(lk),dist=obj$dist,lk=lk,mu.v=mu.v,sg.v=sg.v)
}


ddmsn<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),
delta=rep(0,length(mu)))
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)
sigmd<-sigma+delta%*%t(delta)
#

if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

obj<-denmvn(y,mu,sigma=sigmd)

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)

inv<-solve(sigmd)
#  mu and variance of V|Y=y (truncated normal)
mu.v<-c(  t(delta)%*%inv%*%(t(y)-mu))
sg.v<-c(1-t(delta)%*%inv%*%delta)

#
# E(I(V>0)|Y) = pvalue
pvalue<-2*pnorm(mu.v/sqrt(sg.v))

#
lk<-obj$lk+log(pvalue)
exp(lk)
}


#-----------------------------------------------------

# 4. density of Multivariate Skew Normal/t Distribution

#-----------------------------------------------------


inverse2<-function(sigma)
{
load.emskew()

sigma<-as.matrix(sigma)
if((p=ncol(sigma))!=nrow(sigma))
stop("check sigma!")

obj<-.Fortran('inverse2',
as.double(sigma),inv=double(p*p),det=double(1),
as.integer(p),error = integer(1),integer(1),integer(p)) 

#if(obj$error !=0) return(NULL)

list(error=obj$error,inv = matrix(obj$inv,ncol=p),det = obj$det) 
}


inverse3<-function(sigma)
{
load.emskew()

sigma<-as.matrix(sigma)
if((p=ncol(sigma))!=nrow(sigma))
stop("check sigma!")

obj<-.Fortran('inverse3',
as.double(sigma),inv=double(p*p),det=double(1),
as.integer(p),error = integer(1)) 

if(obj$error !=0) return(NULL)

list(error=obj$error,inv = matrix(obj$inv,ncol=p),det = obj$det) 
}




denmst<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),
delta=rep(0,length(mu)),dof=1)
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)
mu=c(mu)

sigmd<-sigma+delta%*%t(delta)
if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

inv<-solve(sigmd)

obj<-denmvt(y,mu,sigma=sigmd,dof)

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)

#  mu and variance of V|Y=y (truncated normal)
mu.v<-c(  t(delta)%*%inv%*%(t(y)-mu))
sg.v<-c(1-t(delta)%*%inv%*%delta)*(dof+obj$dist)/(dof+p)

#
# E(I(V>0)|Y) = pvalue
pvalue<-2*pt(mu.v/sqrt(sg.v),df=dof+p)

#
lk<-obj$lk+log(pvalue)
list(den=exp(lk),dist=obj$dist,loglik=sum(lk),lk=lk)
}

ddmst<-function(y,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),
delta=rep(0,length(mu)),dof=1)
{
y<-as.matrix(y)
sigma<-as.matrix(sigma)
mu=c(mu)

sigmd<-sigma+delta%*%t(delta)
if((det=det(sigma))==0) 
{
warning("singular covariance matrix")
return(NULL)
}

inv<-solve(sigmd)

obj<-denmvt(y,mu,sigma=sigmd,dof)

if((p=ncol(sigma))>1 && (ncol(y)==1)) y<-t(y)

#  mu and variance of V|Y=y (truncated normal)
mu.v<-c(  t(delta)%*%inv%*%(t(y)-mu))
sg.v<-c(1-t(delta)%*%inv%*%delta)*(dof+obj$dist)/(dof+p)

#
# E(I(V>0)|Y) = pvalue
pvalue<-2*pt(mu.v/sqrt(sg.v),df=dof+p)

#
lk<-obj$lk+log(pvalue)
exp(lk)
}




#-------------------------------------------------------------------
# Multivariate skew (mvn,mvt,mst,msn) mixture - density values
# 
# Parameters
#
# x - points to compute density at 
# mu - matrix of means
# Sigma - matrix of covariance matrices 
# props - vector of mixing proportions 
# delta -  skew parameters
# dof -- degree of freedom
# dist -- distribution,  'mvt' (2) ,'mvn'(1) ,'mst'(4) ,'msn'(3)
#
# Returns
#
# Density values from the normal mixture (at x)
#-------------------------------------------------------------------

dmvskew.mix2 <- function(x, ng, dist, props, mu, sigma, dof, delta)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one\n")


  # multiple component mixture

  dens <- 0

  # sum of each normal density value from each component at x  
if(ng>1)
  for (h in 1:ng)
      dens <- (dens + props[h]*switch(dist,
      msn=  ddmsn(x,mu=mu[,h],sigma=sigma[,,h],           delta=delta[,h]),
      mst=  ddmst(x,mu=mu[,h],sigma=sigma[,,h],dof=dof[h],delta=delta[,h]),
      mvn=  ddmvn(x,mu=mu[,h],sigma=sigma[,,h]                           ),
      mvt=  ddmvt(x,mu=mu[,h],sigma=sigma[,,h],dof=dof[h]               )))
else
      dens <- (dens + switch(dist,
      msn=  ddmsn(x,mu=mu,sigma=sigma,        delta=delta),
      mst=  ddmst(x,mu=mu,sigma=sigma,dof=dof,delta=delta),
      mvn=  ddmvn(x,mu=mu,sigma=sigma                           ),
      mvt=  ddmvt(x,mu=mu,sigma=sigma,dof=dof               )))


  return(dens)
}   


#loads emskew library
load.emskew <- function()
{
    isWindows <- Sys.info()[["sysname"]]=="Windows"
    if(isWindows)
    {
        if(!is.loaded('emskew.dll'))
            dyn.load('emskew.dll')
    }
    else
    {
        if(!is.loaded('emskew.so'))
            dyn.load('emskew.so')
    }
}



#---------------------------------------

# random number

#---------------------------------------

# based on R

rdmvn<-function (n, mu = rep(0, nrow(sigma)), 
sigma = diag(length(mu)), 
method = c("svd", "chol")) 
{
sigma<-as.matrix(sigma)
    
if (nrow(sigma) != ncol(sigma)) {
        stop("sigma must be a square matrix")
    }
    if (length(mu) != nrow(sigma)) {
        stop("mu and sigma have non-conforming size")
    }
    method <- match.arg(method)
    if (method == "svd") {
        ev <- eigen(sigma, sym = TRUE)$values
        if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1]))) {
            warning("sigma is numerically not positive definite")
        }
        sigsvd <- svd(sigma)
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    }
    if (method == "chol") {
        retval <- chol(sigma, pivot = T)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mu, "+")
    retval
}

#random number of multivariate t

rdmvt<-function(n,mu = rep(0, nrow(sigma)), sigma = diag(length(mu)),dof=1)
{
sigma<-as.matrix(sigma)
u<-rgamma(n,dof/2,dof/2)
t(t(rdmvn(n,sigma=sigma)/sqrt(u))+mu)
}


rdmsn<-function(n,mu = rep(0, nrow(sigma)), 
sigma = diag(length(mu)),delta=rep(0,length(mu)))
{

x<-rdmvn(n,mu,sigma)
z<-abs(rnorm(n))
as.matrix(z%*%t(delta)+x)

}



rdmst<-function(n,mu = rep(0, nrow(sigma)), 
sigma = diag(length(mu)),dof=10,
delta = rep(0,length(mu)))
{

sigma<-as.matrix(sigma)
u<-rgamma(n,dof/2,dof/2)
x<-t(t(rdmvn(n,sigma=sigma)/sqrt(u))+mu)

#x<-rdmvt(n,mu=mu,sigma=sigma,dof=dof)
z<-abs(rnorm(n)/sqrt(u))
as.matrix(z%*%t(delta)+x)
}


rmvskew.mix<-function(n,p,ng,dist,
mu,sigma,dof,delta)
{

if(length(c(n))!=ng) stop("n should be a vector")

dat<-array(0,c(10,p))
mvrand<-function(n,dist,mu,sigma,dof,delta)
{
switch(dist,
mvn=rdmvn(n,mu,sigma                    ),
mvt=rdmvt(n,mu,sigma,dof=dof            ),
msn=rdmsn(n,mu,sigma,        delta=delta),
mst=rdmst(n,mu,sigma,dof=dof,delta=delta))
}

if(ng>1)
for(h in 1:ng)
{
if(n[h]>0)
dat<-rbind(dat,mvrand(n[h],dist,mu[,h],sigma[,,h],dof[h],delta[,h]))

}
else
dat<-rbind(dat,mvrand(n,dist,mu,sigma,dof,delta))



dat[-(1:10),]
}



#-------------------------------------------------------------------
# Functions for contour plotting
#-------------------------------------------------------------------

skewmix.contour.2d <- function(dat, props, mu, sigma, dof, delta, clust, dist="mst",
grid=0.1, nrand=1000,levels=seq(10,90,by=10),gridsize, xlim,ylim,xlab,ylab,...)
{
if(ncol(dat)<2) stop("at least 2-D data needed!")
  
dat<-as.matrix(dat)[,1:2]
ng<-length(props)
p<-2


rx<-range(dat[,1])+c(-0.1,0.1)
ry<-range(dat[,2])+c(-0.1,0.1)


if(1)
{
if (missing(xlim))
   xlim<-rx
if (missing(ylim))
   ylim<-ry
}
else
{
maxSigmas <- 4*max(sigma)
if (missing(xlim))
xlim <- c(min(mu[,1]) - maxSigmas, max(mu[,1]) + maxSigmas)
if (missing(ylim))
ylim <- c(min(mu[,2]) - maxSigmas, max(mu[,2]) + maxSigmas)
}


if (missing(gridsize))
{
x <- seq(xlim[1], xlim[2], by=grid) 
y <- seq(ylim[1], ylim[2], by=grid) 
}
else
{
x <- seq(xlim[1], xlim[2], length=gridsize[1])
y <- seq(ylim[1], ylim[2], length=gridsize[2])
}


#----------------------------------
#  get bins
#----------------------------------

nx <- length(x)
ny <- length(y)
xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=TRUE)))
X <- matrix(xoy, nx*ny, 2, byrow=FALSE)


#----------------------------------
#  get density values at grid points  
#----------------------------------

dens<- dmvskew.mix(X, ng, dist, props, mu, sigma, dof, delta)
dens.mat <- matrix(dens,nx,ny)

#--------

n <- round(nrand*props)

k<-which.min(n)

if( sum(n)< nrand ) 
n[k]<-nrand-sum(n[-k])


rand <- rmvskew.mix(n,p,ng,dist,mu,sigma,dof,delta)
rand.den<-dmvskew.mix(rand, ng, dist, props, mu, sigma, dof, delta)
cont <- quantile(rand.den, prob=levels/100)

#---------------------------------------------
#  plot the 2-D contour
#---------------------------------------------
  if (missing(xlab)) xlab <- "x"
  if (missing(ylab)) ylab <- "y"

  x11()
  plot(x = 0, y = 0,type = "n", xlim = rx, ylim = ry, xlab =xlab, ylab =ylab,main=paste(dist,": g = ",ng) )
  points(dat,pch=20,col=clust)
  
  contour(x, y, dens.mat, nlevel=ncont,levels=cont, add=TRUE,drawlabels=FALSE,
         lty=1,col = 4,...)
    
}



#--------------------------
skewmix.contour.3d <- function(dat, props, mu, sigma, dof, delta, clust, dist="mst",
grid=0.1,nrand=10000,levels=seq(10,90,by=10), gridsize, add=FALSE,origin=c(0,0,0),
xlim, ylim, zlim, xlab, ylab, zlab, alphavec, colors)
{

require(rgl)
require(misc3d)

if(ncol(dat)<3) stop("at least 3-D data needed!")

dat<-as.matrix(dat)[,1:3]
ng<-length(props)

if(1)
{
if (missing(xlim))
   xlim<-range(dat[,1])+c(-0.1,0.1)

if (missing(ylim))
   ylim<-range(dat[,2])+c(-0.1,0.1)

if (missing(zlim))
   zlim<-range(dat[,3])+c(-0.1,0.1)

}
else
{
  maxSigmas <- 3.7*max(sigma)

  if (missing(xlim))
    xlim <- c(min(mu[,1]) - maxSigmas, max(mu[,1]) + maxSigmas)
  if (missing(ylim))
    ylim <- c(min(mu[,2]) - maxSigmas, max(mu[,2]) + maxSigmas)
  if (missing(zlim))
    zlim <- c(min(mu[,3]) - maxSigmas, max(mu[,3]) + maxSigmas)

}  


#----------------------------------
#  get bins
#----------------------------------



if(missing(gridsize))
{
x <- seq(xlim[1], xlim[2], by=grid) #length=gridsize[1])
y <- seq(ylim[1], ylim[2], by=grid) #length=gridsize[2])
z <- seq(zlim[1], zlim[2], by=grid) #length=gridsize[3])
}
else
{
x <- seq(xlim[1], xlim[2], length=gridsize[1])
y <- seq(ylim[1], ylim[2], length=gridsize[2])
z <- seq(zlim[1], zlim[2], length=gridsize[3])
}

nx <- length(x)
ny <- length(y)
nz <- length(z)

xoy <- cbind(rep(x,ny), as.vector(matrix(y,nx,ny,byrow=TRUE)))
X <- matrix(xoy, nx*ny, 2, byrow=FALSE)

#----------------------------------
#  get density values at grid points  
#----------------------------------

 
  dens.array <- array(0, dim=c(nx,ny,nz))
  
  for (i in 1:nz)
  {
      
    dens<-dmvskew.mix(cbind(X,z[i]), ng, dist, props, mu, sigma, dof, delta)
    dens.mat <- matrix(dens, nx,ny,byrow=FALSE)
    dens.array[,,i] <- dens.mat
  }
  

rand<-rmvskew.mix(nrand*props,3,ng,dist,mu,sigma,dof,delta)
rand.den<-dmvskew.mix(rand, ng, dist, props, mu, sigma, dof, delta)
cont <- quantile(rand.den, prob=levels/100)


#---------------------------------------------
#  plot the 3-D contour
#---------------------------------------------
  cont<- c(cont)
  nc  <- length(cont)

  if (missing(colors))
    colors <- rev(heat.colors(nc))
  
  if (missing(alphavec))
    alphavec <- seq(0.1,0.5,length=nc)

  if (missing(xlab)) xlab <- "x"
  if (missing(ylab)) ylab <- "y"
  if (missing(zlab)) zlab <- "z"


#----
  objcont<-contour3d(dens.array, level=cont, x, y, z, draw=FALSE)
#----

 
open3d()
plot3d(dat, type = "p", col=clust, add=add,xlab = xlab, ylab = ylab, zlab = zlab)
#decorate3d(xlim, ylim, zlim, xlab = xlab, ylab = ylab, zlab = zlab, 
#box = FALSE, axes = TRUE, main = NULL, sub = NULL, top = TRUE, aspect = FALSE) 

for(h in seq(along = objcont ))
triangles3d(zipTriangles(objcont[[h]]),col=colors[h],front="line",back="line")

}


#  zipTriangles() is from R Package: misc3d.

zipTriangles <- function(tris) 
{
    n <- nrow(tris$v1)
    if (nrow(tris$v2) != n || nrow(tris$v3) != n)
        stop("vertex arrays must have the same number of rows")
    v <- matrix(0, nrow = 3 * n, ncol = 3)
    v[3 * (1 : n) - 2,] <- tris$v1
    v[3 * (1 : n) - 1,] <- tris$v2
    v[3 * (1 : n),] <- tris$v3
    v
}

# end