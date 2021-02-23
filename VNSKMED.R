VNSKMED<-function(arquivo,nclust=2,NSI=20,imax=6,vmax=3,MAXITER=10,range_s=10,range_bl=30)
{
#nclust = k - Number of clusters
#NSI = Numbber of initial solutions generated
#vmax = maximum neighborhoods
#MAXITER = maximum iterations without improvement
#range_s =  p = shaking - maximum neighborhoods
#range_bl = l = local search -  maximum neighborhoods 
  
  
  
Fobj<-function(medoids,D,nclust)
{ 
n<-nrow(D)
Y<-1:n
X<-setdiff(Y,medoids)
Y[medoids]<-1:nclust
Fm<-rep(0,nclust)
Gm<-replicate(nclust,list(NULL))
Y[X]<-rowMins(D[X,medoids])
for(j in 1:nclust) 
  {Gj<-which(Y==j)
   Fm[j]<-sum(D[medoids[j],Gj])
   Gm[[j]]<-Gj
  }
return(list(Me=medoids,Fmed=Fm,Fobj=sum(Fm)/n,C=Gm))
}  
  
Build<-function(NSI,nclust,D)
{ n<-row(D)
  Me<-t(replicate(NSI,sample(n,nclust)))
  zx=apply(Me,1,function(x) ifelse(length(unique(x))==nclust,1,0))
  zx=which(zx==1)
  Me<-Me[zx,]
  XO<-t(apply(Me,1,function(medoids) Update_medoids(medoids,D,nclust)))
  ixMin<-which.min(sapply(XO,function(x) unlist(x$Fobj)))
  Me=Me[ixMin,]
  C=XO[[ixMin]]$C
  Fobj=XO[[ixMin]]$Fobj
  Fmed=XO[[ixMin]]$Fmed
  return(list(Me=Me,C=C,Fobj=Fobj,Fmed=Fmed))
} 


Nk1<-function(x,range_s,nclust)
{
  medoids<-x$Me
  nk<-sapply(x$C,length)
  Near_medoids<-rep(list(NULL),nclust)
  for(i in 1:nclust)
    {nsec<-min(nk[i],range_s)
     Near_medoids[[i]]<-NG[medoids[i],1:nsec]
    }
  return(Near_medoids)
  
}

  
  
Shaking<-function(x,D,k,range_s,nclust)
{
Nk<-Nk1(x,range_s,nclust)
cluster_select<-sample(nclust,min(k,nclust))
feasible<-FALSE
while(feasible==FALSE)
 { medoids_shake<-x$Me
  for(i in cluster_select) 
    {nrc<-setdiff(1:nclust,i)
     if (nclust>1) {j<-sample(nrc,1)} else {j=nrc}
     medoids_shake[i]<-sample(Nk[[j]],1)
  }
 feasible<-ifelse(length(unique(medoids_shake))==nclust,TRUE,FALSE)
}

xshaking<-Update_medoids(medoids_shake,D,nclust)
return(xshaking)

}  
  
  
Update_medoids<-function(x,D,nclust)
  {
    x<-Fobj(x,D,length(x))

    n<-nrow(D)
    Me<-rep(0,nclust)
    Fobj<-0
    Fmed<-rep(0,nclust)
    C=x$C
    for(i in 1:nclust)
    {
      if (length(C[[i]])>1) 
      {z=cluster::pam(D[C[[i]],C[[i]]],1,diss=TRUE)
      Me[i]=as.numeric(z$medoids)
      Fmed[i]=z$objective[2]*length(C[[i]])
      } else {Me[i]<-C[[i]]} 
      
    }
    fobj<-sum(Fmed)/n
    return(list(Me=Me,Fmed=Fmed,Fobj=fobj,C=C))
    
  }
  
  
  
  
  
Lsearch<-function(x1,D,range_bl,nclust)  
{
medoids<-x1$Me
fobjx1<-x1$Fobj
Nk<-Nk1(x1,range_bl,nclust)
L<-unlist(Nk)
reduziu=TRUE
fbest<-Inf
xbest<-NULL
while (reduziu==TRUE)
  {
  
  reduziu<-FALSE
  for(i in 1:nclust)
   {mi<-matrix(rep(medoids[-i],length(L)),ncol=nclust-1,byrow=TRUE)
    mbl<-cbind(L,mi) 
    ix<-which(apply(mbl,1,function(x) all(duplicated(x)==FALSE))==TRUE)
    mbl<-mbl[ix,]
    if (is.matrix(mbl)==FALSE) {mbl<-t(as.matrix(mbl))}
    vx=apply(mbl,1,function(x) Fobj(x,D,nclust))
    f=which.min(sapply(vx,function(y) y$Fobj))
    if (vx[[f]]$Fobj<fbest)
        { fbest=vx[[f]]$Fobj
          xbest=vx[[f]]
          xbest<-Update_medoids(xbest$Me,D,nclust)
          fbest<-xbest$Fobj
        }
   }
   if (fbest<fobjx1)
     {x1<-xbest
      medoids<-x1$Me
      fobjx1<-x1$Fobj
      Nk<-Nk1(x1,range_bl,nclust)
      L<-unlist(Nk)
      reduziu=TRUE
     
     }
 }  
return(x1)
} 
  
  
  
VNS<-function(NSI,D,vmax,range_s,range_bl,MAXITER,nclust)
{ niter.red<-0
  x<-Build(NSI,nclust,D)
  nclust<-length(x$Me)
  while (niter.red<MAXITER)
    {Ng<-1
     niter.red<-niter.red+1
     while (Ng<vmax)
      { x1<-Shaking(x,D,Ng,range_s*Ng,nclust)
        x2<-Lsearch(x1,D,range_bl*Ng,nclust)
       
        if (x2$Fobj<x$Fobj)
          {x<-x2
           Ng<-1
           niter.red<-0
          }
        else {Ng<-Ng+1}
      }
      
    }
    return(x)  
}  
  
  
###############################################################  
########## Leitura dos Dados ################################  
setwd("G:\\Trabs\\VNS_MEDOIDS\\Dados\\")   
B<-read.table(arquivo,sep=";",dec=",") 
D<-as.matrix(dist(matrix(scale(B),ncol=ncol(B))))
  

library(cluster)
library(Rfast)
options(digits=7)
library(parallel)
nucleos<-detectCores(logical=TRUE)
clust<-makeCluster(nucleos)
clusterEvalQ(clust, library(cluster))
clusterEvalQ(clust, library(Rfast))

#############################################################
###Medoids Iniciais$########################################

tempo<-proc.time()

i<-1
Fbest<-Inf
NG<-rowOrder(D)
NG<-NG[,-1]

NSI<-rep(NSI,imax)
s=clusterApplyLB(cl=clust,NSI,function(y) VNS(y,D,vmax,range_s,range_bl,MAXITER,nclust))
ixmin<-which.min(sapply(s,function(x) x$Fobj))
xbest<-s[[ixmin]]
stopCluster(clust)
tempo<-(proc.time()-tempo)[3]
Fobj<-xbest$Fobj
Med<-xbest$Me
Gks<-xbest$C
clustering<-rep(0,length(unlist(Gks)))
for(kx in 1:nclust) {clustering[Gks[[kx]]]<-kx}
cat("VNS Solution ",Fobj," CPU Time = ",tempo,"\n")
return(list(Fobj=Fobj,Med=Med,clustering=clustering,tempo=tempo))
return(s)

}  




