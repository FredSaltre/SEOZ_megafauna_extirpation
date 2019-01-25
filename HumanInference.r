## estimation des dates d'arrivee par solow;
## principe = faire une vraisemblance ponderee, les ponderations dependant
## de la distance entre point d'interet et observations.
## sa base est exponentielle, le param regle par calculs biais/variance
## et minimisation de l'imse par double noyau
#
## hyp bio du solow = l'age du fossile en x est uniforme entre 0 et T(x)=date d'arrivee en x 
##
##

##solow.fct = estimation des dates d'extinction aux points (xx,yy)
##        en sortie Test= date estimee
##                    biais=biais estime
##		    sd   =ecart-type estime
##sg.fct = calcul des coordonnées (xx,yy) des points d'interet
##         nx = gere le pas de grille de depart en long x lat (grills nx x nx)
##         dmmax =distance min entre deux points
##         on elimine les points qui tombent dans l'eau
##         puis les points qui sont à moins de dmmax d'un autre point
##nbcoeur=nombre de coeurs demandes 
solow.fct <- function(xx=xx,yy=yy,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge,nbcoeur=nbcoeur)
{
  cl <- makeCluster(nbcoeur,outfile="")
  registerDoParallel(cl)
  paramsig <- lm(log(SdAge)~Age)$coef ## pour simu plus loin
  distan <- matrix(0,length(xx),length(Lon))
  for( i in 1:length(xx)) {
  for(j in 1:length(Lon)) {
  distan[i,j] <- gdist(xx[i],yy[i],Lon[j],Lat[j],units="km")
   }}
  distan[is.na(distan)]<- max(distan,na.rm=T)
  iappartenance <- foreach (i = 1:ncol(distan),.combine='c') %dopar%
       { order(distan[,i])[1]}
  stopCluster(cl) 
  estim.fct <- function(yobs,sdobs,distan,pas)
   {
    lmv.fct <- function(tmax,yobs=yobs,sdobs=sdobs,pond=pond)
     {
      if(tmax<0){return(10**6)}
      yobs <- yobs/10000
      sdobs <- sdobs/10000
      lm0 <- pnorm(yobs-tmax,mean=0,sd=sdobs,lower.tail=F,log.p=F)
      lm1 <- pnorm(yobs,mean=0,sd=sdobs,lower.tail=F,log.p=F)
      lm2 <- lm0-lm1
      lm2[lm2<=0] <- 10**(-150)
      lm2[is.na(lm2)] <- 10**(-150)
      lm2[!is.finite(lm2)] <- 10**(-150)
      lm2 <- log(lm2)
      lmv <- sum(pond)*log(tmax)-sum(pond*lm2)
      lmv[is.na(lmv)] <- 10**6
      if(!is.finite(lmv)){lmv <- 10**6}
      return(lmv)
     }
   Test <- NULL
     for(i in 1:nrow(distan))
     {
      pond   <- pas*distan[i,]/max(distan)
      pond0 <- exp(-pond*pond)
      pond0[pond*pond>300] <-0
      uu <- nlm(lmv.fct,max(yobs)/10000,yobs=yobs,sdobs=sdobs,pond=pond0,
         stepmax=0.5,iterlim=1000)
      vv <- nlm(lmv.fct,uu$est,yobs=yobs,sdobs=sdobs,pond=pond0,iterlim=1000)
      Test <- c(Test,vv$est)
     }
     Test <- Test*10**4
     return(Test)
    }

   simu.fct <- function(Test,iappartenance,paramsig)
    {
     ysim <- NULL
     for(i in 1:length(iappartenance))
       {
       Tm <- Test[iappartenance[i]]
       ysim <- c(ysim,Tm*runif(1))
       }
     sdsim <- exp(paramsig[1]+paramsig[2]*ysim)
     ysim <- ysim+rnorm(length(ysim),rep(0,length(ysim)),sdsim)
     return(list(ysim=ysim,sdsim=sdsim))
    }

  estpond.fct <- function(Test,iappartenance,paramsig,distan,pas,estim.fct)
   {
    Tsim <- NULL
    for(i in 1:100)
    {
     cdsim <- simu.fct(Test,iappartenance,paramsig)
     Tsim <- cbind(Tsim,estim.fct(cdsim$ysim,cdsim$sdsim,distan,pas))
    }
    msim <- apply(Tsim,1,mean)
    sdsim <- apply(Tsim,1,sd)
    biais <- msim-Test
    mse <- biais*biais+sdsim*sdsim
    imse <- sum(mse)
    return(list(Test=Test,biais=biais,sd=sdsim,imse=imse,pas=pas))
  }

  Test <- estim.fct(Age,SdAge,distan,30)
   cl <- makeCluster(nbcoeur,outfile="")
  registerDoParallel(cl)
  voir <- foreach(iteration=1:100) %dopar% {
    prop <- iteration
 #   print(c(iteration,prop))
    estpond.fct(Test,iappartenance,paramsig,distan,prop,estim.fct)
   }
  stopCluster(cl)
  iv <- NULL
  for(i in 1:100){iv <- c(iv,voir[[i]]$imse)}
  pas <- voir[[order(iv)[1]]]$pas
  Testfinal <- estim.fct(Age,SdAge,distan,pas)
  result <- estpond.fct(Testfinal,iappartenance,paramsig,distan,pas,estim.fct)
  result$pas <- pas
  result$iv <- iv
  result$voir <- voir
  return(result)
}             ##fin solow.fct



sousgrille.fct <- function(nx,dmmax,nbcoeur=nbcoeur)
 {
  xx1 <- rep(1:nx,nx)/(nx+1)
  yy1 <- sort(xx1)
  xx1 <- -180 + xx1*360
  yy1 <-  -59 +yy1*(74+59)
  jj <- map("world",fill=T,plot=F)
  selmap <- map.where(jj,xx1,yy1)
  xx1 <- xx1[!is.na(selmap)]
  yy1 <- yy1[!is.na(selmap)]
  
  cl <- makeCluster(nbcoeur,outfile="")
  registerDoParallel(cl)

  ## chercher a eliminer les points tres proches
  dd <- matrix(0,length(xx1),length(xx1))
  dd <- foreach(i = 1:length(xx1),.combine='rbind',.packages='Imap') %dopar%  {
   ifor <- NULL
   for(j in 1:length(xx1))  {
    ifor <- c(ifor, gdist(xx1[i],yy1[i],xx1[j],yy1[j],units="km"))
   }
   ifor
   }
  dd[is.na(dd)]<- max(dd,na.rm=T)
  diag(dd) <- max(dd)
  dm <- apply(dd,1,min)
  dm[dm < dmmax] <- dmmax
  dd[dd > dmmax] <- dmmax+10


  ichang <- 1
  while(ichang==1)
   {
    iu <- foreach (i = 1:length(unique(yy1)),.combine='c') %dopar% {
       max(apply(dd[yy1 !=unique(yy1)[i],],2,min)/dm)
        }
   ichang <- 0
   if(min(iu) <= 1.00000001) {
    i1 <- order(iu)[1]
    dd <- dd[yy1 != unique(yy1)[i1],]
    xx1 <- xx1[yy1 != unique(yy1)[i1]]
    yy1 <- yy1[yy1 != unique(yy1)[i1]]
    ichang <- 1
    }
    }

  ichang <- 1
  while(ichang==1)
   {
    iu <- foreach (i = 1:length(xx1),.combine='c') %dopar% {
      max(apply(dd[-i,],2,min)/dm)
      } 
     ichang <- 0
    if(min(iu) <= 1.0000001) {
      i1 <- order(iu)[1]
      xx1 <- xx1[-i1]
      yy1 <- yy1[-i1]
      dd <- dd[-i1,]
      ichang <- 1
    }
  }
stopCluster(cl)
return(list(xx1=xx1,yy1=yy1))
}


Full.dat <- read.csv("Humandataset(mostreliable)_run.csv",
header=T,sep=",",dec=".",na.strings="na")

library(maps)
library(doParallel)
library(foreach)
library(Imap)


selobs <- !is.na(Full.dat$Lon) & !is.na(Full.dat$Lat)
Lon <- Full.dat$Lon[selobs]
Lat <- Full.dat$Lat[selobs]
Age <- Full.dat$Calibrated.age[selobs]
SdAge <- Full.dat$Calibrated.error[selobs]


nx <-200
dmmax <- 500
nbcoeur <- 28
sg <- sousgrille.fct(nx,dmmax,nbcoeur=nbcoeur)
xx1 <- sg$xx1
yy1 <- sg$yy1
est_solow <- solow.fct(xx=xx1,yy=yy1,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge,nbcoeur=nbcoeur)
#print("fin estim")
#q()


##hist(est_solow$Test)
##plot(est_solow$Test,est_solow$sd)
##map("world",ylim=c(3,80))
##points(xx1,yy1,pch=19,col=grey(est_solow$Test/75000),cex=0.5)
##points(Full.dat$Lon,Full.dat$Lat,col=2,cex=0.5,pch=19)

sel <- (est_solow$Test < 65000)
##map("world")
##points(xx1[sel],yy1[sel],pch=19,col=grey(est_solow$Test[sel]/65000),cex=0.5)
##points(Full.dat$Lon,Full.dat$Lat,col=2,cex=0.5,pch=19)

res <- as.data.frame(list(Lon=xx1,Lat=yy1,estim=est_solow$Test,biais=est_solow$biais,sd=est_solow$sd))

write.table(res,"Estimates_HumanSolow.csv",sep=";")


