
solow.fct <- function(xx=xx,yy=yy,Lon=Lon,Lat=Lat,Age=Age,
    SdAge=SdAge,tmax=tmax)
{
etape <- 1
plot(Lon,Lat,pch=19,col=grey((Age-min(Age))/(max(Age)-min(Age))))
paramsig <- lm(log(SdAge)~Age)$coef
plot(Age,log(SdAge))
abline(paramsig,col=2,lwd=2) 
distan <- matrix(0,length(xx),length(Lon))
for( i in 1:length(xx)) {
for(j in 1:length(Lon)) {
distan[i,j] <- gdist(xx[i],yy[i],Lon[j],Lat[j],units="km")
}}
distan[is.na(distan)]<- max(distan,na.rm=T)

iappartenance <- foreach (i = 1:ncol(distan),.combine='c') %dopar%
       { order(distan[,i])[1]}

estim.fct <- function(yobs,sdobs,tmax,distan,pas)
{
 tmax <- tmax/10000
 lmv.fct <- function(tmin,tmax=tmax,yobs=yobs,sdobs=sdobs,pond=pond,trace=0)
  {
  if(tmin > tmax){return(10**6)}
  yobs <- yobs/10000
  sdobs <- sdobs/10000
  lm0 <- pnorm(tmax-yobs,mean=0,sd=sdobs,lower.tail=T,log.p=F)
  lm1 <- pnorm(tmin-yobs,mean=0,sd=sdobs,lower.tail=T,log.p=F)
  lm2 <- lm0-lm1
  lm2[lm2<=0] <- 10**(-150)
  lm2[is.na(lm2)] <- 10**(-150)
  lm2[!is.finite(lm2)] <- 10**(-150)
  lm2 <- log(lm2)
  lm2 <- lm2[pond>0]
  pond <- pond[pond>0]
  lmv <- sum(pond)*log(tmax-tmin)-sum(pond*lm2)
  lmv[is.na(lmv)] <- 10**6
  if(!is.finite(lmv)){lmv <- 10**6}
    return(lmv)
  }                   ### fin lmv
 Test <- NULL
 for(i in 1:nrow(distan))
  {
   trace <- 0
   if(i==100){trace <- 1}
   pond   <- pas*distan[i,]/max(distan)
   pond0 <- exp(-pond*pond)
   pond0[pond*pond>300] <- 0
   #uu <- nlm(lmv.fct,min(yobs)/10000,tmax=tmax,yobs=yobs,
   #         sdobs=sdobs,pond=pond,trace=0,
   #	          stepmax=0.5,iterlim=1000)
   #vv <- list(est=uu$est)
   #vv <- nlm(lmv.fct,uu$est,tmax=tmax,yobs=yobs,
   #       sdobs=sdobs,pond=pond,trace=trace,iterlim=1000)
   vv <- optimize(lmv.fct,c(0.2,7),tmax=tmax,yobs=yobs,
          sdobs=sdobs,pond=pond0)
   vv$est <- vv$min	  
   Test <- c(Test,vv$est)
  }
  Test <- Test*10**4
  return(Test)
}     ### fin estim

simu.fct <- function(Test,iappartenance,paramsig,tmax)
{
ysim <- NULL
for(i in 1:length(iappartenance))
{
  Tm <- Test[iappartenance[i]]
  ysim <- c(ysim,Tm+(tmax-Tm)*runif(1))
}
  sdsim <- exp(paramsig[1]+paramsig[2]*ysim)
  ysim <- ysim+rnorm(length(ysim),rep(0,length(ysim)),sdsim)
  ysim[ysim < 1000] <- 1000
  ysim[ysim > 70000] <- 70000
  return(list(ysim=ysim,sdsim=sdsim))
}

estpond.fct <- function(Test,iappartenance,paramsig,tmax,distan,pas,estim.fct)
{
Tsim <- NULL
for(i in 1:100)
{
cdsim <- simu.fct(Test,iappartenance,paramsig,tmax)
Tsim <- cbind(Tsim,estim.fct(cdsim$ysim,cdsim$sdsim,tmax,distan,pas))
}
msim <- apply(Tsim,1,mean)
sdsim <- apply(Tsim,1,sd)
biais <- msim-Test
mse <- biais*biais+sdsim*sdsim
imse <- sum(mse)
return(list(Test=Test,biais=biais,sd=sdsim,imse=imse))
} ## fin estpond

print(c("etape",etape))
if(etape==1) {

Test <- estim.fct(Age,SdAge,tmax,distan,30)
voir <- foreach(iteration=1:100) %dopar% {
prop <- iteration
print(c(iteration,prop))
estpond.fct(Test,iappartenance,paramsig,tmax,distan,prop,estim.fct)
}


} ## fin etape==1

iv <- NULL
for(i in 1:100){iv <- c(iv,voir[[i]]$imse)}
plot(iv,type="l",xlab="pas",ylab="imse")
pas <- order(iv)[1]

Testfinal <- estim.fct(Age,SdAge,tmax,distan,pas)
result <- estpond.fct(Testfinal,iappartenance,paramsig,tmax,
          distan,pas,estim.fct)
result$pas <- pas
result$iv <- iv
result$voir <- voir
return(result)
}             ##fin solow.fct
########################
########################

library(maps)
library(doParallel)
library(foreach)
library(Imap)

Full.dat <- read.csv("AllMegafauna_Finalcutoff_v2.csv",header=T,sep=",",dec=".",na.strings="na")
Full.coord <- read.csv("GridCoordinate_v2.csv",header=T,sep=",",dec=".",na.strings="na")

xx1 <- Full.coord$Lon
yy1 <- Full.coord$Lat
jj <- map("world","australia",fill=T)
selmap <- map.where(jj,xx1,yy1)
xx1 <- xx1[!is.na(selmap)]
yy1 <- yy1[!is.na(selmap)]


cl <- makeCluster(28,outfile="")
registerDoParallel(cl)

selobs <- !is.na(Full.dat$Lon) & !is.na(Full.dat$Lat)
Lon <- Full.dat$Lon[selobs]
Lat <- Full.dat$Lat[selobs]
Age <- Full.dat$Age[selobs]
SdAge <- Full.dat$SdAge[selobs]
tmax <- 90000
est_solow <- solow.fct(xx=xx1,yy=yy1,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge,
tmax=tmax)

#library(maps)
#hist(est_solow$Test)
#plot(est_solow$Test,est_solow$sd)
#map("world","australia")
#points(xx1,yy1,pch=19,col=grey(est_solow$Test/75000),cex=0.5)
#points(Full.dat$Lon,Full.dat$Lat,col=2,cex=0.5,pch=19)


nx <- 300

xca <- rep(1:nx,nx)/(nx+1)
yca <- sort(xca)
xca <- 107 + xca*50
yca <-  (-46)+yca*40
zval <- rep(0,length(xca))

iu <- foreach(i =1:length(xca),.combine=c) %dopar%
{
print(c(i,length(xca)))
library(Imap)
dv <- NULL
for(j in 1:length(xx1))
{
dv <- c(dv,gdist(xca[i],yca[i],xx1[j],yy1[j],unit="km"))
}
iu <- order(dv)[1]
}
stopCluster(cl)

zval <- est_solow$Test[iu]

#library(maps)
#jj <- map("world","australie",fill=T)
#sela <- map.where(jj,xca,yca)
#zval[is.na(sela)] <- (-100000)
#map("world","australie",ylim=c(35,80))
#zval <- matrix(zval,ncol=nx,byrow=F)
#image(sort(unique(xca)),sort(unique(yca)),zval,add=T,zlim=c(0,70000))

resbrut <- list(lon=xx1,lat=yy1,date=est_solow$Test,date=est_solow$sd)
write.csv(resbrut,"Outputs_AllMegafauna_brut_v3.csv")
resfin <- list(lon=xca,lat=yca,date=zval)

