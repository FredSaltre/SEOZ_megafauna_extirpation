


solow_euclide.fct <- function(xx=xx,yy=yy,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge)
{
etape <- 1
plot(Lon,Lat,pch=19,col=grey((Age-min(Age))/(max(Age)-min(Age))))
paramsig <- lm(log(SdAge)~Age)$coef
plot(Age,log(SdAge))
abline(paramsig,col=2,lwd=2) 


distan <- outer(xx,Lon,"-")**2 +outer(yy,Lat,"-")**2
distan <- sqrt(distan)


iappartenance <- foreach (i = 1:ncol(distan),.combine='c') %dopar%
       { order(distan[,i])[1]}

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
   pond <- exp(-pond*pond)
   uu <- nlm(lmv.fct,max(yobs)/10000,yobs=yobs,sdobs=sdobs,pond=pond,
         stepmax=0.5,iterlim=1000)
   vv <- nlm(lmv.fct,uu$est,yobs=yobs,sdobs=sdobs,pond=pond,iterlim=1000)
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
vu <- range(Test)
vu <- vu[1]+(vu[2]-vu[1])*(0:5)/5
leg <- as.character(round(vu,3)*10000)
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
return(list(Test=Test,biais=biais,sd=sdsim,imse=imse))
}


if(etape==1) {

Test <- estim.fct(Age,SdAge,distan,20)
voir <- foreach(iteration=1:100) %dopar% {
prop <- iteration
print(c(iteration,prop))
estpond.fct(Test,iappartenance,paramsig,distan,prop,estim.fct)
}


} ## fin etape==1

iv <- NULL
for(i in 1:100){iv <- c(iv,voir[[i]]$imse)}
plot(iv,type="l",xlab="pas",ylab="imse")
pas <- order(iv)[1]

Testfinal <- estim.fct(Age,SdAge,distan,pas)
result <- estpond.fct(Testfinal,iappartenance,paramsig,distan,pas,estim.fct)
result$pas <- pas
result$iv <- iv
result$voir <- voir
return(result)
}             ##fin solow.fct
########################
########################






library(doParallel)
library(foreach)

source("humgen.r")
mat <- humgen.fct(500,1,1)

nx <- 30
xx1 <- rep(1:nx,nx)/(nx+1)
yy1 <- sort(xx1)

Lon <- mat[,1]
Lat <- mat[,2]
Age <- mat[,3]
SdAge <- mat[,4]


cl <- makeCluster(18,outfile="")
registerDoParallel(cl)



est_solow <- solow_euclide.fct(xx=xx1,yy=yy1,Lon=Lon,Lat=Lat,Age=Age,SdAge=SdAge)

stopCluster(cl)

res.est<-cbind(xx1,yy1,est_solow$Test);res.bias<-cbind(xx1,yy1,est_solow$bias);
res.sd<-cbind(xx1,yy1,est_solow$sd);res.imse<-est_solow$imse;

write.csv(res.est,file="FinalEstimates.csv",row.names = FALSE)
write.csv(res.bias,file="Spatialbias.csv",row.names = FALSE)
write.csv(res.sd,file="Spatialvariance.csv",row.names = FALSE)
write.csv(res.imse,file="Imse.csv",row.names = FALSE)
write.csv(mat,file="Simulated_dataInitial.csv",row.names = FALSE)