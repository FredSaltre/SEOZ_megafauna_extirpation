# 1. read files
# 3. GLS with spatial autocorrelation. Assuming gaussian correlation kernel ("corGaus")

rm(list=ls(all=TRUE))
#library(rgdal)
library(fields)
library(nlme)
library(ncf)
library(AICcmodavg)
library(modEvA)

rm(list=ls(all=TRUE))
data <- read.table("Megafauna_(Climate)Bearing_v4(0kaLag_withSD).csv", sep=",", header=T) ### Modifier ici en header=T (ligne suivante aussi)
mdat<-mean(data$Extinction);SCR<-sum((data$Extinction - mdat)**2)#calcul de la somme des residue sur les donnees
output<-matrix(0,29,4);Vargls<-matrix(0,29,1);i<-1;Cand.models <- list( );i

#1
Cand.models[[i]] <-gls.Tx <- gls(Extinction ~ Temp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Tx, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Tx);output[i,3]<-BIC(gls.Tx);output[i,4]<-gls.Tx$logLik;
sdat<-predict(gls.Tx);SCR2<-sum(gls.Tx$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#2
Cand.models[[i]] <-gls.TxP <- gls(Extinction ~ Temp + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxP, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxP);output[i,3]<-BIC(gls.TxP);output[i,4]<-gls.TxP$logLik; 
sdat<-predict(gls.TxP);SCR2<-sum(gls.TxP$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#3
Cand.models[[i]] <-gls.TxPDD <- gls(Extinction ~ Temp + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDD);output[i,3]<-BIC(gls.TxPDD);output[i,4]<-gls.TxPDD$logLik;
sdat<-predict(gls.TxPDD);SCR2<-sum(gls.TxPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#4
Cand.models[[i]] <-gls.TxPDDFD <- gls(Extinction ~ Temp + Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDDFD);output[i,3]<-BIC(gls.TxPDDFD);output[i,4]<-gls.TxPDDFD$logLik;
sdat<-predict(gls.TxPDDFD);SCR2<-sum(gls.TxPDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#5
Cand.models[[i]] <-gls.TxPDDFDGS <- gls(Extinction ~ Temp + Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDDFDGS);output[i,3]<-BIC(gls.TxPDDFDGS);output[i,4]<-gls.TxPDDFDGS$logLik;
sdat<-predict(gls.TxPDDFDGS);SCR2<-sum(gls.TxPDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#6
Cand.models[[i]] <-gls.TxDD <- gls(Extinction ~ Temp + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDD);output[i,3]<-BIC(gls.TxDD);output[i,4]<-gls.TxDD$logLik; 
sdat<-predict(gls.TxDD);SCR2<-sum(gls.TxDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#7
Cand.models[[i]] <-gls.TxDDFD <- gls(Extinction ~ Temp + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDDFD);output[i,3]<-BIC(gls.TxDDFD);output[i,4]<-gls.TxDDFD$logLik; 
sdat<-predict(gls.TxDDFD);SCR2<-sum(gls.TxDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#8
Cand.models[[i]] <-gls.TxDDFDGS <- gls(Extinction ~ Temp + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDDFDGS);output[i,3]<-BIC(gls.TxDDFDGS);output[i,4]<-gls.TxDDFDGS$logLik;
sdat<-predict(gls.TxDDFDGS);SCR2<-sum(gls.TxDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#9
Cand.models[[i]] <-gls.TxFD <- gls(Extinction ~ Temp + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxFD);output[i,3]<-BIC(gls.TxFD);output[i,4]<-gls.TxFD$logLik;
sdat<-predict(gls.TxFD);SCR2<-sum(gls.TxFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#10
Cand.models[[i]] <-gls.TxFDGS <- gls(Extinction ~ Temp + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxFDGS);output[i,3]<-BIC(gls.TxFDGS);output[i,4]<-gls.TxFDGS$logLik;
sdat<-predict(gls.TxFDGS);SCR2<-sum(gls.TxFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#11
Cand.models[[i]] <-gls.TxPFD <- gls(Extinction ~ Temp + Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPFD);output[i,3]<-BIC(gls.TxPFD);output[i,4]<-gls.TxPFD$logLik;
sdat<-predict(gls.TxPFD);SCR2<-sum(gls.TxPFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#12
Cand.models[[i]] <-gls.TxPFDGS <- gls(Extinction ~ Temp + Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPFDGS);output[i,3]<-BIC(gls.TxPFDGS);output[i,4]<-gls.TxPFDGS$logLik;
sdat<-predict(gls.TxPFDGS);SCR2<-sum(gls.TxPFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#13
Cand.models[[i]] <-gls.TxPGS <- gls(Extinction ~ Temp + Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPGS);output[i,3]<-BIC(gls.TxPGS);output[i,4]<-gls.TxPGS$logLik; 
sdat<-predict(gls.TxPGS);SCR2<-sum(gls.TxPGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#14
Cand.models[[i]] <-gls.P <- gls(Extinction ~ Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.P, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.P);output[i,3]<-BIC(gls.P);output[i,4]<-gls.P$logLik;
sdat<-predict(gls.P);SCR2<-sum(gls.P$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#15
Cand.models[[i]] <-gls.PDD <- gls(Extinction ~ Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDD);output[i,3]<-BIC(gls.PDD);output[i,4]<-gls.PDD$logLik;
sdat<-predict(gls.PDD);SCR2<-sum(gls.PDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#16
Cand.models[[i]] <-gls.PDDFD <- gls(Extinction ~ Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL) ;output[i,2]<-AIC(gls.PDDFD);output[i,3]<-BIC(gls.PDDFD);output[i,4]<-gls.PDDFD$logLik;
sdat<-predict(gls.PDDFD);SCR2<-sum(gls.PDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#17
Cand.models[[i]] <-gls.PDDFDGS <- gls(Extinction ~ Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDDFDGS);output[i,3]<-BIC(gls.PDDFDGS);output[i,4]<-gls.PDDFDGS$logLik;
sdat<-predict(gls.PDDFDGS);SCR2<-sum(gls.PDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#18
Cand.models[[i]] <-gls.PDDGS <- gls(Extinction ~ Prec + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDDGS);output[i,3]<-BIC(gls.PDDGS);output[i,4]<-gls.PDDGS$logLik; 
sdat<-predict(gls.PDDGS);SCR2<-sum(gls.PDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#19
Cand.models[[i]] <-gls.PFD <- gls(Extinction ~ Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PFD);output[i,3]<-BIC(gls.PFD);output[i,4]<-gls.PFD$logLik;
sdat<-predict(gls.PFD);SCR2<-sum(gls.PFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#20
Cand.models[[i]] <-gls.PFDGS <- gls(Extinction ~ Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PFDGS);output[i,3]<-BIC(gls.PFDGS);output[i,4]<-gls.PFDGS$logLik; 
sdat<-predict(gls.PFDGS);SCR2<-sum(gls.PFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#21
Cand.models[[i]] <-gls.PGS <- gls(Extinction ~ Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PGS);output[i,3]<-BIC(gls.PGS);output[i,4]<-gls.PGS$logLik;
sdat<-predict(gls.PGS);SCR2<-sum(gls.PGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#22
Cand.models[[i]] <-gls.DD <- gls(Extinction ~ EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DD);output[i,3]<-BIC(gls.DD);output[i,4]<-gls.DD$logLik; 
sdat<-predict(gls.DD);SCR2<-sum(gls.DD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#23
Cand.models[[i]] <-gls.DDFD <- gls(Extinction ~ EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDFD);output[i,3]<-BIC(gls.DDFD);output[i,4]<-gls.DDFD$logLik; 
sdat<-predict(gls.DDFD);SCR2<-sum(gls.DDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#24
Cand.models[[i]] <-gls.DDFDGS <- gls(Extinction ~ EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDFDGS);output[i,3]<-BIC(gls.DDFDGS);output[i,4]<-gls.DDFDGS$logLik;
sdat<-predict(gls.DDFDGS);SCR2<-sum(gls.DDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#25
Cand.models[[i]] <-gls.DDGS <- gls(Extinction ~ EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDGS);output[i,3]<-BIC(gls.DDGS);output[i,4]<-gls.DDGS$logLik;
sdat<-predict(gls.DDGS);SCR2<-sum(gls.DDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#26
Cand.models[[i]] <-gls.FD <- gls(Extinction ~ Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.FD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.FD);output[i,3]<-BIC(gls.FD);output[i,4]<-gls.FD$logLik; 
sdat<-predict(gls.FD);SCR2<-sum(gls.FD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#27
Cand.models[[i]] <-gls.FDGS <- gls(Extinction ~ Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.FDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.FDGS);output[i,3]<-BIC(gls.FDGS);output[i,4]<-gls.FDGS$logLik;
sdat<-predict(gls.FDGS);SCR2<-sum(gls.FDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#28
Cand.models[[i]] <-gls.GS <- gls(Extinction ~ Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.GS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.GS);output[i,3]<-BIC(gls.GS);output[i,4]<-gls.GS$logLik;
sdat<-predict(gls.GS);SCR2<-sum(gls.GS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#29
Cand.models[[i]] <-gls.Int <- gls(Extinction ~ 1, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Int, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Int);output[i,3]<-BIC(gls.Int);output[i,4]<-gls.Int$logLik; 
sdat<-predict(gls.Int);SCR2<-sum(gls.Int$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; 


write.table(output,"AICweight_statsBearingGrad(Climate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
stat<-aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = TRUE)
write.table(stat,"ModelClassification_statsBearingGrad(Climate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
write.table(Vargls,"VarianceExplainedBearingGrad(Climate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);

