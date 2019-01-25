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
data <- read.table("Megafauna_(ClimateHuman)Bearing_v4(0kaLag_withSD).csv", sep=",", header=T) ### Modifier ici en header=T (ligne suivante aussi)
mdat<-mean(data$Extinction);SCR<-sum((data$Extinction - mdat)**2)#calcul de la somme des residue sur les donnees
output<-matrix(0,60,4);Vargls<-matrix(0,60,1);i<-1;Cand.models <- list( );i

#Mod 1
Cand.models[[i]] <- gls.Tn <- gls(Extinction ~ Arrival, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Tn, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Tn);output[i,3]<-BIC(gls.Tn);output[i,4]<-gls.Tn$logLik;
sdat<-predict(gls.Tn);SCR2<-sum(gls.Tn$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 2
Cand.models[[i]] <-gls.TnTx <- gls(Extinction ~ Arrival + Temp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTx, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTx);output[i,3]<-BIC(gls.TnTx);output[i,4]<-gls.TnTx$logLik;
sdat<-predict(gls.TnTx);SCR2<-sum(gls.TnTx$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 3
Cand.models[[i]] <-gls.TnTxP <- gls(Extinction ~ Arrival + Temp + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxP, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxP);output[i,3]<-BIC(gls.TnTxP);output[i,4]<-gls.TnTxP$logLik;
sdat<-predict(gls.TnTxP);SCR2<-sum(gls.TnTxP$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 4
Cand.models[[i]] <-gls.TnTxPDD <- gls(Extinction ~ Arrival + Temp + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPDD);output[i,3]<-BIC(gls.TnTxPDD);output[i,4]<-gls.TnTxPDD$logLik;
sdat<-predict(gls.TnTxPDD);SCR2<-sum(gls.TnTxPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 5
Cand.models[[i]] <-gls.TnTxPDDFD <- gls(Extinction ~ Arrival + Temp + Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPDDFD);output[i,3]<-BIC(gls.TnTxPDDFD);output[i,4]<-gls.TnTxPDDFD$logLik;
sdat<-predict(gls.TnTxPDDFD);SCR2<-sum(gls.TnTxPDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 6
Cand.models[[i]] <-gls.TnTxPDDFDGS <- gls(Extinction ~ Arrival + Temp + Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPDDFDGS);output[i,3]<-BIC(gls.TnTxPDDFDGS);output[i,4]<-gls.TnTxPDDFDGS$logLik;
sdat<-predict(gls.TnTxPDDFDGS);SCR2<-sum(gls.TnTxPDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 7
Cand.models[[i]] <-gls.TnTxDD <- gls(Extinction ~ Arrival + Temp + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxDD);output[i,3]<-BIC(gls.TnTxDD);output[i,4]<-gls.TnTxDD$logLik;
sdat<-predict(gls.TnTxDD);SCR2<-sum(gls.TnTxDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 8
Cand.models[[i]] <-gls.TnTxDDFD <- gls(Extinction ~ Arrival + Temp + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxDDFD);output[i,3]<-BIC(gls.TnTxDDFD);output[i,4]<-gls.TnTxDDFD$logLik;
sdat<-predict(gls.TnTxDDFD);SCR2<-sum(gls.TnTxDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 9
Cand.models[[i]] <-gls.TnTxDDFDGS <- gls(Extinction ~ Arrival + Temp + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxDDFDGS);output[i,3]<-BIC(gls.TnTxDDFDGS);output[i,4]<-gls.TnTxDDFDGS$logLik; 
sdat<-predict(gls.TnTxDDFDGS);SCR2<-sum(gls.TnTxDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 10
Cand.models[[i]] <-gls.TnTxPDDGS <- gls(Extinction ~ Arrival + Temp + Prec + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPDDGS);output[i,3]<-BIC(gls.TnTxPDDGS);output[i,4]<-gls.TnTxPDDGS$logLik;
sdat<-predict(gls.TnTxPDDGS);SCR2<-sum(gls.TnTxPDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 11
Cand.models[[i]] <-gls.TnTxPFD <- gls(Extinction ~ Arrival + Temp + Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPFD);output[i,3]<-BIC(gls.TnTxPFD);output[i,4]<-gls.TnTxPFD$logLik;
sdat<-predict(gls.TnTxPFD);SCR2<-sum(gls.TnTxPFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 12
Cand.models[[i]] <-gls.TnTxPFDGS <- gls(Extinction ~ Arrival + Temp + Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPFDGS);output[i,3]<-BIC(gls.TnTxPFDGS);output[i,4]<-gls.TnTxPFDGS$logLik;
sdat<-predict(gls.TnTxPFDGS);SCR2<-sum(gls.TnTxPFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 13
Cand.models[[i]] <-gls.TnTxFD <- gls(Extinction ~ Arrival + Temp + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxFD);output[i,3]<-BIC(gls.TnTxFD);output[i,4]<-gls.TnTxFD$logLik;
sdat<-predict(gls.TnTxFD);SCR2<-sum(gls.TnTxFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 14
Cand.models[[i]] <-gls.TnTxFDGS <- gls(Extinction ~ Arrival + Temp + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxFDGS);output[i,3]<-BIC(gls.TnTxFDGS);output[i,4]<-gls.TnTxFDGS$logLik;
sdat<-predict(gls.TnTxFDGS);SCR2<-sum(gls.TnTxFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 15
Cand.models[[i]] <-gls.TnTxGS <- gls(Extinction ~ Arrival + Temp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxGS);output[i,3]<-BIC(gls.TnTxGS);output[i,4]<-gls.TnTxGS$logLik; 
sdat<-predict(gls.TnTxGS);SCR2<-sum(gls.TnTxGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 16
Cand.models[[i]] <-gls.TnTxDDGS <- gls(Extinction ~ Arrival + Temp + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxDDGS);output[i,3]<-BIC(gls.TnTxDDGS);output[i,4]<-gls.TnTxDDGS$logLik;
sdat<-predict(gls.TnTxDDGS);SCR2<-sum(gls.TnTxDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 17
Cand.models[[i]] <-gls.TnTxPGS <- gls(Extinction ~ Arrival + Temp + Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnTxPGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnTxPGS);output[i,3]<-BIC(gls.TnTxPGS);output[i,4]<-gls.TnTxPGS$logLik;
sdat<-predict(gls.TnTxPGS);SCR2<-sum(gls.TnTxPGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 18
Cand.models[[i]] <-gls.TnPFD <- gls(Extinction ~ Arrival + Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPFD);output[i,3]<-BIC(gls.TnPFD);output[i,4]<-gls.TnPFD$logLik;
sdat<-predict(gls.TnPFD);SCR2<-sum(gls.TnPFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 19
Cand.models[[i]] <-gls.TnPFDGS <- gls(Extinction ~ Arrival + Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPFDGS);output[i,3]<-BIC(gls.TnPFDGS);output[i,4]<-gls.TnPFDGS$logLik;
sdat<-predict(gls.TnPFDGS);SCR2<-sum(gls.TnPFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 20
Cand.models[[i]] <-gls.TnDDGS <- gls(Extinction ~ Arrival + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDDGS);output[i,3]<-BIC(gls.TnDDGS);output[i,4]<-gls.TnDDGS$logLik;
sdat<-predict(gls.TnDDGS);SCR2<-sum(gls.TnDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 21
Cand.models[[i]] <-gls.TnPGS <- gls(Extinction ~ Arrival + Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPGS);output[i,3]<-BIC(gls.TnPGS);output[i,4]<-gls.TnPGS$logLik; 
sdat<-predict(gls.TnPGS);SCR2<-sum(gls.TnPGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 22
Cand.models[[i]] <-gls.TnP <- gls(Extinction ~ Arrival + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnP, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnP);output[i,3]<-BIC(gls.TnP);output[i,4]<-gls.TnP$logLik;
sdat<-predict(gls.TnP);SCR2<-sum(gls.TnP$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 23
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 24
Cand.models[[i]] <-gls.TnPDDFD <- gls(Extinction ~ Arrival + Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDDFD);output[i,3]<-BIC(gls.TnPDDFD);output[i,4]<-gls.TnPDDFD$logLik;
sdat<-predict(gls.TnPDDFD);SCR2<-sum(gls.TnPDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 25
Cand.models[[i]] <-gls.TnPDDFDGS <- gls(Extinction ~ Arrival + Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDDFDGS);output[i,3]<-BIC(gls.TnPDDFDGS);output[i,4]<-gls.TnPDDFDGS$logLik;
sdat<-predict(gls.TnPDDFDGS);SCR2<-sum(gls.TnPDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 26
Cand.models[[i]] <-gls.TnDD <- gls(Extinction ~ Arrival + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDD);output[i,3]<-BIC(gls.TnDD);output[i,4]<-gls.TnDD$logLik; 
sdat<-predict(gls.TnDD);SCR2<-sum(gls.TnDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 27
Cand.models[[i]] <-gls.TnDDFD <- gls(Extinction ~ Arrival + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDDFD);output[i,3]<-BIC(gls.TnDDFD);output[i,4]<-gls.TnDDFD$logLik;
sdat<-predict(gls.TnDDFD);SCR2<-sum(gls.TnDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 28
Cand.models[[i]] <-gls.TnDDFDGS <- gls(Extinction ~ Arrival + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDDFDGS);output[i,3]<-BIC(gls.TnDDFDGS);output[i,4]<-gls.TnDDFDGS$logLik;
sdat<-predict(gls.TnDDFDGS);SCR2<-sum(gls.TnDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 29
Cand.models[[i]] <-gls.TnFD <- gls(Extinction ~ Arrival + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnFD);output[i,3]<-BIC(gls.TnFD);output[i,4]<-gls.TnFD$logLik;
sdat<-predict(gls.TnFD);SCR2<-sum(gls.TnFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 30
Cand.models[[i]] <-gls.TnFDGS <- gls(Extinction ~ Arrival + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnFDGS);output[i,3]<-BIC(gls.TnFDGS);output[i,4]<-gls.TnFDGS$logLik; 
sdat<-predict(gls.TnFDGS);SCR2<-sum(gls.TnFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 31
Cand.models[[i]] <-gls.TnGS <- gls(Extinction ~ Arrival + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnGS);output[i,3]<-BIC(gls.TnGS);output[i,4]<-gls.TnGS$logLik;
sdat<-predict(gls.TnGS);SCR2<-sum(gls.TnGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 32
Cand.models[[i]] <-gls.Tx <- gls(Extinction ~ Temp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Tx, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Tx);output[i,3]<-BIC(gls.Tx);output[i,4]<-gls.Tx$logLik;
sdat<-predict(gls.Tx);SCR2<-sum(gls.Tx$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 33
Cand.models[[i]] <-gls.TxP <- gls(Extinction ~ Temp + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxP, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxP);output[i,3]<-BIC(gls.TxP);output[i,4]<-gls.TxP$logLik; 
sdat<-predict(gls.TxP);SCR2<-sum(gls.TxP$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 34
Cand.models[[i]] <-gls.TxPDD <- gls(Extinction ~ Temp + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDD);output[i,3]<-BIC(gls.TxPDD);output[i,4]<-gls.TxPDD$logLik;
sdat<-predict(gls.TxPDD);SCR2<-sum(gls.TxPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 35
Cand.models[[i]] <-gls.TxPDDFD <- gls(Extinction ~ Temp + Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDDFD);output[i,3]<-BIC(gls.TxPDDFD);output[i,4]<-gls.TxPDDFD$logLik;
sdat<-predict(gls.TxPDDFD);SCR2<-sum(gls.TxPDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 36
Cand.models[[i]] <-gls.TxPDDFDGS <- gls(Extinction ~ Temp + Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPDDFDGS);output[i,3]<-BIC(gls.TxPDDFDGS);output[i,4]<-gls.TxPDDFDGS$logLik;
sdat<-predict(gls.TxPDDFDGS);SCR2<-sum(gls.TxPDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 37
Cand.models[[i]] <-gls.TxDD <- gls(Extinction ~ Temp + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDD);output[i,3]<-BIC(gls.TxDD);output[i,4]<-gls.TxDD$logLik; 
sdat<-predict(gls.TxDD);SCR2<-sum(gls.TxDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 38
Cand.models[[i]] <-gls.TxDDFD <- gls(Extinction ~ Temp + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDDFD);output[i,3]<-BIC(gls.TxDDFD);output[i,4]<-gls.TxDDFD$logLik; 
sdat<-predict(gls.TxDDFD);SCR2<-sum(gls.TxDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 39
Cand.models[[i]] <-gls.TxDDFDGS <- gls(Extinction ~ Temp + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxDDFDGS);output[i,3]<-BIC(gls.TxDDFDGS);output[i,4]<-gls.TxDDFDGS$logLik;
sdat<-predict(gls.TxDDFDGS);SCR2<-sum(gls.TxDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 40
Cand.models[[i]] <-gls.TxFD <- gls(Extinction ~ Temp + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxFD);output[i,3]<-BIC(gls.TxFD);output[i,4]<-gls.TxFD$logLik;
sdat<-predict(gls.TxFD);SCR2<-sum(gls.TxFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 41
Cand.models[[i]] <-gls.TxFDGS <- gls(Extinction ~ Temp + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxFDGS);output[i,3]<-BIC(gls.TxFDGS);output[i,4]<-gls.TxFDGS$logLik;
sdat<-predict(gls.TxFDGS);SCR2<-sum(gls.TxFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 42
Cand.models[[i]] <-gls.TxPFD <- gls(Extinction ~ Temp + Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPFD);output[i,3]<-BIC(gls.TxPFD);output[i,4]<-gls.TxPFD$logLik;
sdat<-predict(gls.TxPFD);SCR2<-sum(gls.TxPFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 43
Cand.models[[i]] <-gls.TxPFDGS <- gls(Extinction ~ Temp + Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPFDGS);output[i,3]<-BIC(gls.TxPFDGS);output[i,4]<-gls.TxPFDGS$logLik;
sdat<-predict(gls.TxPFDGS);SCR2<-sum(gls.TxPFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 44
Cand.models[[i]] <-gls.TxPGS <- gls(Extinction ~ Temp + Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TxPGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TxPGS);output[i,3]<-BIC(gls.TxPGS);output[i,4]<-gls.TxPGS$logLik; 
sdat<-predict(gls.TxPGS);SCR2<-sum(gls.TxPGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 45
Cand.models[[i]] <-gls.P <- gls(Extinction ~ Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.P, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.P);output[i,3]<-BIC(gls.P);output[i,4]<-gls.P$logLik;
sdat<-predict(gls.P);SCR2<-sum(gls.P$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 46
Cand.models[[i]] <-gls.PDD <- gls(Extinction ~ Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDD);output[i,3]<-BIC(gls.PDD);output[i,4]<-gls.PDD$logLik;
sdat<-predict(gls.PDD);SCR2<-sum(gls.PDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 47
Cand.models[[i]] <-gls.PDDFD <- gls(Extinction ~ Prec + EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL) ;output[i,2]<-AIC(gls.PDDFD);output[i,3]<-BIC(gls.PDDFD);output[i,4]<-gls.PDDFD$logLik;
sdat<-predict(gls.PDDFD);SCR2<-sum(gls.PDDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 48
Cand.models[[i]] <-gls.PDDFDGS <- gls(Extinction ~ Prec + EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDDFDGS);output[i,3]<-BIC(gls.PDDFDGS);output[i,4]<-gls.PDDFDGS$logLik;
sdat<-predict(gls.PDDFDGS);SCR2<-sum(gls.PDDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 49
Cand.models[[i]] <-gls.PDDGS <- gls(Extinction ~ Prec + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PDDGS);output[i,3]<-BIC(gls.PDDGS);output[i,4]<-gls.PDDGS$logLik; 
sdat<-predict(gls.PDDGS);SCR2<-sum(gls.PDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 50
Cand.models[[i]] <-gls.PFD <- gls(Extinction ~ Prec + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PFD);output[i,3]<-BIC(gls.PFD);output[i,4]<-gls.PFD$logLik;
sdat<-predict(gls.PFD);SCR2<-sum(gls.PFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 51
Cand.models[[i]] <-gls.PFDGS <- gls(Extinction ~ Prec + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PFDGS);output[i,3]<-BIC(gls.PFDGS);output[i,4]<-gls.PFDGS$logLik; 
sdat<-predict(gls.PFDGS);SCR2<-sum(gls.PFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 52
Cand.models[[i]] <-gls.PGS <- gls(Extinction ~ Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.PGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.PGS);output[i,3]<-BIC(gls.PGS);output[i,4]<-gls.PGS$logLik;
sdat<-predict(gls.PGS);SCR2<-sum(gls.PGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 53
Cand.models[[i]] <-gls.DD <- gls(Extinction ~ EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DD);output[i,3]<-BIC(gls.DD);output[i,4]<-gls.DD$logLik; 
sdat<-predict(gls.DD);SCR2<-sum(gls.DD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 54
Cand.models[[i]] <-gls.DDFD <- gls(Extinction ~ EminP + Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDFD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDFD);output[i,3]<-BIC(gls.DDFD);output[i,4]<-gls.DDFD$logLik; 
sdat<-predict(gls.DDFD);SCR2<-sum(gls.DDFD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 55
Cand.models[[i]] <-gls.DDFDGS <- gls(Extinction ~ EminP + Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDFDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDFDGS);output[i,3]<-BIC(gls.DDFDGS);output[i,4]<-gls.DDFDGS$logLik;
sdat<-predict(gls.DDFDGS);SCR2<-sum(gls.DDFDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 56
Cand.models[[i]] <-gls.DDGS <- gls(Extinction ~ EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.DDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.DDGS);output[i,3]<-BIC(gls.DDGS);output[i,4]<-gls.DDGS$logLik;
sdat<-predict(gls.DDGS);SCR2<-sum(gls.DDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 57
Cand.models[[i]] <-gls.FD <- gls(Extinction ~ Npp, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.FD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.FD);output[i,3]<-BIC(gls.FD);output[i,4]<-gls.FD$logLik; 
sdat<-predict(gls.FD);SCR2<-sum(gls.FD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 58
Cand.models[[i]] <-gls.FDGS <- gls(Extinction ~ Npp + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.FDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.FDGS);output[i,3]<-BIC(gls.FDGS);output[i,4]<-gls.FDGS$logLik;
sdat<-predict(gls.FDGS);SCR2<-sum(gls.FDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 59
Cand.models[[i]] <-gls.GS <- gls(Extinction ~ Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.GS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.GS);output[i,3]<-BIC(gls.GS);output[i,4]<-gls.GS$logLik;
sdat<-predict(gls.GS);SCR2<-sum(gls.GS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 60
Cand.models[[i]] <-gls.Int <- gls(Extinction ~ 1, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Int, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Int);output[i,3]<-BIC(gls.Int);output[i,4]<-gls.Int$logLik; 
sdat<-predict(gls.Int);SCR2<-sum(gls.Int$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; 

write.table(output,"AICweight_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
stat<-aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = TRUE)
write.table(stat,"ModelClassification_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
write.table(Vargls,"VarianceExplained_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);


##############################################################################################################
## REFINE USING MODEL WITH DEVIANCE EXPLAINED > 80%
##############################################################################################################
rm(list=ls(all=TRUE))
#library(rgdal)
library(fields)
library(nlme)
library(ncf)
library(AICcmodavg)
library(modEvA)

rm(list=ls(all=TRUE))
data <- read.table("Megafauna_(ClimateHuman)Bearing_v4(0kaLag_withSD).csv", sep=",", header=T) ### Modifier ici en header=T (ligne suivante aussi)
mdat<-mean(data$Extinction);SCR<-sum((data$Extinction - mdat)**2)#calcul de la somme des residue sur les donnees
output<-matrix(0,7,4);Vargls<-matrix(0,7,1);i<-1;Cand.models <- list( );i

#Mod 1
Cand.models[[i]] <- gls.Tn <- gls(Extinction ~ Arrival, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.Tn, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.Tn);output[i,3]<-BIC(gls.Tn);output[i,4]<-gls.Tn$logLik;
sdat<-predict(gls.Tn);SCR2<-sum(gls.Tn$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 2
Cand.models[[i]] <-gls.TnDDGS <- gls(Extinction ~ Arrival + EminP + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDDGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDDGS);output[i,3]<-BIC(gls.TnDDGS);output[i,4]<-gls.TnDDGS$logLik;
sdat<-predict(gls.TnDDGS);SCR2<-sum(gls.TnDDGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 3
Cand.models[[i]] <-gls.TnPGS <- gls(Extinction ~ Arrival + Prec + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPGS);output[i,3]<-BIC(gls.TnPGS);output[i,4]<-gls.TnPGS$logLik; 
sdat<-predict(gls.TnPGS);SCR2<-sum(gls.TnPGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 4
Cand.models[[i]] <-gls.TnP <- gls(Extinction ~ Arrival + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnP, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnP);output[i,3]<-BIC(gls.TnP);output[i,4]<-gls.TnP$logLik;
sdat<-predict(gls.TnP);SCR2<-sum(gls.TnP$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 5
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 6
Cand.models[[i]] <-gls.TnDD <- gls(Extinction ~ Arrival + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnDD);output[i,3]<-BIC(gls.TnDD);output[i,4]<-gls.TnDD$logLik; 
sdat<-predict(gls.TnDD);SCR2<-sum(gls.TnDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i

#Mod 7
Cand.models[[i]] <-gls.TnGS <- gls(Extinction ~ Arrival + Dfrac, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnGS, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnGS);output[i,3]<-BIC(gls.TnGS);output[i,4]<-gls.TnGS$logLik;
sdat<-predict(gls.TnGS);SCR2<-sum(gls.TnGS$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 


write.table(output,"M80%_AICweight_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
stat<-aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = TRUE)
write.table(stat,"M80%_ModelClassification_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
write.table(Vargls,"M80%_VarianceExplained_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);

##############################################################################################################
## REFINE USING BEST MODEL WITH DEVIANCE EXPLAINED > 80% + INTERACTIONS
##############################################################################################################
rm(list=ls(all=TRUE))
#library(rgdal)
library(fields)
library(nlme)
library(ncf)
library(AICcmodavg)
library(modEvA)

rm(list=ls(all=TRUE))
data <- read.table("Megafauna_(ClimateHuman)Bearing_v4(0kaLag_withSD).csv", sep=",", header=T) ### Modifier ici en header=T (ligne suivante aussi)
mdat<-mean(data$Extinction);SCR<-sum((data$Extinction - mdat)**2)#calcul de la somme des residue sur les donnees
output<-matrix(0,5,4);Vargls<-matrix(0,5,1);i<-1;Cand.models <- list( );i

#Mod 1
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 2
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP + (Arrival * Prec), data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 3
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP + (Arrival * EminP), data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 4
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP + (Prec * EminP), data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 5
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP + (Prec * EminP * Arrival), data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

write.table(output,"BestInteract_AICweight_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
stat<-aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = TRUE)
write.table(stat,"BestInteract_ModelClassification_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
write.table(Vargls,"BestInteract_VarianceExplained_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);

##############################################################################################################
## RELATIVE IMPORTANCE OF VARIABLE FOR EACH PREDICTOR IN BEST MODEL WITH DEVIANCE EXPLAINED > 80% + INTERACTIONS
##############################################################################################################
## CASE #1 => BASED ON AICc OUTPUTS (modele 1)
rm(list=ls(all=TRUE))
#library(rgdal)
library(fields)
library(nlme)
library(ncf)
library(AICcmodavg)
library(modEvA)

rm(list=ls(all=TRUE))
data <- read.table("Megafauna_(ClimateHuman)Bearing_v4(0kaLag_withSD).csv", sep=",", header=T) ### Modifier ici en header=T (ligne suivante aussi)
mdat<-mean(data$Extinction);SCR<-sum((data$Extinction - mdat)**2)#calcul de la somme des residue sur les donnees
output<-matrix(0,4,4);Vargls<-matrix(0,4,1);i<-1;Cand.models <- list( );i

#Mod 1
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 2
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Prec + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 3
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + EminP, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

#Mod 4
Cand.models[[i]] <-gls.TnPDD <- gls(Extinction ~ Arrival + Prec, data=data, correlation=corGaus(form=~Lon+Lat), method="ML",verbose=T);output[i,1]<-AICc(gls.TnPDD, return.K = FALSE, second.ord = TRUE, nobs = NULL);output[i,2]<-AIC(gls.TnPDD);output[i,3]<-BIC(gls.TnPDD);output[i,4]<-gls.TnPDD$logLik;
sdat<-predict(gls.TnPDD);SCR2<-sum(gls.TnPDD$residuals)**2;Vargls[i,1]<-100*(SCR-SCR2)/SCR; i<- i + 1;i 

write.table(output,"RelativeImportance_AICweight_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
stat<-aictab(Cand.models, modnames = NULL, second.ord = TRUE, nobs = NULL,sort = TRUE)
write.table(stat,"RelativeImportance_ModelClassification_stats_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
write.table(Vargls,"RelativeImportance_VarianceExplained_Bearing(HumanClimate)_v2.csv",sep = " ",dec = ".",row.names = FALSE);
