#######################################################################################
humgen.fct <- function(mpts,maxr,maxc)
{
  curve=c(51.7,0.00007);
  mtheo<-matrix(0,mpts,5);
#[,1]=row, [,2]=column, [,3]=age, [,4]=sd, [,5]theo of the gridcell

  mtheo[,1] <- runif(mpts)*maxr
  mtheo[,2] <- runif(mpts)*maxc
  xv <- mtheo[,1]/maxr+mtheo[,2]/maxc
  mtheo[,5] <- xv*40000
  mtheo[,3] <- runif(mpts)*mtheo[,5]
  mtheo[,4]=curve[1]*exp(curve[2]*mtheo[,3]);
  mtheo[,3]=mtheo[,3]+rnorm(nrow(mtheo),0,mtheo[,4]);
  for (i in 1:nrow(mtheo)){ mtheo[i,3]=max(mtheo[i,3],1)};

  return(mtheo);
}

####
## ajout
humgencond <- function(rtheo,ind)
  # pour un vecteur de valeurs theoriques rtheo,
  # pour une suite d'indices ind,
  # simulation d'ages et sd ages
{
  curve=c(51.7,0.00007)
  age <- runif(length(ind))*(-rtheo[ind]) 
  sdage <- curve[1]*exp(curve[2]*age)
  return(list(age=-age,sdage=sdage))
}
## END
##############################################

humgen.comp      <-  humgen.fct
humgencond.comp  <-  humgencond 
