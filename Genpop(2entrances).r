#######################################################################################
humgen.fct <- function(mpts,maxr,maxc)
{
  curve=c(51.7,0.00007);
mtheo<-matrix(0,mpts,5);
#[,1]=row, [,2]=column, [,3]=age, [,4]=sd, [,5]theo of the gridcell

  mtheo[,1] <- runif(mpts)*maxr
  mtheo[,2] <- runif(mpts)*maxc
  xv <- dnorm(mtheo[,1],0.2,0.2)*dnorm(mtheo[,2],0.2,0.2)
  xv <- xv+dnorm(mtheo[,1],0.7,0.2)*dnorm(mtheo[,2],0.7,0.2)
  mtheo[,5] <- xv*80000/(dnorm(0,0,0.2)**2)
  mtheo[,3] <- runif(mpts)*mtheo[,5]
  mtheo[,4]=curve[1]*exp(curve[2]*mtheo[,3]);
  mtheo[,3]=mtheo[,3]+rnorm(nrow(mtheo),0,mtheo[,4]);
  for (i in 1:nrow(mtheo)){ mtheo[i,3]=max(mtheo[i,3],1)};
  return(mtheo);
  
}

