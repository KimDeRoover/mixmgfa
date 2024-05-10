#' Model selection plots for mixture multigroup factor analyses
#' @param x Outputobject$overview (class "mixmgfaMS") resulting from mixmgfa function
#' @param AIC Specify AIC = 1 if you would like to plot the AIC instead of the BIC_G in the first plot


#' @export
plot.mixmgfaMS <- function(x,AIC=0, ...){
  nrows=nrow(x)
  par(mfrow=c(2,1))
  par(mar=c(4.1,4.1,4.1,2.1))
  if(AIC==0){
    plot(x[,1],x[,5],type="b",xlab="number of clusters",ylab="BIC_G")
  } else {
    plot(x[,1],x[,6],type="b",xlab="number of clusters",ylab="AIC")
  }

  #plot(x[,1],x[,2],type="b",xlab="number of clusters",ylab="loglikelihood")
  if(nrows>2){
    convexhull=which(!is.na(x[,7]))
    convexhull=c(1,convexhull,nrows)
    notonhull=1:nrows
    notonhull=notonhull[notonhull!=convexhull]
    plot(x[convexhull,3],x[convexhull,2],type="b",xlab="number of free parameters",ylab="loglikelihood")
    points(x[notonhull,3],x[notonhull,2])
    prefix=seq(x[1,1],x[nrows,1])
    suffix="clusters"
    pos_vector <- rep(3, nrows)
    pos_vector[1] <- 4
    pos_vector[nrows-1] <- 1
    pos_vector[nrows] <- 2
    labels<-paste(prefix,suffix,sep=" ")
    if(x[1,1]==1){
      labels[1]="1 cluster"
    }
    text(x[,3],x[,2],labels,pos=pos_vector)
  }
  mtext("Model selection plots for mixture multigroup factor analyses", side = 3, line = -2, outer = TRUE)

}
