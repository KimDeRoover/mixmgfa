#' Model selection plots for mixture multigroup factor analyses
#' @param x Outputobject$overview (class "mixmgfaMS") resulting from mixmgfa function


#' @export
plot.mixmgfaMS <- function(x, ...){
  nrows=nrow(x)
  par(mfrow=c(2,1))
  par(mar=c(4.1,4.1,4.1,2.1))
  plot(x[,1],x[,5],type="b",xlab="number of clusters",ylab="BIC_G")
  #plot(x[,1],x[,2],type="b",xlab="number of clusters",ylab="loglikelihood")
  convexhull=which(!is.na(x[,6]))
  convexhull=c(x[1,1],convexhull,x[nrows,1])
  notonhull=x[1,1]:x[nrows,1]
  notonhull[notonhull!=convexhull]
  plot(x[convexhull,3],x[convexhull,2],type="b",xlab="number of free parameters",ylab="loglikelihood")
  points(x[notonhull,3],x[notonhull,2])
  prefix=seq(x[1,1]:x[nrows,1])
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
  mtext("Model selection plots for mixture multigroup factor analyses", side = 3, line = -2, outer = TRUE)

}
