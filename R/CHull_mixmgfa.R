#' Convex Hull procedure for Mixture Multigroup Factor Analysis
# --------------------------------------------------------------------------------------
# Code written by Kim De Roover

#' @description
#' Performs the CHull for an overview table or a merge of muliple overview tables (e.g., when adding more numbers of clusters or repeating some analyses with more starts to avoid local maxima). The function automatically selects the most optimal solution for each number of clusters, reorders the rows and recomputes the CHull scree ratios.
#' @param x Outputobject$overview (class "mixmgfa") resulting from the mixmgfa function, or a merge of multiple such overview tables (vertically concatenated). The table should have 8 or 9 columns (so excluding or including the scree ratios). The rows don't have to be ordered according to the numbers of clusters.
#'
#'
#' @export
CHull_mixmgfa <- function(x){
  if(ncol(x)==9){
    overview=x[,c(1:6,8:9)] # remove current screeratios column, which needs to be re-done
  } else if(ncol(x)==8) {
    overview=x
  } else {
    stop("The number of columns is incorrect. Only input an overview table resulting from a mixmgfa analysis, or a vertical concatenation of multiple overview tables.")
  }
  nsclust=unique(overview[,1])
  indmax=matrix(0,length(nsclust),1) # for each unique number of clusters, retain the best fitting solution
  for(nclust in nsclust[1]:nsclust[length(nsclust)]){
    rows=which(overview[,1]==nclust)
    indmax[nclust-nsclust[1]+1]=rows[which.max(overview[rows,2])]
  }
  overview=overview[indmax,] # this sorts at the same time

  nsclust=c(nsclust[1],nsclust[length(nsclust)])
  if((max(overview))>0){
    nrows=nrow(overview)
    if(nrows>2){
      screeratios=matrix(NA,nsclust[2]-nsclust[1]+1,1)
      CHullcheck=0
      for(nclust in (nsclust[1]+1):(nsclust[2]-1)){
        LL_nclust=overview[nclust-nsclust[1]+1,2]
        npar_nclust=overview[nclust-nsclust[1]+1,3]
        LL_nclustmin1=overview[nclust-nsclust[1],2]
        npar_nclustmin1=overview[nclust-nsclust[1],3]
        LL_nclustplus1=overview[nclust-nsclust[1]+2,2]
        npar_nclustplus1=overview[nclust-nsclust[1]+2,3]
        # determine whether intermediate point is part of the convex hull
        slope=(LL_nclustplus1-LL_nclustmin1)/(npar_nclustplus1-npar_nclustmin1)
        point_line=LL_nclustmin1+slope*(npar_nclust-npar_nclustmin1)
        #diff_fit=(LL_nclust-LL_nclustmin1)/abs(LL_nclust)
        if(isTRUE(LL_nclust>=(point_line-.01))){ # && diff_fit>.0001)
          screeratios[nclust-nsclust[1]+1]=((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
        }
      }
      #screeratios[CHullcheck<0]=NA
      convexhull=which(!is.na(screeratios))
      #nrows=nrow(overview)
      convexhull=c(1,convexhull,nrows)
      nrhull=length(convexhull)
      change=0
      if(nrhull<nrows){
        change=1
      }
      while(nrhull>2 && change==1){ # check again whether intermediate points are on the convex hull
        nsclusthull=overview[convexhull,1]
        change=0
        for(indhull in 2:(nrhull-1)){
          if(!identical(convexhull[(indhull-1):(indhull+1)],c(convexhull[indhull]-1,convexhull[indhull],convexhull[indhull]+1))){
            LL_nclust=overview[convexhull[indhull],2]
            npar_nclust=overview[convexhull[indhull],3]
            LL_nclustmin1=overview[convexhull[indhull-1],2]
            npar_nclustmin1=overview[convexhull[indhull-1],3]
            LL_nclustplus1=overview[convexhull[indhull+1],2]
            npar_nclustplus1=overview[convexhull[indhull+1],3]
            # determine whether intermediate point is part of the convex hull
            slope=(LL_nclustplus1-LL_nclustmin1)/(npar_nclustplus1-npar_nclustmin1)
            point_line=LL_nclustmin1+slope*(npar_nclust-npar_nclustmin1)
            if(LL_nclust>=(point_line-.01)){
              # when the subset of three points spans across a point not on the hull, this is the corrected scree ratio (comparing with the previous and next point ON THE HULL)
              screeratios[convexhull[indhull]]=((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
            } else {
              screeratios[convexhull[indhull]]=NA
              change=1
            }
          }
        }
        convexhull=which(!is.na(screeratios))
        convexhull=c(1,convexhull,nrows)
        nrhull=length(convexhull)
      }
      overview=cbind(overview[,1:6],screeratios,overview[,7:8])
    }
    if(nrows>2){
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","AIC","screeratios","convergence","nr.activated.constraints")
    } else {
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","AIC","convergence","nr.activated.constraints")
    }
  }
  class(overview)<-"mixmgfaMS"
  return(overview)
}
