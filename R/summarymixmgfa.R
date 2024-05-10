#' Creates summary of clustering and measurement parameters for estimated MMG-FA model
#' @param x Outputobject$MMGFAsolutions (class "mixmgfa") resulting from mixmgfa function
#' @param nclust Number of clusters you want to create the summary for
#'
#'
#' @export
summary.mixmgfa <- function(x,nclust=0, ...){
  options(scipen = 6)
  if(nclust==0){
    cat("Please specify a number of clusters as the second argument of 'summary'")
  } else {
    # if(x$overview[nclust,7]==1){
    #   cat("Estimation converged within maximal number of iterations.")
    # } else {
    #   cat("Estimation did not converge within maximal number of iterations: please increase 'maxiter'.")
    # }
    # cat("\n")
    # cat("Loglikelihood: ")
    # cat(round(x$overview[nclust,2],1))
    # cat("\n")
    cat("cluster memberships of the groups:")
    cat("\n")
    print(round(x[[nclust]]$clustermemberships,4))
    cat("\n")
    cat("\n")
    ngroups=nrow(x[[nclust]]$clustermemberships)

    modalassign=apply(x[[nclust]]$clustermemberships,1,which.max)
    for(k in 1:nclust){
      cl_names<-names(which(modalassign==k))#names(which(x[[nclust]]$clustermemberships[,k]>(1/nclust)))
      cat(paste("Groups modally assigned to cluster ",k,":",sep=""))
      cat("\n")
      cat(cl_names)
      cat("\n")
    }
    cat("\n")
    cat("\n")

    cat("cluster proportions:")
    cat("\n")
    cat("\n")
    print(round(x[[nclust]]$clusterproportions,4))
    cat("\n")
    cat("\n")

    if(!is.null(x[[nclust]]$clusterspecific.loadings)){
      cat("cluster-specific loadings:")
      cat("\n")
      cat("\n")
      #Loadingstable=matrix(unlist(x[[nclust]]$clusterspecific.loadings),nrow = nrow(x[[nclust]]$clusterspecific.loadings[[1]]),ncol=ncol(x[[nclust]]$clusterspecific.loadings[[1]])*nclust)
      #rownames(Loadingstable)<-rownames(x[[nclust]]$clusterspecific.loadings[[1]])
      nvar<-nrow(x[[nclust]]$clusterspecific.loadings[[1]])
      nfactors<-ncol(x[[nclust]]$clusterspecific.loadings[[1]])
      prefix="Fact"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep=""))
      for(k in 1:nclust){
        loadings_cl <- cbind.data.frame(round(x[[nclust]]$clusterspecific.loadings[[k]],4),"  " = c(rep("  ",nvar)))
        suffix=paste("Cl",k,sep="")
        factorlabels2<-t(paste(factorlabels,suffix,sep=""))
        labels<-cbind(factorlabels2,"  ")
        colnames(loadings_cl)<-labels
        if(k==1){
          Loadingstable <- loadings_cl
        } else {
          Loadingstable <- cbind.data.frame(Loadingstable,loadings_cl)
        }
      }
      print(Loadingstable)
    } else {
      cat("invariant loadings:")
      cat("\n")
      cat("\n")
      nvar<-nrow(x[[nclust]]$invariant.loadings)
      nfactors<-ncol(x[[nclust]]$invariant.loadings)
      prefix="Fact"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep=""))
      Loadingstable <- as.data.frame(round(x[[nclust]]$invariant.loadings,4))
      colnames(Loadingstable)<-factorlabels
      print(Loadingstable)
    }
    cat("\n")
    cat("\n")


    if(!is.null(x[[nclust]]$groupspecific.means)){
      cat("group-specific intercepts (= means):")
      cat("\n")
      cat("\n")
      print(round(x[[nclust]]$groupspecific.means,4))
    } else if(!is.null(x[[nclust]]$clusterspecific.intercepts)){
      cat("cluster-specific intercepts:")
      cat("\n")
      cat("\n")
      print(round(x[[nclust]]$clusterspecific.intercepts,4))
    } else if(!is.null(x[[nclust]]$invariant.intercepts)){
      cat("invariant intercepts:")
      cat("\n")
      cat("\n")
      print(round(x[[nclust]]$invariant.intercepts,4))
    }
    cat("\n")
    cat("\n")

    if(!is.null(x[[nclust]]$groupspecific.uniquevariances)){
      cat("group-specific unique variances:")
      cat("\n")
      cat("\n")
      uniquevariances=x[[nclust]]$groupspecific.uniquevariances
    } else if (!is.null(x[[nclust]]$clusterspecific.uniquevariances)){
      cat("cluster-specific unique variances:")
      cat("\n")
      cat("\n")
      uniquevariances=x[[nclust]]$clusterspecific.uniquevariances
    }
    if (sum(uniquevariances[[1]]!=0)==nvar){ # no residual covariances
      #print(round(t(as.data.frame(lapply(uniquevariances,diag))),4))
      l=length(uniquevariances)
      uniquevariances=lapply(uniquevariances,diag)
      varnames=names(uniquevariances[[1]])
      uniquevariances=lapply(uniquevariances,round,digits=4)
      uniquevariances=matrix(unlist(uniquevariances,use.names = FALSE),nrow=l,ncol=nvar,byrow = TRUE)
      colnames(uniquevariances)=varnames
      if(!is.null(x[[nclust]]$groupspecific.uniquevariances)){
        groupnames=rownames(x[[nclust]]$clustermemberships)
        rownames(uniquevariances)=groupnames
      }
      print(as.data.frame(uniquevariances))
    }
    cat("\n")
  }
}
