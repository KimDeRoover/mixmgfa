#' @export
do_mixmgfa_interceptsresiduals <- function(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = 5000,start = 1,nruns = 25,design=design,rotation=0,preselect=10,targetT=0,targetW=0,parcomp=0){
  nvar <- ncol(data2$covariances[[1]])
  G=length(N_gs)
  if(nclust==1){
    if(parcomp==0){
      cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
      cat("\n")
    }
    output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
  }
  else {
    if(parcomp==0){
      cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
      cat("\n")
    }
    if(nclust==G){
      output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
    }
    else{
      output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
    }
  }
  # Perform rotations and organize and label the output
  if(sum(design)==0 && rotation!=0){
    if(rotation=="varimax"){
      rot<-GPForth(output_nclust$Lambda, method="varimax")
    } else if(rotation=="oblimin"){
      rot<-GPFoblq(output_nclust$Lambda, method="oblimin")
    } else if(rotation=="target"){
      rot<-GPFoblq(output_nclust$Lambda, method="pst",methodArgs =list(W=targetW,Target=targetT))
    }
    rotatedloadings=rot$loadings
    if(rotation=="varimax"){
      T_matrix=rot$Th
    } else {
      T_matrix=t(solve(rot$Th))
    }
    if(rotation!="target"){ # reflect in case of rotated EFA (not for target rotation)
      strongloadcutoff=mean(apply(abs(rotatedloadings),2,max))/2
      strongloadind=which(abs(rotatedloadings)>=strongloadcutoff)
      strongload=matrix(0,nrow=nvar,ncol=nfactors)
      strongload[strongloadind]=rotatedloadings[strongloadind]
      toreflect=colSums(strongload>0)<colSums(strongload<0)
      refl=matrix(1,nrow=1,ncol=nfactors)
      refl[toreflect]=-1
      rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
      rotatedloadings_cl1=rotatedloadings
      Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
    }
    output_nclust$Lambda=rotatedloadings

    invT=solve(T_matrix)
    for(g in 1:G){ # counter-rotate all sets of factor (co)variances and factor means
      crotatedFcov=invT%*%output_nclust$Phi_gs[[g]]%*%t(invT)
      if(rotation!="target"){
        crotatedFcov=Rm*crotatedFcov*t(Rm)
      }
      output_nclust$Phi_gs[[g]]=crotatedFcov
      for(k in 1:nclust){
        crotatedFmeans=output_nclust$alpha_gks[[g,k]]%*%t(invT)
        if(rotation!="target"){
          crotatedFmeans=refl*crotatedFmeans
        }
        output_nclust$alpha_gks[[g,k]]=crotatedFmeans
      }
    }
  }
  if(sum(design)>0 || nfactors==1 || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # in case of CFA or zero-approximating target rotation, only reflect
    loadings=output_nclust$Lambda
    if(rotation=="target"){
      strongloadcutoff=mean(apply(abs(loadings),2,max))/2
      strongloadind=which(abs(loadings)>=strongloadcutoff)
      strongload=matrix(0,nrow=nvar,ncol=nfactors)
      strongload[strongloadind]=loadings[strongloadind]
    } else {
      strongload=loadings # for CFA, all non-zero loadings are considered strong loadings
    }
    toreflect=colSums(strongload>0)<colSums(strongload<0)
    refl=matrix(1,nrow=1,ncol=nfactors)
    refl[toreflect]=-1
    reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
    reflloadings_cl1=reflloadings
    output_nclust$Lambda=reflloadings

    Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
    for(g in 1:G){ # reflect all corresponding sets of factor (co)variances and factor means
      Fcov=output_nclust$Phi_gs[[g]]
      Fcov=Rm*Fcov*t(Rm)
      output_nclust$Phi_gs[[g]]=Fcov
      for(k in 1:nclust){
        Fmeans=output_nclust$alpha_gks[[g,k]]
        Fmeans=refl*Fmeans
        output_nclust$alpha_gks[[g,k]]=Fmeans
      }
    }
  }

  prefix="Cluster"
  suffix=seq(1,nclust)
  colnames(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  names(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
  names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  if(is.null(grouplabels)==FALSE){
    rownames(output_nclust$z_gks)<-grouplabels
    rownames(output_nclust$alpha_gks)<-grouplabels
    rownames(output_nclust$Phi_gs)<-grouplabels
  }
  if(is.null(varlabels)==FALSE){
    colnames(output_nclust$tau_ks)<-varlabels
    rownames(output_nclust$Lambda)<-varlabels
    for(k in 1:nclust){
      colnames(output_nclust$Psi_ks[[k]])<-varlabels
      rownames(output_nclust$Psi_ks[[k]])<-varlabels
    }
  }
  prefix="Factor"
  suffix=seq(1,nfactors)
  factorlabels=noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$Lambda)<-factorlabels
  for(g in 1:G){
    colnames(output_nclust$Phi_gs[[g]])<-factorlabels
    rownames(output_nclust$Phi_gs[[g]])<-factorlabels
    for(k in 1:nclust){
      colnames(output_nclust$alpha_gks[[g,k]])<-factorlabels
    }
  }
  output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,invariant.loadings=output_nclust$Lambda,groupspecific.factorcovariances=output_nclust$Phi_gs,clusterspecific.uniquevariances=output_nclust$Psi_ks,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks,loglik=output_nclust$bestloglik,nrpars=output_nclust$nrpars,convergence=output_nclust$convergence>0,nractivatedconstraints=output_nclust$nractivatedconstraints)
  return(output_nclust2)
}
