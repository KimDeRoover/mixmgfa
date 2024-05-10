#' @export
do_mixmgfa_loadingsresiduals <- function(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = 5000,start = 1,nruns = 25,design=design,rotation=0,preselect=10,targetT=0,targetW=0,rescov=0,parcomp=0){
  mean_gs=data2$means
  data2=data2$covariances
  nvar=ncol(data2[[1]])
  G=length(N_gs)
  if(nclust==1){
    if(parcomp==0){
      cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
      cat("\n")
    }
    output_nclust <- mixmgfa_loadingsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect, rescov=rescov)
  }
  else {
    if(parcomp==0){
      cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
      cat("\n")
    }
    if(nclust==G){
      output_nclust <- mixmgfa_loadingsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect, rescov=rescov)
    }
    else{
      output_nclust <- mixmgfa_loadingsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect, rescov=rescov)
    }
  }
  # Perform rotations and organize and label the output
  if(sum(design)==0 && rotation!=0){
    for(k in 1:nclust){
      if(rotation=="varimax"){
        rot<-GPForth(output_nclust$Lambda_ks[[k]], method="varimax")
      } else if(rotation=="oblimin"){
        rot<-GPFoblq(output_nclust$Lambda_ks[[k]], method="oblimin")
      } else if(rotation=="target"){
        rot<-GPFoblq(output_nclust$Lambda_ks[[k]], method="pst",methodArgs =list(W=targetW,Target=targetT))
      }
      rotatedloadings=rot$loadings
      if(rotation=="varimax"){
        T_matrix=rot$Th
      } else {
        T_matrix=t(solve(rot$Th))
      }
      if(rotation!="target"){
        if(k==1){
          # ssq=colSums(rotatedloadings^2)
          # perm<-sort(ssq,decreasing=TRUE,index.return=TRUE)
          # perm<-perm$ix
          # rotatedloadings=rotatedloadings[,perm]
          strongloadcutoff=mean(apply(abs(rotatedloadings),2,max))/2
          strongloadind=which(abs(rotatedloadings)>=strongloadcutoff)
          strongload=matrix(0,nrow=nvar,ncol=nfactors)
          strongload[strongloadind]=rotatedloadings[strongloadind]
          toreflect=colSums(strongload>0)<colSums(strongload<0)
          refl=matrix(1,nrow=1,ncol=nfactors)
          refl[toreflect]=-1
          rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
          rotatedloadings_cl1=rotatedloadings
        } else {
          agreem_cl1=matrix(0,nfactors,nfactors)
          agreem_cl1_refl=agreem_cl1
          for(q in 1:nfactors){
            load=rotatedloadings[,q,drop=FALSE]
            for(q1 in 1:nfactors){
              agreem_cl1[q,q1]=sum(((load)-rotatedloadings_cl1[,q1])^2)
              agreem_cl1_refl[q,q1]=sum(((-1*load)-rotatedloadings_cl1[,q1])^2)
            }
          }
          mi=apply(agreem_cl1,2,sort)
          mi_ind=apply(agreem_cl1,2,order)
          mi_refl=apply(agreem_cl1_refl,2,sort)
          mi_refl_ind=apply(agreem_cl1_refl,2,order)
          refl=matrix(1,nrow=1,ncol=nfactors)
          refl[mi[1,]>mi_refl[1,]]=-1
          perm=mi_ind[1,]
          perm[refl==-1]=mi_refl_ind[1,refl==-1]
          if (length(unique(perm))<length(perm)){
            nonuniq=which(tabulate(perm)>1)
            missing=which(is.element(1:nfactors,perm)==FALSE)
            #perm[is.element(perm,nonuniq)]
            for(nu in 1:length(nonuniq)){
              fnu=nonuniq[nu]
              pos=which(perm==fnu)
              pos=pos[order(mi[1,perm==fnu])]
              perm[pos[2:length(pos)]]=0
            }
            perm[perm==0]=missing # replacement is not optimized at this time
          }
          rotatedloadings=rotatedloadings[,perm]
          rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
        }
        Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
      }
      output_nclust$Lambda_ks[[k]]=rotatedloadings

      invT=solve(T_matrix)
      for(g in 1:G){ # counter-rotate all sets of factor (co)variances and factor means
        crotatedFcov=invT%*%output_nclust$Phi_gks[[g,k]]%*%t(invT)
        if(k>1 && rotation!="target"){
          crotatedFcov=crotatedFcov[perm,perm]
        }
        if(rotation!="target"){
          crotatedFcov=Rm*crotatedFcov*t(Rm)
        }
        output_nclust$Phi_gks[[g,k]]=crotatedFcov
      }
    }
  }
  if(sum(design)>0 || nfactors==1 || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # in case of CFA or zero-approximating target rotation, only reflect
    for(k in 1:nclust){
      loadings=output_nclust$Lambda_ks[[k]]
      if(k==1){
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
      } else {
        agreem_cl1=matrix(0,nfactors,nfactors)
        agreem_cl1_refl=agreem_cl1
        for(q in 1:nfactors){
          load=loadings[,q,drop=FALSE]
          agreem_cl1[q,q]=sum(((load)-reflloadings_cl1[,q])^2)
          agreem_cl1_refl[q,q]=sum(((-1*load)-reflloadings_cl1[,q])^2)
        }
        refl=matrix(1,nrow=1,ncol=nfactors)
        refl[diag(agreem_cl1)>diag(agreem_cl1_refl)]=-1
        reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
      }
      output_nclust$Lambda_ks[[k]]=reflloadings

      Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
      for(g in 1:G){ # reflect all corresponding sets of factor (co)variances and factor means
        Fcov=output_nclust$Phi_gks[[g,k]]
        Fcov=Rm*Fcov*t(Rm)
        output_nclust$Phi_gks[[g,k]]=Fcov
      }
    }
  }
  output_nclust$mu_gs=mean_gs
  prefix="Cluster"
  suffix=seq(1,nclust)
  colnames(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
  names(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  names(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$Phi_gks)<-noquote(paste(prefix,suffix,sep="_"))
  colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
  names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
  if(!is.null(grouplabels)){
    rownames(output_nclust$z_gks)<-grouplabels
    rownames(output_nclust$mu_gs)<-grouplabels
    rownames(output_nclust$Phi_gks)<-grouplabels
  }
  if(!is.null(varlabels)){
    colnames(output_nclust$mu_gs)<-varlabels
    for(k in 1:nclust){
      rownames(output_nclust$Lambda_ks[[k]])<-varlabels
      colnames(output_nclust$Psi_ks[[k]])<-varlabels
      rownames(output_nclust$Psi_ks[[k]])<-varlabels
    }
  }
  prefix="Factor"
  suffix=seq(1,nfactors)
  factorlabels=noquote(paste(prefix,suffix,sep="_"))
  for(k in 1:nclust){
    colnames(output_nclust$Lambda_ks[[k]])<-factorlabels
    for(g in 1:G){
      colnames(output_nclust$Phi_gks[[g,k]])<-factorlabels
      rownames(output_nclust$Phi_gks[[g,k]])<-factorlabels
    }
  }
  output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,clusterspecific.uniquevariances=output_nclust$Psi_ks,groupspecific.means=output_nclust$mu_gs,loglik=output_nclust$bestloglik,nrpars=output_nclust$nrpars,convergence=output_nclust$convergence>0,nractivatedconstraints=output_nclust$nractivatedconstraints)
  return(output_nclust2)
}
