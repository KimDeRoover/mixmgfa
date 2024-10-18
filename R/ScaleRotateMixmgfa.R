#' Re-scaling and rotation function for output objects of mixture multigroup factor analysis
# -----------------------------------------------------------------------------------------
# Code written by Kim De Roover

#' @description
#' Allows to rotate and/or re-scale (to marker variable scale) mixture multigroup factor analysis solutions. In this way, you don't need to re-estimate the solutions when you want to change the scale or rotations.
#' It also counter-rotates the factor (co)variances and factor means.

#' @param OutputObject Output object resulting from using the mixmgfa function, either the entire output object or the MMGFAsolutions part. This may be rotated or unrotated.
#' @param N_gs Vector with number of subjects (sample size) for each group (in the same order as they appear in the data).
#' @param cluster.spec Measurement parameters you clustered the groups on: "loadings", "intercepts", "residuals", c("loadings","intercepts"), c("intercepts","residuals"), c("loadings","residuals"), or c("loadings","intercepts","residuals").
#' @param nsclust Vector of length one or two indicating the number(s) of clusters you want to rotate or re-scale for. In case of length two, the vector indicates the minimal and maximal number of clusters. If left unspecified, nsclust is derived from OutputObject.
#' @param design For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).
#' @param rotation Rotation criterion to use in case of EFA; currently either "oblimin", "geomin", "varimax" or "target" (i.e., semi-specified oblique Procrustes rotation), whereas 0 = no rotation. (Note: The GPArotation package is loaded or installed for rotation.)
#' @param rescale Equal to 1 when solutions need to be rescaled.
#' @param markers Matrix (with ncol = nfactors and nrow = number of items) indicating position of marker variable loadings with '1', wheres the rest of the matrix should be equal to '0'. Note that marker variable re-scaling does not make a lot of sense for an unrotated EFA solution.
#' @param targetT Target matrix to use when rotation = "target". Note that targetT can be one target matrix (with nrow = number of items) or a vertical concatenation of cluster-specific target matrices (with nrow = number of items x number of clusters).
#' Note that, when nsclust contains >1 numbers of clusters, the same target matrix (i.e., the first one) is used for all clusters.
#' @param targetW Weights to be used when rotation = "target". You can set an entire row to zero to make sure that the simple structure is not messed up by a 'complex variable' (i.e., with strong loadings for multiple factors). Set all weights equal to '1' if you prefer fully specified target rotation.
#' When left unspecified while rotation = target, the default is a targetW where the zeros in the target get a weight of '1' and the non-zeros in the target get a weight of '0'. If this results in all zero weights or too many zeros for the rotation to be identified, a fully specified target rotation (all weights equal to '1') is used instead.

#' @export
ScaleRotateMixmgfa <- function(OutputObject,N_gs,cluster.spec,nsclust=c(),design=0,rescale=0,markers=0,rotation=0,targetT=0,targetW=0){

  if(length(OutputObject)==2){
    Output=OutputObject$MMGFAsolutions
  } else {
    Output=OutputObject
  }

  if(is.null(nsclust)){
    nclust_max=length(Output[[length(Output)]]$clusterproportions)
    nclust_min=length(Output[[1]]$clusterproportions)
    nsclust=c(nclust_min,nclust_max)
  }


  # check whether 'cluster.specific' is correctly specified
  if(is.element("loadings",cluster.spec) && is.null(Output[[1]]$clusterspecific.loadings)){
    stop("The specification of cluster.spec does not match the OutputObject. Please specify cluster.spec correctly and try again.")
  }
  if(is.element("intercepts",cluster.spec) && is.null(Output[[1]]$clusterspecific.intercepts)){
    stop("The specification of cluster.spec does not match the OutputObject. Please specify cluster.spec correctly and try again.")
  }

  if(is.element("loadings",cluster.spec)){
    nfactors=ncol(Output[[1]]$clusterspecific.loadings[[1]])
  } else {
    nfactors=ncol(Output[[1]]$invariant.loadings)
  }

  if(max(design)!=0){
    if(ncol(design)!=nfactors){
      stop("The specified design matrix does not match the number of factors.")
    }
  }
  if(sum(design)==0){
    EFA=1
  } else {
    EFA=0
  }

  if(rotation==1 && (nfactors==1 || EFA==0)){
    cat("You requested rotation for a one-factor or CFA solution. Rotation will not be performed.")
    cat("\n")
    rotation=0
  }

  if(rotation!=0){
    if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
      rotation="varimax"
    } else if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
      rotation="oblimin"
    } else if(rotation=="geomin" || rotation=="Geomin" || rotation=="GEOMIN"){
      rotation="geomin"
    } else if(rotation=="target" || rotation=="Target" || rotation=="TARGET"){
      rotation="target"
    } else {
      cat("You seem to have misspelled your rotation criterion. Oblimin rotation will be performed.")
      cat("\n")
      rotation="oblimin"
    }
  }
  if(rotation!=0){
    if(rotation=="target"){
      if(max(targetT)==0){
        stop("You failed to specify a targetT matrix or it contains only zeros. Please specify a valid targetT matrix")
      } else if(max(targetT)!=0){
        if(ncol(targetT)!=nfactors){
          stop("The specified targetT matrix does not match the specified number of factors.")
        }
      }
      if(max(targetW)==0){
        targetW=1-ceiling(targetT)
        cat("Note: targetW is unspecified. A weight matrix is derived from targetT so that the zeros in targetT are approximated.")
        cat("\n")
        if(sum(targetW)<(nfactors*(nfactors-1))){
          targetW=matrix(1,ncol(targetT),nfactors)
          cat("Note: targetW is unspecified. A valid weight matrix could not be derived from targetT, so fully specified target rotation will be performed.")
          cat("\n")
        }
      }
    }
    if(!require(GPArotation)){ # test if package is loaded
      loadtest<-try(library(GPArotation),silent=TRUE) # test if package can be loaded (if it is installed)
      if(class(loadtest)=="try-error"){
        install.packages("GPArotation")
      }
    }
  }

  if(rescale==1 && max(markers)==0){
    if(EFA==0){
      cat("You requested rescaling for a CFA solution, but did not specify marker variables. Default marker variables will be derived from the design matrix.")
      cat("\n")
      markers=matrix(0,nrow(design),ncol(design))
      ssitems=rowSums(design)==1 # items with simple structure in the design
      ssitems=which(ssitems==1)
      markerind=ssitems[apply(design[ssitems,],2,which.max)]
      markers[markerind,seq_len(nfactors)]=1
      markers=markers*design
    } else {
      cat("You requested rescaling for an EFA solution, but did not specify marker variables. Rescaling will not be performed.")
      cat("\n")
      rescale=0
    }
  } else if(rescale==1 && (EFA==0 || max(targetT)!=0)){ # checking whether markers match design or targetT
    markeritems=which(rowSums(markers)>0)
    if(EFA==0){
      indhl=apply(design[markeritems,],1,which.max) # index of factor each marker item has non-zero loading on according to 'index'
    } else if (max(targetT)!=0){
      indhl=apply(targetT[markeritems,],1,which.max) # index of factor each marker item has non-zero loading on according to 'index'
    }
    ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
    permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
    markers=markers[,permm]
  }

  if(is.element("loadings",cluster.spec)){
    ngroup=nrow(Output[[1]]$group.and.clusterspecific.factorcovariances)
  } else {
    ngroup=nrow(Output[[1]]$groupspecific.factorcovariances)
  }

  if(length(nsclust)>2){ # nsclust needs to be of length two
    cat("You specified more than two numbers of clusters. The minimum and maximum of these numbers will be used.")
    cat("\n")
    nsclust=c(min(nsclust),max(nsclust))
  }
  if(length(nsclust)==2){
    if(nsclust[2]>ngroup){
      cat("Note: the number of clusters cannot exceed the number of groups.")
      cat("\n")
      nsclust[2]=ngroup
    }
  } else if(length(nsclust)==1){
    if(nsclust[1]>ngroup){
      cat("Note: the number of clusters cannot exceed the number of groups.")
      cat("\n")
      nsclust[1]=ngroup
    }
    nsclust=c(nsclust,nsclust)
  }
  # Check whether the specified numbers of clusters match the ones in the Output object
  nclust_max=length(Output[[length(Output)]]$clusterproportions)
  nclust_min=length(Output[[1]]$clusterproportions)
  correction=0
  if(nsclust[1]<nclust_min){
    nsclust[1]=nclust_min
    correction=1
  }
  if(nsclust[2]>nclust_max){
    nsclust[2]=nclust_max
    correction=1
  }
  if(correction==1){
    cat("Note: the specified number(s) of clusters exceed the number(s) of clusters in the Output object. This is corrected.")
    cat("\n")
  }

  N=sum(N_gs)
  # start of loop over numbers of clusters
  for(nclust in nsclust[1]:nsclust[2]){
    ind_nclust=nclust-nclust_min+1
    # if rotation is requested in case of EFA, make factors orthogonal
    # per cluster or overall
    if(EFA==1 && rotation!=0){
      if(is.element("loadings",cluster.spec)){
        z_gks=Output[[ind_nclust]]$clustermemberships
        pi_ks=Output[[ind_nclust]]$clusterproportions
        N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
        N_ks=colSums(N_gks)
        for(k in 1:nclust){
          if(N_ks[k]>0){
            phi_k=matrix(0,nfactors,nfactors)
            for(g in 1:ngroup){
              phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
              phi_k=phi_k+(N_gks[g,k]/N_ks[k])*phi_gk;
            }
            # find matrix square root via eigenvalue decomposition
            ed=eigen(phi_k)
            sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
            invsqrtFscale=solve(sqrtFscale) # This normally has a diagonal of ones, so that it only orthogonalizes, without rescaling
            for(g in 1:ngroup){
              phi_gk=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
              phi_gk=invsqrtFscale%*%phi_gk%*%invsqrtFscale;
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry
            }
            Output[[ind_nclust]]$clusterspecific.loadings[[k]]=Output[[ind_nclust]]$clusterspecific.loadings[[k]]%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
            if(is.element("intercepts",cluster.spec)){
              for(g in 1:ngroup){ # transform the factor means accordingly
                Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%invsqrtFscale
              }
            }
          }
        }
      } else {
        phi=matrix(0,nfactors,nfactors)
        for(g in 1:ngroup){
          phi_g=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
          phi=phi+(N_gs[g]/N)*phi_g
        }
        # find matrix square root via eigenvalue decomposition
        ed=eigen(phi)
        sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
        invsqrtFscale=solve(sqrtFscale)
        for(g in 1:ngroup){
          phi_g=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
          phi_g=invsqrtFscale%*%phi_g%*%invsqrtFscale;
          Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=((phi_g+t(phi_g))/2); # enforce perfect symmetry
        }
        Output[[ind_nclust]]$invariant.loadings=Output[[ind_nclust]]$invariant.loadings%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
        if(is.element("intercepts",cluster.spec)){
          for(g in 1:ngroup){ # transform the factor means accordingly
            for(k in 1:nclust){
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%invsqrtFscale
            }
          }
        } else {
          for(g in 1:ngroup){ # transform the factor means accordingly
            Output[[ind_nclust]]$groupspecific.factormeans[[g]]=Output[[ind_nclust]]$groupspecific.factormeans[[g]]%*%invsqrtFscale
          }
        }
      }
    }



    if(is.element("loadings",cluster.spec)){ # rotation/rescaling needs to be performed per cluster
      if(sum(design)==0){ # EFA
        for(k in 1:nclust){
          if(rotation!=0){
            if(rotation=="varimax"){
              rot<-GPForth(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="varimax")
            } else if(rotation=="oblimin"){
              rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="oblimin")
            } else if(rotation=="geomin"){
              rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="geomin")
            } else if(rotation=="target"){
              rot<-GPFoblq(Output[[ind_nclust]]$clusterspecific.loadings[[k]], method="pst",methodArgs =list(W=targetW,Target=targetT))
            }
            rotatedloadings=rot$loadings
            nvar=nrow(rotatedloadings)
            if(rotation=="varimax"){
              T_matrix=rot$Th
            } else {
              T_matrix=t(solve(rot$Th))
            }
            invT=solve(T_matrix)
            if(rotation!="target"){ # reflect and permute in case of rotated EFA (but not target rotation)
              if(k==1){
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
            } else {
              # refl en rm?
            }
          } else {
            rotatedloadings=Output[[ind_nclust]]$clusterspecific.loadings[[k]]
            refl=matrix(1,nrow=1,ncol=nfactors)
            Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
            invT=diag(nfactors)
          }

          if(rescale==1){ # determine rescaling constant per factor, taking into account reflection and permutation
            if(k==1){
              if(rotation!="target"){ # also when rotation=0 (object can/should be previously rotated)
                markeritems=which(rowSums(markers)>0)
                indhl=apply(rotatedloadings[markeritems,],1,which.max) # index of factor each marker item has highest loading on
                ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
                permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
                if (length(unique(permm))<length(permm)){ # check whether this includes no duplicate elements and thus pertains to a proper permutation
                  nonuniq=which(tabulate(permm)>1)
                  missing=which(is.element(1:nfactors,permm)==FALSE)
                  for(nu in 1:length(nonuniq)){
                    fnu=nonuniq[nu] # non-unique element
                    pos=which(permm==fnu) # positions of non-unique element
                    load=rotatedloadings*markers[,permm] # corresponding loading values
                    load=load[markeritems[fnu],]
                    pos=pos[order(load[permm==fnu],decreasing=TRUE)]
                    permm[pos[2:length(pos)]]=0 # only keep element at position corresponding to highest loading
                  }
                  permm[permm==0]=missing # replacement is not optimized at this time
                }
                markers2=markers[,permm] # Note that markers2 is also used for k>1
              } else { # for target rotation, the permutation of 'markers' is already checked at the beginning, in correspondence to targetT
                markers2=markers
              }
            }
            markload=rotatedloadings*markers2 # rescaling also happens for fully specified target rotation
            markload=markload[markload!=0]
            sqrtFscale=diag(1/markload)
            invsqrtFscale=diag(markload)
            rotatedloadings=rotatedloadings%*%sqrtFscale
          }
          Output[[ind_nclust]]$clusterspecific.loadings[[k]]=rotatedloadings

          if(is.element("intercepts",cluster.spec)){ # factor means also need to be counterrotated/rescaled
            for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
              crotatedFcov=invT%*%Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]%*%t(invT)
              crotatedFmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%t(invT)
              if(k>1 && rotation!="target" && rotation!=0){
                crotatedFcov=crotatedFcov[perm,perm]
                crotatedFmeans=crotatedFmeans[,perm]
              }
              if(rotation!="target"){
                crotatedFcov=Rm*crotatedFcov*t(Rm)
                crotatedFmeans=refl*crotatedFmeans
              }
              if(rescale==1){ # correct for rescaling as well
                crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
                crotatedFmeans=crotatedFmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=crotatedFcov
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=crotatedFmeans
            }
          } else {
            for(g in 1:ngroup){ # counter-rotate all corresponding sets of factor (co)variances
              crotatedFcov=invT%*%Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]%*%t(invT)
              if(k>1 && rotation!="target" && rotation!=0){
                crotatedFcov=crotatedFcov[perm,perm]
              }
              if(rotation!="target"){
                crotatedFcov=Rm*crotatedFcov*t(Rm)
              }
              if(rescale==1){
                crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=crotatedFcov
            }
          }
        }
      }
      if(EFA==0 || nfactors==1 || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # in case of CFA or zero-approximating target rotation, only reflect
        for(k in 1:nclust){
          loadings=Output[[ind_nclust]]$clusterspecific.loadings[[k]]
          nvar=nrow(loadings)
          if(k==1){ # reflect to have positive strong loadings
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
          } else { # reflect to agreement with cluster 1
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
          if(rescale==1 && EFA==0){ # determine rescaling constant per factor, taking into account reflection
            # for EFA with target rotation, rescaling already happened above
            markload=reflloadings*markers
            markload=markload[markload!=0]
            sqrtFscale=diag(1/markload)
            invsqrtFscale=diag(markload)
            reflloadings=reflloadings%*%sqrtFscale
          }
          Output[[ind_nclust]]$clusterspecific.loadings[[k]]=reflloadings
          Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
          for(g in 1:ngroup){ # reflect all corresponding sets of factor (co)variances
            Fcov=Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]
            Fcov=Rm*Fcov*t(Rm)
            if(rescale==1 && EFA==0){
              Fcov=invsqrtFscale%*%Fcov%*%invsqrtFscale
            }
            Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]]=Fcov
            if(is.element("intercepts",cluster.spec)){ # factor means also need to be counterrotated/rescaled
              Fmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]
              Fmeans=refl*Fmeans
              if(rescale==1 && EFA==0){
                Fmeans=Fmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Fmeans
            }
          }
        }
      }

      prefix="Factor"
      suffix=seq(1,nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      for(k in 1:nclust){
        colnames(Output[[ind_nclust]]$clusterspecific.loadings[[k]])<-factorlabels
        for(g in 1:ngroup){
          colnames(Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]])<-factorlabels
          rownames(Output[[ind_nclust]]$group.and.clusterspecific.factorcovariances[[g,k]])<-factorlabels
          if(is.element("intercepts",cluster.spec)){
            colnames(Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]])<-factorlabels
          }
        }
      }

    } else { # rotation/rescaling needs to be performed once, across all groups
      if(sum(design)==0){
        if(rotation!=0){
          if(rotation=="varimax"){
            rot<-GPForth(Output[[ind_nclust]]$invariant.loadings, method="varimax")
          } else if(rotation=="oblimin"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="oblimin")
          } else if(rotation=="geomin"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="geomin")
          } else if(rotation=="target"){
            rot<-GPFoblq(Output[[ind_nclust]]$invariant.loadings, method="pst",methodArgs =list(W=targetW,Target=targetT))
          }
          rotatedloadings=rot$loadings
          nvar=nrow(rotatedloadings)
          if(rotation=="varimax"){
            T_matrix=rot$Th
          } else {
            T_matrix=t(solve(rot$Th))
          }
          invT=solve(T_matrix)
          if(rotation!="target" || (rotation=="target" & sum(targetT==0 & targetW==1)==sum(targetW))){ # reflect in case of rotated EFA (not for target rotation, unless it only approximates zeros)
            strongloadcutoff=mean(apply(abs(rotatedloadings),2,max))/2
            strongloadind=which(abs(rotatedloadings)>=strongloadcutoff)
            strongload=matrix(0,nrow=nvar,ncol=nfactors)
            strongload[strongloadind]=rotatedloadings[strongloadind]
            toreflect=colSums(strongload>0)<colSums(strongload<0)
            refl=matrix(1,nrow=1,ncol=nfactors)
            refl[toreflect]=-1
            rotatedloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*rotatedloadings
          } else {
            refl=matrix(1,nrow=1,ncol=nfactors)
          }
        } else {
          rotatedloadings=Output[[ind_nclust]]$invariant.loadings
          refl=matrix(1,nrow=1,ncol=nfactors)
          invT=diag(nfactors)
        }
        Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
        if(rescale==1){
          if(rotation!="target"){
            markeritems=which(rowSums(markers)>0)
            indhl=apply(rotatedloadings[markeritems,],1,which.max) # index of factor each marker item has highest loading on
            ind1=apply(markers[markeritems,],1,which.max) # index of factor each item has a 1 for in 'markers'
            permm=match(ind1,indhl) # match the two to check whether permutation of 'markers' is needed
            if (length(unique(permm))<length(permm)){ # check whether this includes no duplicate elements and thus pertains to a proper permutation
              nonuniq=which(tabulate(permm)>1)
              missing=which(is.element(1:nfactors,permm)==FALSE)
              for(nu in 1:length(nonuniq)){
                fnu=nonuniq[nu] # non-unique element
                pos=which(permm==fnu) # positions of non-unique element
                load=rotatedloadings*markers[,permm] # corresponding loading values
                load=load[markeritems[fnu],]
                pos=pos[order(load[permm==fnu],decreasing=TRUE)]
                permm[pos[2:length(pos)]]=0 # only keep element at position corresponding to highest loading
              }
              permm[permm==0]=missing # replacement is not optimized at this time
            }
            markers2=markers[,permm]
          } else {
            markers2=markers
          }
          markload=rotatedloadings*markers2
          markload=markload[markload!=0]
          sqrtFscale=diag(1/markload)
          invsqrtFscale=diag(markload)
          rotatedloadings=rotatedloadings%*%sqrtFscale
        }
        Output[[ind_nclust]]$invariant.loadings=rotatedloadings

        if(is.element("intercepts",cluster.spec)){ # intercepts are cluster-specific, so factor means are group-and-cluster-specific
          for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
            crotatedFcov=invT%*%Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]%*%t(invT)
            crotatedFcov=Rm*crotatedFcov*t(Rm)
            if(rescale==1){
              crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
            }
            Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=crotatedFcov
            for(k in 1:nclust){
              crotatedFmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]%*%t(invT)
              crotatedFmeans=refl*crotatedFmeans
              if(rescale==1){
                crotatedFmeans=crotatedFmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=crotatedFmeans
            }
          }
        } else { # intercepts are invariant, so factor means are group-specific
          for(g in 1:ngroup){ # counter-rotate all sets of factor (co)variances and factor means
            crotatedFcov=invT%*%Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]%*%t(invT)
            crotatedFmeans=Output[[ind_nclust]]$groupspecific.factormeans[[g]]%*%t(invT)
            crotatedFcov=Rm*crotatedFcov*t(Rm)
            crotatedFmeans=refl*crotatedFmeans
            if(rescale==1){
              crotatedFcov=invsqrtFscale%*%crotatedFcov%*%invsqrtFscale
              crotatedFmeans=crotatedFmeans%*%invsqrtFscale
            }
            Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=crotatedFcov
            Output[[ind_nclust]]$groupspecific.factormeans[[g]]=crotatedFmeans
          }
        }
      }
      if(sum(design)>0 || nfactors==1){ # in case of CFA, only reflect
        loadings=Output[[ind_nclust]]$invariant.loadings
        nvar=nrow(loadings)
        strongload=loadings # for CFA, all non-zero loadings are considered strong loadings
        toreflect=colSums(strongload>0)<colSums(strongload<0)
        refl=matrix(1,nrow=1,ncol=nfactors)
        refl[toreflect]=-1
        reflloadings=matrix(refl,nrow=nvar,ncol=nfactors,byrow=TRUE)*loadings
        if(rescale==1){
          markload=reflloadings*markers
          markload=markload[markload!=0]
          sqrtFscale=diag(1/markload)
          invsqrtFscale=diag(markload)
          reflloadings=reflloadings%*%sqrtFscale
        }
        Output[[ind_nclust]]$invariant.loadings=reflloadings

        Rm=matrix(refl,nrow=nfactors,ncol=nfactors,byrow=TRUE)
        for(g in 1:ngroup){ # reflect all corresponding sets of factor (co)variances and factor means
          Fcov=Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]
          Fcov=Rm*Fcov*t(Rm)
          if(rescale==1){
            Fcov=invsqrtFscale%*%Fcov%*%invsqrtFscale
          }
          Output[[ind_nclust]]$groupspecific.factorcovariances[[g]]=Fcov
          if(is.element("intercepts",cluster.spec)){
            for(k in 1:nclust){
              Fmeans=Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]
              Fmeans=refl*Fmeans
              if(rescale==1){
                Fmeans=Fmeans%*%invsqrtFscale
              }
              Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]]=Fmeans
            }
          } else {
            Fmeans=Output[[ind_nclust]]$groupspecific.factormeans[[g]]
            Fmeans=refl*Fmeans
            Output[[ind_nclust]]$groupspecific.factormeans[[g]]=Fmeans
          }
        }
      }

      prefix="Factor"
      suffix=seq(1,nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      colnames(Output[[ind_nclust]]$invariant.loadings)<-factorlabels
      for(g in 1:ngroup){
        colnames(Output[[ind_nclust]]$groupspecific.factorcovariances[[g]])<-factorlabels
        rownames(Output[[ind_nclust]]$groupspecific.factorcovariances[[g]])<-factorlabels
        if(is.element("intercepts",cluster.spec)){
          for(k in 1:nclust){
            colnames(Output[[ind_nclust]]$group.and.clusterspecific.factormeans[[g,k]])<-factorlabels
          }
        } else {
          colnames(Output[[ind_nclust]]$groupspecific.factormeans[[g]])<-factorlabels
        }
      }
    }
  }

  if(length(OutputObject)==2){
    OutputObject$MMGFAsolutions=Output
  } else {
    OutputObject=Output
  }
  return(OutputObject)
}
