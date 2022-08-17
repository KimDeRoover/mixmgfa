#' Mixture Multigroup Factor Analysis
# -------------------------------------
#' @description
#' Perform mixture multigroup factor analyses (MMG-FA) with multiple numbers of clusters. MMG-FA aims to find clusters of groups with equivalent measurement model parameters (i.e., measurement invariance).
#'
# Written by Kim De Roover
# Please cite publications: https://doi.apa.org/doi/10.1037/met0000355
# and https://www.tandfonline.com/doi/full/10.1080/10705511.2020.1866577
# See these papers for recommendations on sample size, number of starts, model selection and how to move on from the results of these analyses.

# INPUT:
#' @param data A list consisting of "$covariances" (a vertically concatenated matrix or list of group-specific (co)variance matrices) and "$means" (a matrix with rows = group-specific means); or a matrix containing the vertically concatenated raw data for all groups.
#' Note: In case of raw data input without specifying N_gs, the first column of the data should contain group IDs. The remaining variables are then factor-analyzed.
#' @param N_gs Vector with number of subjects (sample size) for each group (in the same order as they appear in the data).
#' If left unspecified in case of raw data input, this vector is derived from the first column of the data matrix.
#' If left unspecified in case of covariance matrix & means input, a warning is issued.
#' @param nfactors Number of factors.
#' @param cluster.spec Measurement parameters you want to cluster the groups on; "loadings", "intercepts", "residuals", c("loadings","intercepts"), c("intercepts","residuals"), or c("loadings","intercepts","residuals").
#' Note: cluster.spec = "intercepts" and cluster.spec = c("intercepts","residuals") impose invariant loadings across all groups, cluster.spec = "residuals" also imposes invariant intercepts across all groups.
#' @param nsclust Vector of length two, indicating the minimal and maximal number of clusters (it is recommended to set the minimal number to one).
#' @param maxiter Maximum number of iterations used in each MMG-FA analysis. Increase in case of non-convergence.
#' @param nruns Number of (preselected) random starts (important for avoiding local maxima in case of few groups and/or small groups).
#' @param design For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).
#'          (Using different design matrices for different clusters is currently not supported.)
#' @param rotation Rotation criterion to use in case of EFA; currently either "oblimin", "varimax" or "target" (i.e., semi-specified oblique Procrustes rotation), whereas 0 = no rotation. (Note: The GPArotation package is loaded or installed for rotation.)
#' @param preselect Percentage of best starts taken in pre-selection of initial partitions (for huge datasets, increase to speed up multistart procedure).
#' @param targetT Target matrix to use when rotation = "target". Note that the same target is used for all clusters. For using cluster-specific target matrices, use rotation = 0 and rotate the cluster-specific factors afterwards with the GPFobl function of GPArotation.
#' @param targetW Weights to be used when rotation = "target". You can set an entire row to zero to make sure that the simple structure is not messed up by a 'complex variable' (i.e., with strong loadings for multiple factors). Set all weights equal to '1' if you prefer fully specified target rotation. When left unspecified while rotation = target, the default is a targetW where the zeros in the target get a weight of '1' and the non-zeros in the target get a weight of '0'. If this results in all zero weights or too many zeros for the rotation to be identified, a fully specified target rotation (all weights equal to '1') is used instead.


# OUTPUT:
#' @return Output object (list) with:
#'
#' $overview = overview of fitted MMG-FA solutions with loglikelihood (loglik), number of parameters (nrpars), BIC_N (using total number of observations as sample size), BIC_G (using number of groups as sample size),
#' CHull screeratios (NA = scree ratio could not be computed due to solution being first, last or not on the hull), convergence (1 = converged) and number of activated constraints on the unique variances.
#'
#' $MMGFAsolutions = list of MMG-FA solutions with different numbers of clusters.
#'          Access parameter values of solution with a specific number of clusters as, for example, OutputObject$MMGFAsolutions$"2.clusters".
#'
#' @references
#' De Roover, K., Vermunt, J. K., & Ceulemans, E. (2022). Mixture multigroup factor analysis for unraveling factor loading noninvariance across many groups. Psychological Methods, 27(3), 281-306.
#'
#' De Roover, K. (2021). Finding clusters of groups with measurement invariance: Unraveling intercept non-invariance with mixture multigroup factor analysis. Structural Equation Modeling: A Multidisciplinary Journal, 28(5), 663-683.
#'

#' @export
mixmgfa <- function(data,N_gs=c(),nfactors=1, cluster.spec = c("loadings","intercepts","residuals"),nsclust = c(1,5),maxiter = 5000,nruns = 25,design=0,rotation=0,preselect = 10,targetT=0,targetW=0){
  if(rotation!=0){
    if(rotation=="varimax" || rotation=="Varimax" || rotation=="VARIMAX"){
      rotation="varimax"
    } else if(rotation=="oblimin" || rotation=="Oblimin" || rotation=="OBLIMIN"){
      rotation="oblimin"
    } else if(rotation=="target" || rotation=="Target" || rotation=="TARGET"){
      rotation="target"
    } else {
      cat("Note: You seem to have misspelled your rotation criterion. Oblimin rotation will be performed.")
      cat("\n")
      rotation="oblimin"
    }
  }

  if(rotation!=0 & nfactors==1){
    rotation=0
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
          targetW=matrix(1,nvar,nfactors)
          cat("Note: targetW is unspecified. A valid weight matrix could not be derived from targetT, so fully specified target rotation will be performed.")

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

  if(max(design)!=0){
    if(ncol(design)!=nfactors){
      stop("The specified design matrix does not match the specified number of factors.")
    }
  }


  options(warn=-1)

  if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){ # input is raw data matrix
    if (is.null(N_gs)){ # create N_gs if not given as input
      T=table(data[,1]) #N_gs<-tabulate(data[,1])
      N_gs=as.numeric(T)
      grouplabels=names(T)
      N_gs=N_gs[N_gs!=0]
      N_gs=as.matrix(N_gs)
      checkfirstcolumn=0
      firstcolID=TRUE
    } else { # if N_gs is given as input, check whether first column of data contains group ID
      checkfirstcolumn=1
      firstcolID=FALSE
    }
    ngroup <- length(N_gs)
    if(nrow(N_gs)!=ngroup || is.null(nrow(N_gs))){ # make sure N_gs is a column vector
      N_gs_colvec=matrix(0,ngroup,1)
      for (g in 1:ngroup){
        N_gs_colvec[g,]=N_gs[g]
      }
      N_gs <- N_gs_colvec
    }
    if (checkfirstcolumn==1){
      T<-table(data[,1])
      checkcol1=as.numeric(T)
      checkcol1=checkcol1[checkcol1!=0]
      checkcol1=matrix(checkcol1,nrow=length(checkcol1),ncol=1)
      firstcolID=all.equal(N_gs,checkcol1)
      firstcolID=(firstcolID[1]==TRUE)
      # if(firstcolID){
      #   grouplabels=names(T)
      # }
    }

    if (firstcolID){
      nvar <- ncol(data)-1
      Xsup=as.matrix(data[,2:(nvar+1)])
      grouplabels=unique(data[,1])
      varlabels=colnames(data[,2:(nvar+1)])
    } else {
      nvar <- ncol(data)
      Xsup=as.matrix(data)
      grouplabels=rownames(N_gs)
      varlabels=colnames(data)
    }

    Ncum <- matrix(0,ngroup,2)
    Ncum[1,1]=1
    Ncum[1,2]=N_gs[1]
    for(g in 2:ngroup){
      Ncum[g,1]=sum(N_gs[1:(g-1)])+1 # Ncum[g,1]: first row in Xsup for group g
      Ncum[g,2]=sum(N_gs[1:g]) # Ncum[g,2]: last row in Xsup for group g
    }

    # compute group-specific observed means, covariances and mean cross-products
    mean_gs <- matrix(0,ngroup,nvar)
    obsS_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      mean_g <- colMeans(X_g)
      mean_gs[g,] <- mean_g
      # X_gc <- sweep(X_g,2,mean_g,check.margin = FALSE)
      CP_g <- (1/N_gs[g])*crossprod(X_g,X_g) #(t(X_g)%*%X_g)
      obsS_gs[[g]] <- CP_g-tcrossprod(mean_g,mean_g) #mean_g%*%t(mean_g)
    }
    if(length(cluster.spec)==1 & cluster.spec=="loadings"){
      data2=obsS_gs
    } else {
      data2=list(covariances=obsS_gs,means=mean_gs)
    }
  } else { # input is list or concatenation of covariance matrices and mean vectors
    if (is.null(N_gs)){
      stop("WARNING: You seem to have used covariances and means as input, without specifying the vector of group-specific sample sizes N_gs.")
    }
    grouplabels=rownames(N_gs)
    obsS_gs=data$covariances
    mean_gs=data$means
    if (is.list(obsS_gs)==FALSE){ # concatenation of covariance matrices should be turned into a list
      nvar=ncol(obsS_gs)
      varlabels=colnames(obsS_gs)
      obsS_gslist <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        obsS_gslist[[g]] <- obsS_gs[nvar*(g-1)+1:nvar*g,]
      }
      obsS_gs=obsS_gslist
    } else {
      nvar=ncol(obsS_gs[[1]])
      varlabels=colnames(obsS_gs[[1]])
      if(is.null(grouplabels)){
        grouplabels=names(obsS_gs)
      }
      if(is.null(grouplabels)){
        grouplabels=rownames(obsS_gs)
      }
    }
    if (is.list(mean_gs)){ # list of mean vectors should be turned into a matrix
      mean_gs=matrix(unlist(mean_gs),ncol=nvar,nrow=ngroup,byrow = TRUE)
    } else {
      mean_gs=matrix(mean_gs,ncol=nvar,nrow=ngroup,byrow = TRUE)
    }
    if(length(cluster.spec)==1 & cluster.spec=="loadings"){
      data2=obsS_gs
    } else {
      data2=list(covariances=obsS_gs,means=mean_gs)
    }
  }

  N=sum(N_gs)
  G=length(N_gs)
  if(is.null(grouplabels)){
    grouplabels=1:G
  }
  if(is.null(varlabels)){
    varlabels=1:nvar
  }

  if(rotation=="target" && max(targetT)!=0){
    if(nrow(targetT)!=nvar){
      stop("The specified targetT matrix does not match the number of variables.")
    }
  }
  if(max(design)!=0){
    if(nrow(design)!=nvar){
      stop("The specified design matrix does not match the number of variables.")
    }
  }

  overview=matrix(0,nsclust[2]-nsclust[1]+1,7)
  MMGFAsolutions=matrix(list(NA),nrow = 1, ncol=nsclust[2]-nsclust[1]+1)
  if (length(cluster.spec)==1 && cluster.spec[1]=="loadings"){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        #print(paste("Fitting MMG-FA with",nclust,"cluster ..."),quote=FALSE)
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_loadings(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        #print(paste("Fitting MMG-FA with",nclust,"clusters ..."),quote=FALSE)
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if (nclust==G){
          output_nclust <- mixmgfa_loadings(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else {
          output_nclust <- mixmgfa_loadings(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
          if(rotation!="target"){ # reflect and permute in case of rotated EFA (but not target rotation)
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
          for(g in 1:G){ # counter-rotate all corresponding sets of factor (co)variances
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
          for(g in 1:G){ # reflect all corresponding sets of factor (co)variances
            Fcov=output_nclust$Phi_gks[[g,k]]
            Fcov=Rm*Fcov*t(Rm)
            output_nclust$Phi_gks[[g,k]]=Fcov
          }
        }
      }
      prefix="Cluster"
      suffix=seq(1:nclust)
      colnames(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Phi_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      if(is.null(grouplabels)==FALSE){
        rownames(output_nclust$Psi_gs)<-grouplabels
        rownames(output_nclust$z_gks)<-grouplabels
        rownames(output_nclust$mu_gs)<-grouplabels
        rownames(output_nclust$Phi_gks)<-grouplabels
      }
      if(is.null(varlabels)==FALSE){
        colnames(output_nclust$mu_gs)<-varlabels
        for(k in 1:nclust){
          rownames(output_nclust$Lambda_ks[[k]])<-varlabels
        }
        for(g in 1:G){
          colnames(output_nclust$Psi_gs[[g]])<-varlabels
          rownames(output_nclust$Psi_gs[[g]])<-varlabels
        }
      }
      prefix="Factor"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      for(k in 1:nclust){
        colnames(output_nclust$Lambda_ks[[k]])<-factorlabels
        for(g in 1:G){
          colnames(output_nclust$Phi_gks[[g,k]])<-factorlabels
          rownames(output_nclust$Phi_gks[[g,k]])<-factorlabels
        }
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,groupspecific.uniquevariances=output_nclust$Psi_gs,groupspecific.means=mean_gs)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
  if (length(cluster.spec)==1 && cluster.spec=="intercepts"){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_intercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if(nclust==G){
          output_nclust <- mixmgfa_intercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else{
          output_nclust <- mixmgfa_intercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
      suffix=seq(1:nclust)
      rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      if(is.null(grouplabels)==FALSE){
        rownames(output_nclust$Psi_gs)<-grouplabels
        rownames(output_nclust$z_gks)<-grouplabels
        rownames(output_nclust$alpha_gks)<-grouplabels
        rownames(output_nclust$Phi_gs)<-grouplabels
      }
      if(is.null(varlabels)==FALSE){
        colnames(output_nclust$tau_ks)<-varlabels
        rownames(output_nclust$Lambda)<-varlabels
        for(g in 1:G){
          colnames(output_nclust$Psi_gs[[g]])<-varlabels
          rownames(output_nclust$Psi_gs[[g]])<-varlabels
        }
      }
      prefix="Factor"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Lambda)<-factorlabels
      for(g in 1:G){
        colnames(output_nclust$Phi_gs[[g]])<-factorlabels
        rownames(output_nclust$Phi_gs[[g]])<-factorlabels
        for(k in 1:nclust){
          colnames(output_nclust$alpha_gks[[g,k]])<-factorlabels
        }
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,invariant.loadings=output_nclust$Lambda,groupspecific.factorcovariances=output_nclust$Phi_gs,groupspecific.uniquevariances=output_nclust$Psi_gs,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
  if (length(cluster.spec)==1 && cluster.spec=="residuals"){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_residuals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if(nclust==G){
          output_nclust <- mixmgfa_residuals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else{
          output_nclust <- mixmgfa_residuals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
        for(g in 1:G){ # counter-rotate all corresponding sets of factor (co)variances
          crotatedFcov=invT%*%output_nclust$Phi_gs[[g]]%*%t(invT)
          crotatedFmeans=output_nclust$alpha_gs[[g]]%*%t(invT)
          if(rotation!="target"){
            crotatedFcov=Rm*crotatedFcov*t(Rm)
            crotatedFmeans=refl*crotatedFmeans
          }
          output_nclust$Phi_gs[[g]]=crotatedFcov
          output_nclust$alpha_gs[[g]]=crotatedFmeans
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
          Fmeans=output_nclust$alpha_gs[[g]]
          Fcov=Rm*Fcov*t(Rm)
          Fmeans=refl*Fmeans
          output_nclust$Phi_gs[[g]]=Fcov
          output_nclust$alpha_gs[[g]]=Fmeans
        }
      }

      prefix="Cluster"
      suffix=seq(1:nclust)
      colnames(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      if(is.null(grouplabels)==FALSE){
        rownames(output_nclust$z_gks)<-grouplabels
        rownames(output_nclust$alpha_gs)<-grouplabels
        rownames(output_nclust$Phi_gs)<-grouplabels
      }
      if(is.null(varlabels)==FALSE){
        colnames(output_nclust$tau)<-varlabels
        rownames(output_nclust$Lambda)<-varlabels
        for(k in 1:nclust){
          colnames(output_nclust$Psi_ks[[k]])<-varlabels
          rownames(output_nclust$Psi_ks[[k]])<-varlabels
        }
      }
      prefix="Factor"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Lambda)<-factorlabels
      for(g in 1:G){
        colnames(output_nclust$Phi_gs[[g]])<-factorlabels
        rownames(output_nclust$Phi_gs[[g]])<-factorlabels
        colnames(output_nclust$alpha_gs[[g]])<-factorlabels
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,invariant.loadings=output_nclust$Lambda,groupspecific.factorcovariances=output_nclust$Phi_gs,clusterspecific.uniquevariances=output_nclust$Psi_ks,invariant.intercepts=output_nclust$tau,groupspecific.factormeans=output_nclust$alpha_gs)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
  if (is.element("loadings",cluster.spec) & is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec)==FALSE){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_loadingsintercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if(nclust==G){
          output_nclust <- mixmgfa_loadingsintercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else{
          output_nclust <- mixmgfa_loadingsintercepts(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
            crotatedFmeans=output_nclust$alpha_gks[[g,k]]%*%t(invT)
            if(k>1 && rotation!="target"){
              crotatedFcov=crotatedFcov[perm,perm]
              crotatedFmeans=crotatedFmeans[,perm]
            }
            if(rotation!="target"){
              crotatedFcov=Rm*crotatedFcov*t(Rm)
              crotatedFmeans=refl*crotatedFmeans
            }
            output_nclust$Phi_gks[[g,k]]=crotatedFcov
            output_nclust$alpha_gks[[g,k]]=crotatedFmeans
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
            Fmeans=output_nclust$alpha_gks[[g,k]]
            Fcov=Rm*Fcov*t(Rm)
            Fmeans=refl*Fmeans
            output_nclust$Phi_gks[[g,k]]=Fcov
            output_nclust$alpha_gks[[g,k]]=Fmeans
          }
        }
      }
      prefix="Cluster"
      suffix=seq(1:nclust)
      colnames(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Phi_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      if(is.null(grouplabels)==FALSE){
        rownames(output_nclust$Psi_gs)<-grouplabels
        rownames(output_nclust$z_gks)<-grouplabels
        rownames(output_nclust$alpha_gks)<-grouplabels
        rownames(output_nclust$Phi_gks)<-grouplabels
      }
      if(is.null(varlabels)==FALSE){
        colnames(output_nclust$tau_ks)<-varlabels
        for(k in 1:nclust){
          rownames(output_nclust$Lambda_ks[[k]])<-varlabels
        }
        for(g in 1:G){
          colnames(output_nclust$Psi_gs[[g]])<-varlabels
          rownames(output_nclust$Psi_gs[[g]])<-varlabels
        }
      }
      prefix="Factor"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      for(k in 1:nclust){
        colnames(output_nclust$Lambda_ks[[k]])<-factorlabels
        for(g in 1:G){
          colnames(output_nclust$alpha_gks[[g,k]])<-factorlabels
          colnames(output_nclust$Phi_gks[[g,k]])<-factorlabels
          rownames(output_nclust$Phi_gks[[g,k]])<-factorlabels
        }
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,groupspecific.uniquevariances=output_nclust$Psi_gs,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
  if (is.element("loadings",cluster.spec) & is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec)){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_loadingsinterceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if(nclust==G){
          output_nclust <- mixmgfa_loadingsinterceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else{
          output_nclust <- mixmgfa_loadingsinterceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
            crotatedFmeans=output_nclust$alpha_gks[[g,k]]%*%t(invT)
            if(k>1 && rotation!="target"){
              crotatedFcov=crotatedFcov[perm,perm]
              crotatedFmeans=crotatedFmeans[,perm]
            }
            if(rotation!="target"){
              crotatedFcov=Rm*crotatedFcov*t(Rm)
              crotatedFmeans=refl*crotatedFmeans
            }
            output_nclust$Phi_gks[[g,k]]=crotatedFcov
            output_nclust$alpha_gks[[g,k]]=crotatedFmeans
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
            Fmeans=output_nclust$alpha_gks[[g,k]]
            Fcov=Rm*Fcov*t(Rm)
            Fmeans=refl*Fmeans
            output_nclust$Phi_gks[[g,k]]=Fcov
            output_nclust$alpha_gks[[g,k]]=Fmeans
          }
        }
      }
      prefix="Cluster"
      suffix=seq(1:nclust)
      colnames(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$Lambda_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$Psi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      rownames(output_nclust$tau_ks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$alpha_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Phi_gks)<-noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$z_gks)<-noquote(paste(prefix,suffix,sep="_"))
      names(output_nclust$pi_ks)<-noquote(paste(prefix,suffix,sep="_"))
      if(!is.null(grouplabels)){
        rownames(output_nclust$z_gks)<-grouplabels
        rownames(output_nclust$alpha_gks)<-grouplabels
        rownames(output_nclust$Phi_gks)<-grouplabels
      }
      if(!is.null(varlabels)){
        colnames(output_nclust$tau_ks)<-varlabels
        for(k in 1:nclust){
          rownames(output_nclust$Lambda_ks[[k]])<-varlabels
          colnames(output_nclust$Psi_ks[[k]])<-varlabels
          rownames(output_nclust$Psi_ks[[k]])<-varlabels
        }
      }
      prefix="Factor"
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      for(k in 1:nclust){
        colnames(output_nclust$Lambda_ks[[k]])<-factorlabels
        for(g in 1:G){
           colnames(output_nclust$alpha_gks[[g,k]])<-factorlabels
           colnames(output_nclust$Phi_gks[[g,k]])<-factorlabels
           rownames(output_nclust$Phi_gks[[g,k]])<-factorlabels
        }
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,clusterspecific.loadings=output_nclust$Lambda_ks,group.and.clusterspecific.factorcovariances=output_nclust$Phi_gks,clusterspecific.uniquevariances=output_nclust$Psi_ks,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
  if (is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec) & is.element("loadings",cluster.spec)==FALSE){
    for(nclust in nsclust[1]:nsclust[2]){
      if(nclust==1){
        cat(paste("Fitting MMG-FA with",nclust,"cluster ..."))
        cat("\n")
        output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
      }
      else {
        cat(paste("Fitting MMG-FA with",nclust,"clusters ..."))
        cat("\n")
        if(nclust==G){
          output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = 1,design=design, preselect=preselect)
        }
        else{
          output_nclust <- mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,maxiter = maxiter,start = 1,nruns = nruns,design=design, preselect=preselect)
        }
      }
      loglik=output_nclust$bestloglik
      nrpars=output_nclust$nrpars
      convergence=output_nclust$convergence>0
      overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(G),convergence,output_nclust$nractivatedconstraints)
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
      suffix=seq(1:nclust)
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
      suffix=seq(1:nfactors)
      factorlabels=noquote(paste(prefix,suffix,sep="_"))
      colnames(output_nclust$Lambda)<-factorlabels
      for(g in 1:G){
        colnames(output_nclust$Phi_gs[[g]])<-factorlabels
        rownames(output_nclust$Phi_gs[[g]])<-factorlabels
        for(k in 1:nclust){
          colnames(output_nclust$alpha_gks[[g,k]])<-factorlabels
        }
      }
      output_nclust2<-list(clustermemberships=output_nclust$z_gks,clusterproportions=output_nclust$pi_ks,invariant.loadings=output_nclust$Lambda,groupspecific.factorcovariances=output_nclust$Phi_gs,clusterspecific.uniquevariances=output_nclust$Psi_ks,clusterspecific.intercepts=output_nclust$tau_ks,group.and.clusterspecific.factormeans=output_nclust$alpha_gks)

      MMGFAsolutions[[nclust-nsclust[1]+1]]=output_nclust2
    }
  }
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
    convexhull=c(overview[1,1],convexhull,overview[nrows,1])
    nrhull=length(convexhull)
    change=0
    if(nrhull<nrows){
      change=1
    }
    while(nrhull>2 && change==1){ # check again whether intermediate points are on the convex hull
      nsclusthull=overview[convexhull,1]
      change=0
      for(nclust in 2:(nrhull-1)){
        if(!identical(convexhull[(nclust-1):(nclust+1)],c(convexhull[nclust]-1,convexhull[nclust],convexhull[nclust]+1))){
          LL_nclust=overview[convexhull[nclust]-nsclust[1]+1,2]
          npar_nclust=overview[convexhull[nclust]-nsclust[1]+1,3]
          LL_nclustmin1=overview[convexhull[nclust-1]-nsclust[1]+1,2]
          npar_nclustmin1=overview[convexhull[nclust-1]-nsclust[1]+1,3]
          LL_nclustplus1=overview[convexhull[nclust+1]-nsclust[1]+1,2]
          npar_nclustplus1=overview[convexhull[nclust+1]-nsclust[1]+1,3]
          # determine whether intermediate point is part of the convex hull
          slope=(LL_nclustplus1-LL_nclustmin1)/(npar_nclustplus1-npar_nclustmin1)
          point_line=LL_nclustmin1+slope*(npar_nclust-npar_nclustmin1)
          if(LL_nclust>=(point_line-.01)){
            # when the subset of three points spans across a point not on the hull, this is the corrected scree ratio (comparing with the previous and next point ON THE HULL)
            screeratios[convexhull[nclust]-nsclust[1]+1]=((LL_nclust-LL_nclustmin1)/(npar_nclust-npar_nclustmin1))/((LL_nclustplus1-LL_nclust)/(npar_nclustplus1-npar_nclust))
          } else {
            screeratios[convexhull[nclust]-nsclust[1]+1]=NA
            change=1
          }
        }
      }
      convexhull=which(!is.na(screeratios))
      convexhull=c(overview[1,1],convexhull,overview[nrows,1])
      nrhull=length(convexhull)
    }
    overview=cbind(overview[,1:5],screeratios,overview[,6:7])
    }
    prefix=seq(nsclust[1]:nsclust[2])
    suffix="clusters"
    sollistnames<-c(paste(prefix,suffix,sep="."))
    names(MMGFAsolutions)<-noquote(sollistnames)
    #overview=as.data.frame(overview,row.names=FALSE)
    if(nrows>2){
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","screeratios","convergence","nr.activated.constraints")
    } else {
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","convergence","nr.activated.constraints")
    }

    # windows(8,6)
    # par(mfrow=c(2,1))
    # par(mar=c(4.1,4.1,4.1,2.1))
    # plot(overview[,1],overview[,5],type="b",xlab="number of clusters",ylab="BIC_G")
    # plot(overview[,1],overview[,2],type="b",xlab="number of clusters",ylab="loglik")
    # if(length(cluster.spec)==1){
    #   mtext(paste("Mixture multigroup factor analyses with cluster-specific",cluster.spec), side = 3, line = -2, outer = TRUE)
    # }
    # if(length(cluster.spec)==2){
    #   mtext(paste("Mixture multigroup factor analyses with cluster-specific",cluster.spec[1],"and",cluster.spec[2]), side = 3, line = -2, outer = TRUE)
    # }
    # if(length(cluster.spec)==3){
    #   mtext(paste("Mixture multigroup factor analyses with cluster-specific",cluster.spec[1],",",cluster.spec[2],"and",cluster.spec[3]), side = 3, line = -2, outer = TRUE)
    # }

    cat("\n")
    print(overview)
    cat("\n")
    cat("Choose the best number of clusters ('K_best') based on the BIC_G and CHull scree ratios and the plots. For plots, use 'plot(OutputObject$overview)'.")
    cat("\n")
    cat("Based on the BIC_G, look for the number of clusters that minimizes the BIC_G or that corresponds to an elbow point in the BIC_G plot (after which the decrease with additional clusters levels off).")
    cat("\n")
    cat("Based on the CHull (scree ratios AND plot), look for the number of clusters that has the maximal scree ratio AND check whether this corresponds to at least a mild elbow point in the lower plot.")
    cat("\n")
    cat("Access the corresponding cluster memberships and parameter estimates by using OutputObject$MMGFAsolutions[[K_best]]$clustermemberships and, for example, OutputObject$MMGFAsolutions[[K_best]]$clusterspecific.loadings.")
    cat("\n")
    cat("The parameter sets are further subdivided in group- and/or cluster-specific parameter sets.")
    cat("\n")
    class(MMGFAsolutions)<-"mixmgfa"
    class(overview)<-"mixmgfaMS"
    output<-list(overview=overview,MMGFAsolutions=MMGFAsolutions)
    #class(output)<-"mixmgfa"

  } else {
    cat("You seem to have made a typo in 'cluster.spec' (make sure you use lowercase letters) or have requested a model that is not supported.")
    cat("\n")
    cat("Please try again.")
    cat("\n")
  }
}
