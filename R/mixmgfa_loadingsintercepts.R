#' Mixture multigroup factor analysis for loadings and intercepts (non-)invariance
# ----------------------------------------------------------------
# Code written by Kim De Roover

#' @description
#' Finds clusters of groups based on their loadings and intercepts, given a user-specified number of clusters.
#' Residual variances are group-specific. Factor (co)variances and factor means are group- and cluster-specific.

# INPUT:
#' @param data A matrix of vertically concatenated, group-specific (co)variance matrices; or a matrix containing the vertically concatenated raw data for all groups (first column should not contain group ID).
#' @param N_gs Vector containing the sample size for each group (in the same order as they appear in the data).
#' @param nclust User-specified number of clusters.
#' @param nfactors User-specified number of factors.
#' @param maxiter Maximum number of iterations. Increase in case of non-convergence.
#' @param start Type of start (start = 1: pre-selected random starts, start = 2: start from a user-specified startpartition).
#' @param nruns Number of starts (based on pre-selected random partitions when start = 1).
#' @param preselect Percentage of best starts taken in pre-selection (for huge data, increase to speed up startprocedure).
#' @param design For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).
#' @param startpartition Partition of groups (vector) to start from (use with start = 2 and nruns = 1).

# OUTPUT:
#' @return Output object with:
#'
#' $z_gks = cluster memberships of groups (posterior classification probabilities)
#'
#' $pi_ks= mixing proportions (prior classification probabilities)
#'
#' $Lambda_ks = cluster-specific loadings, access loadings of cluster k via Lambda_ks[[k]]
#'
#' $Psi_gs = group-specific unique variances, access loadings of group g via Psi_gs[[g]]
#'
#' $Phi_gks = group- and cluster-specific factor (co)variances, access (co)variances of group g in cluster k via Phi_gks[[g,k]]
#'
#' $tau_ks = cluster-specific intercepts, access intercepts of cluster k via tau_ks[k,]
#'
#' $alpha_gks = group- and cluster-specific factor means, access factor means of group g in cluster k via alpha_gks[[g,k]]
#'
#' $bestloglik = loglikelihood of best start
#'
#' $logliks = loglikelihoods (first column) and number of activated constraints on the unique variances (second column) for all starts
#'
#' $nrpars = number of free parameters, to be used for model selection in combination with bestloglik
#'
#' $convergence = 2 if best start converged on loglikelihood, 1 if converged on parameter changes, 0 if not converged
#'
#' $nractivatedconstraints = number of constraints on the unique variances (across groups, for best start) to avoid unique variances approaching zero

#' @export
mixmgfa_loadingsintercepts <- function(data,N_gs,nclust,nfactors=1,maxiter = 5000,start = 1,nruns = 25,design = 0,preselect = 10,startpartition){

  ngroup <- length(N_gs)
  if(nrow(N_gs)!=ngroup || is.null(nrow(N_gs))){ # make sure N_gs is a column vector
    N_gs_colvec=matrix(0,ngroup,1)
    for (g in 1:ngroup){
      N_gs_colvec[g,]=N_gs[g]
    }
    N_gs <- N_gs_colvec
  }
  N <- sum(N_gs);
  IM <- diag(nclust)

  if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){ # input is raw data matrix
    nvar <- ncol(data)
    Xsup=as.matrix(data)
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
    CP_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      mean_g <- colMeans(X_g)
      mean_gs[g,] <- mean_g
      # X_gc <- sweep(X_g,2,mean_g,check.margin = FALSE)
      CP_g <- (1/N_gs[g])*crossprod(X_g,X_g)
      CP_gs[[g]] <- CP_g
      obsS_gs[[g]] <- CP_g-tcrossprod(mean_g,mean_g) #(1/N_gs[g])*(t(X_gc)%*%X_gc)
    }
  } else { # input is list or concatenation of covariance matrices and mean vectors
    obsS_gs=data$covariances
    mean_gs=data$means
    if (is.list(obsS_gs)==FALSE){ # concatenation of covariance matrices should be turned into a list
      nvar=ncol(obsS_gs)
      obsS_gslist <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        obsS_gslist[[g]] <- obsS_gs[nvar*(g-1)+1:nvar*g,]
      }
      obsS_gs=obsS_gslist
    } else { # input is a list of (1) a list of covariance matrices and (2) a list of mean vectors
      nvar=ncol(obsS_gs[[1]])
    }
    if (is.list(mean_gs)){ # list of mean vectors should be turned into a matrix
      mean_gs=matrix(unlist(mean_gs),ncol=nvar,nrow=ngroup,byrow = TRUE)
    } else {
      mean_gs=matrix(mean_gs,ncol=nvar,nrow=ngroup,byrow = TRUE)
    }

    # compute group-specific mean cross-products
    CP_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      mean_g <- mean_gs[g,]
      CP_g <- obsS_gs[[g]]+tcrossprod(mean_g,mean_g)
      CP_gs[[g]] <- CP_g
    }
  }

  if (sum(design)==0){ # if design is unspecified, EFA is used
    design=matrix(1,nvar,nfactors)
    EFA=1
  } else {
    EFA=0
  }


  if(start==1){
    if(nclust>1 && nclust<ngroup){
      # pre-selection of random partitions
      nrtrialstarts=min(nruns*round(100/preselect),stirlingnr2(ngroup,nclust)) # generate 'nruns'*(100/preselect) different random partitions
      randpartvecs=matrix(0,nrtrialstarts,ngroup);
      for (trialstart in 1:nrtrialstarts){
        aris=1;
        #iterrnd=0;
        while (sum(aris==1)>0){ # && iterrnd<5){
          cl=0;
          while(length(cl)<nclust){
            randpartvec <- sample(1:nclust,ngroup,replace=TRUE) # generate random partition
            cl=unique(randpartvec)
          }
          nstartprev=trialstart-1
          aris=matrix(0,nstartprev,1)
          if (nstartprev>0){
            for (r in 1:nstartprev){
              prevpartvec=randpartvecs[r,]
              aris[r]<-adjrandindex(prevpartvec,randpartvec)
            }
            #iterrnd=iterrnd+1;
          }
        }
        randpartvecs[trialstart,]=randpartvec
      }
      if (preselect<100){
        ODLLs_trialstarts=rep(0,nrtrialstarts,1)
        for (trialstart in 1:nrtrialstarts){ # select 'preselect'% (default 10%) best fitting random partitions
          randpartvec=randpartvecs[trialstart,];
          if(nclust>1){
            z_gks=IM[randpartvec,]
            pi_ks=(1/ngroup)*colSums(z_gks)
          }
          else {
            z_gks=t(randpartvec)
            pi_ks=1
          }
          N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
          N_ks=colSums(N_gks)

          # tau_ks <- matrix(0,nclust,nvar)
          # for(k in 1:nclust){
          #   Xsup_k=matrix(0,N_ks[k],nvar)
          #   for(g in 1:ngroup){
          #     if(randpartvec[g]==k){
          #       X=Xsup[Ncum[g,1]:Ncum[g,2],]
          #       Xsup_k[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=X
          #     }
          #   }
          #   tau_ks[k,] <- colMeans(Xsup_k)
          # }
          # # center data with initialized intercepts
          # Xsupcent <- matrix(0,N,nvar)
          # for(g in 1:ngroup){
          #   k=randpartvec[g]
          #   X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
          #   Xsupcent[Ncum[g,1]:Ncum[g,2],]=sweep(X_g,2,tau_ks[k,],check.margin = FALSE)
          # }

          tau_ks <- matrix(0,nclust,nvar)
          for(k in 1:nclust){
            for(g in 1:ngroup){
              if(randpartvec[g]==k){
                tau_ks[k,]=tau_ks[k,]+(N_gks[g,k]/N_ks[k])*mean_gs[g,]
              }
            }
          }

          Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
          uniq_ks <- matrix(0,nclust,nvar)
          for(k in 1:nclust){
            S_k=matrix(0,nvar,nvar)
            for(g in 1:ngroup){
              if(randpartvec[g]==k){
                mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
                mu_gk=matrix(tau_ks[k,],nrow=1,ncol=nvar)
                S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
                S_k=S_k+N_gs[g]*S_gk
              }
            }
            S_k=(1/N_ks[k])*S_k
            ed<-eigen(S_k, symmetric=TRUE, only.values = FALSE)
            val<-ed$values
            u<-ed$vectors
            totalerror=sum((val[-seq_len(nfactors)]))
            meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
            Uniq=rep(meanerror,nvar)
            lambda_k=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
            if (EFA==0){
              lambda_k <- procr(lambda_k,design)
              lambda_k=lambda_k*design # non-zero loadings should be indicated with '1' for this to work properly
            }
            Lambda_ks[[k]]=lambda_k
            uniq_ks[k,]=Uniq
          }



          Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
          Phi_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # initialize group- and cluster-specific factor covariances
          alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
          for(g in 1:ngroup){
            k=randpartvec[g]
            Psi_gs[[g]]=diag(uniq_ks[k,])
            if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
              X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
            }
            for(k2 in 1:nclust){
              Phi_gks[[g,k2]]=diag(nfactors)
              lambda_k=Lambda_ks[[k2]]
              if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
                #X_c=X_g-(matrix(tau_ks[k2,],ncol=nvar,nrow=N_gs[g],byrow=TRUE))
                X_c=sweep(X_g,2,tau_ks[k2,],check.margin = FALSE)
                udv=svd(X_c%*%lambda_k)
                U=udv$u
                V=udv$v
                Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
                Fvar=apply(Fscores,2,var)
                Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
                alpha_gks[[g,k2]] <- colMeans(Fscores)
              } else {
                tlambda_k=t(lambda_k)
                S_g=obsS_gs[[g]]
                tlambdak_Sg=tlambda_k%*%solve(S_g)
                alpha_gks[[g,k2]]=(mean_gs[g,]-tau_ks[k2,])%*%t(solve(((tlambdak_Sg)%*%lambda_k),(tlambdak_Sg)))
              }
            }
          }


          Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
          invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
          for(k in 1:nclust){
            lambda_k=Lambda_ks[[k]]
            tlambda_k=t(lambda_k)
            for(g in 1:ngroup){
              psi_g=Psi_gs[[g]]
              invPsi_g=diag(1/diag(psi_g))
              invPsi_g_lambda_k=invPsi_g%*%lambda_k
              phi_gk=Phi_gks[[g,k]]
              invPhi_gk=phi_gk # phi_gk is still identity matrix
              sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
              Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
              invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk+tlambda_k%*%invPsi_g_lambda_k)
              invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk_tLambdainvPsi_g_Lambda+t(invPhi_gk_tLambdainvPsi_g_Lambda))*(1/2)
              invSigma_gks[[g,k]]=invPsi_g-invPsi_g_lambda_k%*%solve(invPhi_gk_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda_k) # Woodbury identity
            }
          }

          # if (is.matrix(data)==FALSE & is.data.frame(data)==FALSE){
          #   # cycle 1 update factor means (see De Roover, 2021)
          #   for(g in 1:ngroup){
          #     for(k in 1:nclust){
          #       if(N_ks[k]>0){
          #         lambda_k=Lambda_ks[[k]]
          #         tlambda_k=t(lambda_k)
          #         invSigma_gk=invSigma_gks[[g,k]]
          #         alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tlambda_k%*%invSigma_gk)%*%lambda_k),(tlambda_k%*%invSigma_gk)))
          #       }
          #     }
          #   }
          # }

          S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
          S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
          alpha_tlambda_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
          for(g in 1:ngroup){
            S_g=matrix(0,nvar,nvar)
            mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
            for(k in 1:nclust){
              lambda_k=Lambda_ks[[k]]
              tlambda_k=t(lambda_k)
              alpha_tlambda_gks[[g,k]]=alpha_gks[[g,k]]%*%tlambda_k
              mu_gk=tau_ks[k,]+alpha_tlambda_gks[[g,k]] # model-implied mean if group g is in cluster k
              S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
              S_gks[[g,k]] <- S_gk
              S_g=S_g+N_gks[g,k]*S_gk
            }
            S_gs[[g]] <- (1/N_gs[g])*S_g
          }

          # compute Beta_gks and theta_gks
          Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
          Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
          meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
          #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
          for(g in 1:ngroup){
            for(k in 1:nclust){
              lambda_k=Lambda_ks[[k]]
              phi_gk=Phi_gks[[g,k]]
              invsigma_gk=invSigma_gks[[g,k]]
              beta_gk=phi_gk%*%crossprod(lambda_k,invsigma_gk)
              tbeta_gk=t(beta_gk)
              Beta_gks[[g,k]]=beta_gk;
              S_gk=S_gks[[g,k]]
              theta_gk=phi_gk-beta_gk%*%(lambda_k%*%phi_gk-S_gk%*%tbeta_gk)
              Theta_gks[[g,k]]=theta_gk
              meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-alpha_tlambda_gks[[g,k]])%*%tbeta_gk
              #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
            }
          }

          Output_Mstep <- mixmgfa_loadingsintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mean_gs,tau_ks,alpha_gks)
          Lambda_ks=Output_Mstep$Lambda_ks
          Psi_gs=Output_Mstep$Psi_gs
          Phi_gks=Output_Mstep$Phi_gks
          tau_ks=Output_Mstep$tau_ks
          alpha_gks=Output_Mstep$alpha_gks
          Sigma_gks=Output_Mstep$Sigma_gks
          invSigma_gks=Output_Mstep$invSigma_gks
          nractivatedconstraints=Output_Mstep$nractivatedconstraints

          S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
          S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
          for(g in 1:ngroup){
            S_g=matrix(0,nvar,nvar)
            mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
            for(k in 1:nclust){
              lambda_k=Lambda_ks[[k]]
              tlambda_k=t(lambda_k)
              mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tlambda_k # model-implied mean if group g is in cluster k
              S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
              S_gks[[g,k]] <- S_gk
              S_g=S_g+N_gks[g,k]*S_gk
            }
            S_gs[[g]] <- (1/N_gs[g])*S_g
          }

          # compute observed-data log-likelihood for start
          ODLL_trialstart=0;
          loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
          loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
          for(g in 1:ngroup){
            # X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
            # for(k in 1:nclust){
            #   logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
            #   invSigma_gk=invSigma_gks[[g,k]]
            #   lambda_k=Lambda_ks[[k]]
            #   Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(lambda_k),ncol=N_gs[g],nrow=nvar)) # centered data per group
            #   loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
            #   for (n in 1:N_gs[g]){
            #     Xc_n=Xc_gk[n, ,drop=FALSE]
            #     loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_gk,Xc_n))
            #   }
            #   loglik_gks[g,k]=loglik_gk
            #   loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
            # }
            for(k in 1:nclust){
              logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
              loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gks[[g,k]]*invSigma_gks[[g,k]])) # sum(S_gks[[g,k]]*invSigma_gks[[g,k]])=sum(diag(S_gks[[g,k]]%*%invSigma_gks[[g,k]]))
              loglik_gks[g,k]=loglik_gk
              loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
            }
            m_i=max(loglik_gksw[g,]);
            for(k in 1:nclust){
              loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
            }
            ODLL_trialstart=ODLL_trialstart+log(sum(loglik_gksw[g,]))+m_i;
          }
          ODLLs_trialstarts[trialstart]=ODLL_trialstart;
        }
      }
      if (nrtrialstarts>nruns){
        sortedstarts <- sort(ODLLs_trialstarts,decreasing=TRUE,index.return=TRUE);
        index_beststarts=sortedstarts$ix[1:nruns]
        randpartvecs=randpartvecs[index_beststarts,]
        if(nruns==1){
          randpartvecs=matrix(randpartvecs)
        }
      }
      else {
        nruns=nrtrialstarts
      }
    }
    else {
      nruns=1
      if(nclust==1){
        randpartvecs=matrix(1,1,ngroup)
      }
      if (nclust==ngroup){
        randpartvecs=matrix(1:ngroup,1,ngroup)
      }
    }
  }

  # start of loop of multiple starts
  convergence <- 1
  logliks <- matrix(0,nruns,2)
  for(run in 1:nruns){
    nractivatedconstraints <- 0
    if(start==1){
      if(nruns>1){
        randpartvec <- randpartvecs[run,]
      }
      else {
        randpartvec <- randpartvecs
      }
    }
    if(start==2){
      randpartvec <- startpartition
    }

    if(nclust>1){
      z_gks=IM[randpartvec,]
      pi_ks=(1/ngroup)*colSums(z_gks)
    }
    else {
      z_gks=t(randpartvec)
      pi_ks=1
    }
    N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=colSums(N_gks)

    tau_ks <- matrix(0,nclust,nvar)
    for(k in 1:nclust){
      for(g in 1:ngroup){
        if(randpartvec[g]==k){
          tau_ks[k,]=tau_ks[k,]+(N_gks[g,k]/N_ks[k])*mean_gs[g,]
        }
      }
    }

    Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
    uniq_ks <- matrix(0,nclust,nvar)
    for(k in 1:nclust){
      S_k=matrix(0,nvar,nvar)
      for(g in 1:ngroup){
        if(randpartvec[g]==k){
          mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
          mu_gk=matrix(tau_ks[k,],nrow=1,ncol=nvar)
          S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)-crossprod(mu_gk,mean_g)+crossprod(mu_gk,mu_gk)
          S_k=S_k+N_gs[g]*S_gk
        }
      }
      S_k=(1/N_ks[k])*S_k
      ed<-eigen(S_k, symmetric=TRUE, only.values = FALSE)
      val<-ed$values
      u<-ed$vectors
      totalerror=sum((val[-seq_len(nfactors)]))
      meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
      Uniq=rep(meanerror,nvar)
      lambda_k=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
      if (EFA==0){
        lambda_k <- procr(lambda_k,design)
        lambda_k=lambda_k*design # non-zero loadings should be indicated with '1' for this to work properly
      }
      Lambda_ks[[k]]=lambda_k
      uniq_ks[k,]=Uniq
    }



    Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
    Phi_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # initialize group- and cluster-specific factor covariances
    alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
    for(g in 1:ngroup){
      k=randpartvec[g]
      Psi_gs[[g]]=diag(uniq_ks[k,])
      if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
        X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      }
      for(k2 in 1:nclust){
        Phi_gks[[g,k2]]=diag(nfactors)
        lambda_k=Lambda_ks[[k2]]
        if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
          #X_c=X_g-(matrix(tau_ks[k2,],ncol=nvar,nrow=N_gs[g],byrow=TRUE))
          X_c=sweep(X_g,2,tau_ks[k2,],check.margin = FALSE)
          udv=svd(X_c%*%lambda_k)
          U=udv$u
          V=udv$v
          Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
          Fvar=apply(Fscores,2,var)
          Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
          alpha_gks[[g,k2]] <- colMeans(Fscores)
        } else {
          tlambda_k=t(lambda_k)
          S_g=obsS_gs[[g]]
          tlambdak_Sg=tlambda_k%*%solve(S_g)
          alpha_gks[[g,k2]]=(mean_gs[g,]-tau_ks[k2,])%*%t(solve(((tlambdak_Sg)%*%lambda_k),(tlambdak_Sg)))
        }
      }
    }


    Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      tlambda_k=t(lambda_k)
      for(g in 1:ngroup){
        psi_g=Psi_gs[[g]]
        invPsi_g=diag(1/diag(psi_g))
        invPsi_g_lambda_k=invPsi_g%*%lambda_k
        phi_gk=Phi_gks[[g,k]]
        invPhi_gk=solve(phi_gk)
        sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
        Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
        invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk+tlambda_k%*%invPsi_g_lambda_k)
        invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk_tLambdainvPsi_g_Lambda+t(invPhi_gk_tLambdainvPsi_g_Lambda))*(1/2)
        invSigma_gks[[g,k]]=invPsi_g-invPsi_g_lambda_k%*%solve(invPhi_gk_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda_k) # Woodbury identity
      }
    }

    # if (is.matrix(data)==FALSE & is.data.frame(data)==FALSE){
    #   # cycle 1 update factor means (see De Roover, 2021)
    #   for(g in 1:ngroup){
    #     for(k in 1:nclust){
    #       if(N_ks[k]>0){
    #         lambda_k=Lambda_ks[[k]]
    #         tlambda_k=t(lambda_k)
    #         invSigma_gk=invSigma_gks[[g,k]]
    #         alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tlambda_k%*%invSigma_gk)%*%lambda_k),(tlambda_k%*%invSigma_gk)))
    #       }
    #     }
    #   }
    # }


    S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        tlambda_k=t(lambda_k)
        mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tlambda_k # model-implied mean if group g is in cluster k
        S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
        S_gks[[g,k]] <- S_gk
        S_g=S_g+N_gks[g,k]*S_gk
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }

    # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
    for(g in 1:ngroup){
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gks[[g,k]]*invSigma_gks[[g,k]])) # sum(S_gks[[g,k]]*invSigma_gks[[g,k]])=sum(diag(S_gks[[g,k]]%*%invSigma_gks[[g,k]]))
        loglik_gks[g,k]=loglik_gk
      }
    }

    iter=0;
    conv1=1;
    conv2=1;
    ODLL=-Inf;
    pars=pi_ks
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      pars=c(pars,lambda_k[design==1])
    }
    pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
    pars=unlist(pars)
    while(min(conv1,conv2)>1e-4 && iter<100){
      prev_ODLL=ODLL
      prev_Lambda_ks=Lambda_ks
      prev_Psi_gs=Psi_gs
      prev_Phi_gks=Phi_gks
      prev_tau_ks=tau_ks
      prev_alpha_gks=alpha_gks
      prev_pars=pars
      iter=iter+1

      # **E-step**: compute the posterior classification probabilities
      if(nclust>1 && nclust<ngroup){
        z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
        pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
      }

      N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
      N_ks=colSums(N_gks)


      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      alpha_tlambda_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          alpha_tlambda_gks[[g,k]]=tcrossprod(alpha_gks[[g,k]],lambda_k) # used below
          S_g=S_g+N_gks[g,k]*S_gks[[g,k]]
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }

      # compute Beta_gks and theta_gks
      Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
      #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
      for(g in 1:ngroup){
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          phi_gk=Phi_gks[[g,k]]
          invsigma_gk=invSigma_gks[[g,k]]
          beta_gk=phi_gk%*%crossprod(lambda_k,invsigma_gk)
          tbeta_gk=t(beta_gk)
          Beta_gks[[g,k]]=beta_gk
          S_gk=S_gks[[g,k]]
          theta_gk=phi_gk-beta_gk%*%(lambda_k%*%phi_gk-S_gk%*%tbeta_gk)
          Theta_gks[[g,k]]=theta_gk
          meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-alpha_tlambda_gks[[g,k]])%*%tbeta_gk
          #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
        }
      }

      Output_Mstep <- mixmgfa_loadingsintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mean_gs,tau_ks,alpha_gks)
      Lambda_ks=Output_Mstep$Lambda_ks
      Psi_gs=Output_Mstep$Psi_gs
      Phi_gks=Output_Mstep$Phi_gks
      tau_ks=Output_Mstep$tau_ks
      alpha_gks=Output_Mstep$alpha_gks
      Sigma_gks=Output_Mstep$Sigma_gks
      invSigma_gks=Output_Mstep$invSigma_gks
      nractivatedconstraints=Output_Mstep$nractivatedconstraints

      S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          tlambda_k=t(lambda_k)
          mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tlambda_k # model-implied mean if group g is in cluster k
          S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
          S_gks[[g,k]] <- S_gk
          S_g=S_g+N_gks[g,k]*S_gk
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }

      # check on change in observed-data log-likelihood
      ODLL=0;
      loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
      loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
      for(g in 1:ngroup){
        for(k in 1:nclust){
          logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
          loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gks[[g,k]]*invSigma_gks[[g,k]])) # sum(S_gks[[g,k]]*invSigma_gks[[g,k]])=sum(diag(S_gks[[g,k]]%*%invSigma_gks[[g,k]]))
          loglik_gks[g,k]=loglik_gk
          loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
        }
        m_i=max(loglik_gksw[g,]);
        for(k in 1:nclust){
          loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
        }
        ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
      }


      pars=pi_ks
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        pars=c(pars,lambda_k[design==1])
      }
      pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
      pars=unlist(pars)
      ind=which(pars!=0 & prev_pars!=0)
      parsdiff <- pars[ind]-prev_pars[ind] #mapply("-",pars,prev_pars)
      parsdiffdiv <- parsdiff/prev_pars[ind] #mapply("/",parsdiff,prev_pars)
      parsabsdiffdiv <- abs(parsdiffdiv)#lapply(parsdiffdiv,function(x){abs(x)})
      conv1 <- sum(parsabsdiffdiv)#do.call(sum,parsabsdiffdiv)
      if(is.na(conv1)){
        conv1=1
      }

      conv2=ODLL-prev_ODLL

      if(ODLL-prev_ODLL<0){
       print(ODLL-prev_ODLL)
      }


    } # end while-loop till convergence
    logliks[run,]=c(ODLL,nractivatedconstraints);


    if (run==1) {
      bestz_gks=z_gks
      bestpi_ks=pi_ks
      bestLambda_ks=Lambda_ks
      bestPsi_gs=Psi_gs
      bestPhi_gks=Phi_gks
      besttau_ks=tau_ks
      bestalpha_gks=alpha_gks
      bestSigma_gks=Sigma_gks
      bestinvSigma_gks=invSigma_gks
      bestloglik=ODLL
      bestloglik_gks=loglik_gks
      bestiter=iter
      bestconv1=conv1
      bestconv2=conv2
    }
    else {
      if (ODLL>bestloglik){
        bestz_gks=z_gks
        bestpi_ks=pi_ks
        bestLambda_ks=Lambda_ks
        bestPsi_gs=Psi_gs
        bestPhi_gks=Phi_gks
        besttau_ks=tau_ks
        bestalpha_gks=alpha_gks
        bestSigma_gks=Sigma_gks
        bestinvSigma_gks=invSigma_gks
        bestloglik=ODLL
        bestloglik_gks=loglik_gks
        bestiter=iter
        bestconv1=conv1
        bestconv2=conv2
      }
    }


  } # end for-loop over multiple starts

  z_gks=bestz_gks
  pi_ks=bestpi_ks
  Lambda_ks=bestLambda_ks
  Psi_gs=bestPsi_gs
  Phi_gks=bestPhi_gks
  tau_ks=besttau_ks
  alpha_gks=bestalpha_gks
  Sigma_gks=bestSigma_gks
  invSigma_gks=bestinvSigma_gks
  ODLL=bestloglik
  loglik_gks=bestloglik_gks
  iter=bestiter
  conv1=bestconv1
  conv2=bestconv2

  pars=pi_ks
  for(k in 1:nclust){
    lambda_k=Lambda_ks[[k]]
    pars=c(pars,lambda_k[design==1])
  }
  pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
  pars=unlist(pars)
  while(min(conv1,conv2)>1e-6 && iter<maxiter+1){ # iterate till convergence for best start
    prev_ODLL=ODLL
    prev_Lambda_ks=Lambda_ks
    prev_Psi_gs=Psi_gs
    prev_Phi_gks=Phi_gks
    prev_tau_ks=tau_ks
    prev_alpha_gks=alpha_gks
    prev_pars=pars
    iter=iter+1

    # **E-step**: compute the posterior classification probabilities
    if(nclust>1 && nclust<ngroup){
      z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
      pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
    }

    N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=colSums(N_gks)


    if(iter==bestiter+1){
      S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      alpha_tlambda_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          tlambda_k=t(lambda_k)
          alpha_tlambda_gks[[g,k]]=alpha_gks[[g,k]]%*%tlambda_k
          mu_gk=tau_ks[k,]+alpha_tlambda_gks[[g,k]] # model-implied mean if group g is in cluster k
          S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
          S_gks[[g,k]] <- S_gk
          S_g=S_g+N_gks[g,k]*S_gk
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }
    } else {
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      alpha_tlambda_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        for(k in 1:nclust){
          alpha_tlambda_gks[[g,k]]=tcrossprod(alpha_gks[[g,k]],Lambda_ks[[k]]) # used below
          S_g=S_g+N_gks[g,k]*S_gks[[g,k]]
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }
    }


    # compute Beta_gks and theta_gks
    Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
    #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
    for(g in 1:ngroup){
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        phi_gk=Phi_gks[[g,k]]
        invsigma_gk=invSigma_gks[[g,k]]
        beta_gk=phi_gk%*%crossprod(lambda_k,invsigma_gk)
        tbeta_gk=t(beta_gk)
        Beta_gks[[g,k]]=beta_gk
        S_gk=S_gks[[g,k]]
        theta_gk=phi_gk-beta_gk%*%(lambda_k%*%phi_gk-S_gk%*%tbeta_gk)
        Theta_gks[[g,k]]=theta_gk
        meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-alpha_tlambda_gks[[g,k]])%*%tbeta_gk
        #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
      }
    }

    Output_Mstep <- mixmgfa_loadingsintercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mean_gs,tau_ks,alpha_gks)
    Lambda_ks=Output_Mstep$Lambda_ks
    Psi_gs=Output_Mstep$Psi_gs
    Phi_gks=Output_Mstep$Phi_gks
    tau_ks=Output_Mstep$tau_ks
    alpha_gks=Output_Mstep$alpha_gks
    Sigma_gks=Output_Mstep$Sigma_gks
    invSigma_gks=Output_Mstep$invSigma_gks
    nractivatedconstraints=Output_Mstep$nractivatedconstraints

    S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        tlambda_k=t(lambda_k)
        mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tlambda_k # model-implied mean if group g is in cluster k
        S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
        S_gks[[g,k]] <- S_gk
        S_g=S_g+N_gks[g,k]*S_gk
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }

    ODLL=0;
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
    loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
    for(g in 1:ngroup){
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gks[[g,k]]*invSigma_gks[[g,k]])) # sum(S_gks[[g,k]]*invSigma_gks[[g,k]])=sum(diag(S_gks[[g,k]]%*%invSigma_gks[[g,k]]))
        loglik_gks[g,k]=loglik_gk
        loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
      }
      m_i=max(loglik_gksw[g,]);
      for(k in 1:nclust){
        loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
      }
      ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
    }

    pars=pi_ks
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      pars=c(pars,lambda_k[design==1])
    }
    pars=c(pars,lapply(Psi_gs,diag),Phi_gks,tau_ks,alpha_gks)
    pars=unlist(pars)
    ind=which(pars!=0 & prev_pars!=0)
    parsdiff <- pars[ind]-prev_pars[ind] #mapply("-",pars,prev_pars)
    parsdiffdiv <- parsdiff/prev_pars[ind] #mapply("/",parsdiff,prev_pars)
    parsabsdiffdiv <- abs(parsdiffdiv)#lapply(parsdiffdiv,function(x){abs(x)})
    conv1 <- sum(parsabsdiffdiv)#do.call(sum,parsabsdiffdiv)
    if(is.na(conv1)){
      conv1=1
    }

    conv2=ODLL-prev_ODLL
    if(ODLL-prev_ODLL<0){
     print(ODLL-prev_ODLL)
    }
    bestloglik=ODLL


  } # end while-loop till convergence

  if (conv2<1e-6){
    convergence=2 # convergence in terms of loglikelihood
  }
  else {
    if (conv1<1e-6) {
      convergence=1 # convergence in terms of parameters
    }
    else {
      convergence=0 # no convergence
    }
  }


  # set scale of factors across groups PER CLUSTER (and, in case of EFA, make them orthogonal)
  for(k in 1:nclust){
    if(N_ks[k]>0){
      theta_k=matrix(0,nfactors,nfactors)
      for(g in 1:ngroup){
        theta_gk=Theta_gks[[g,k]]
        theta_k=theta_k+(N_gks[g,k]/N_ks[k])*theta_gk;
      }
      if(nfactors>1){
        if(EFA==1){
          # find matrix square root via eigenvalue decomposition
          ed=eigen(theta_k)
          sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
          invsqrtFscale=solve(sqrtFscale)
        }
        else {
          sqrtFscale=diag(diag(theta_k^(1/2)))
          invsqrtFscale=diag(diag((1/theta_k)^(1/2)))
        }
      }
      else {
        sqrtFscale=sqrt(theta_k)
        invsqrtFscale=1/sqrtFscale
      }

      for(g in 1:ngroup){
        phi_gk=Phi_gks[[g,k]]
        phi_gk=invsqrtFscale%*%phi_gk%*%invsqrtFscale;
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry
      }
      Lambda_ks[[k]]=Lambda_ks[[k]]%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
      for(g in 1:ngroup){ # transform the factor means accordingly
        alpha_gks[[g,k]]=alpha_gks[[g,k]]%*%invsqrtFscale
      }
    }
  }

  # identification of factor means and intercepts
  mean_alpha_ks=matrix(0,nclust,nfactors)
  for(g in 1:ngroup){
    for(k in 1:nclust){
      if(N_gks[g,k]>0){
        mean_alpha_ks[k,]=mean_alpha_ks[k,]+N_gks[g,k]/N_ks[k]*alpha_gks[[g,k]]
      }
    }
  }
  # Translation of factor means per cluster and cycle 1 update indicator intercepts (see De Roover, 2021)
  for(k in 1:nclust){
    for(g in 1:ngroup){
      alpha_gks[[g,k]]=alpha_gks[[g,k]]-mean_alpha_ks[k,]
    }
    if(N_ks[k]>0){
      lambda_k=Lambda_ks[[k]]
      tlambda_k=t(lambda_k)
      suminvSigma=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvSigma=matrix(0,1,nvar)
      for(g in 1:ngroup){
        if(N_gks[g,k]>0){
          invSigma_gk=invSigma_gks[[g,k]]
          summeansminusalphaLambdainvSigma=summeansminusalphaLambdainvSigma+N_gks[g,k]*(mean_gs[g,]-alpha_gks[[g,k]]%*%tlambda_k)%*%invSigma_gk
          suminvSigma=suminvSigma+N_gks[g,k]*invSigma_gk
        }
      }
      tau_ks[k,]=t(solve(suminvSigma,t(summeansminusalphaLambdainvSigma)))
    }
  }



  if(EFA==1){
    nrpars=nclust-1+(nvar*nfactors-(nfactors*(nfactors-1)*(1/2)))*nclust+(nfactors*(nfactors+1)/2)*(ngroup-nclust)+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }
  else {
    nrpars=nclust-1+sum(design)*nclust+(nfactors*(nfactors+1)/2)*ngroup-nclust*nfactors+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }


  output_list <- list(z_gks=z_gks,pi_ks=pi_ks,Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,tau_ks=tau_ks,alpha_gks=alpha_gks,bestloglik=bestloglik,logliks=logliks,nrpars=nrpars,convergence=convergence,nractivatedconstraints=nractivatedconstraints)

  return(output_list)
} # end main function



# functions for E-step/posterior classification probabilities and M-step


# Update the cluster-membership probabilities z_gk
# Reuses the loglik_gks to save time

UpdPostProb <- function(pi_ks, loglik_gks, ngroup, nclust, v=1){
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow=ngroup,ncol=nclust)

  for(g in 1:ngroup){
    for(k in 1:nclust){
      z_gks[g,k] <- v*log(pi_ks[k])+v*loglik_gks[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclust))
  }

  # divide by the rowwise sum of the above calculated part
  z_gks <- sweep(z_gks,1/rowSums(z_gks),MARGIN=1,'*',check.margin = FALSE)#diag(1/rowSums(z_gks))%*%z_gks
  z_gks <- round(z_gks,digits=16)
  z_gks <- sweep(z_gks,1/rowSums(z_gks),MARGIN=1,'*',check.margin = FALSE)#diag(1/rowSums(z_gks))%*%z_gks

  return(z_gks)
}



mixmgfa_loadingsintercepts_Mstep <- function(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,meanexpEta_gks,Lambda_ks,Psi_gs,Phi_gks,mean_gs,tau_ks,alpha_gks){
  nractivatedconstraints <- 0
  ngroup <- length(N_gs)
  N_ks=colSums(N_gks)

  if(sum(design)==nvar*nfactors){ # if design contains only '1's, EFA is used
    EFA=1
  } else {
    EFA=0
  }

  invPsi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
  for(g in 1:ngroup){
    psi_g=Psi_gs[[g]]
    invPsi_g=diag(1/diag(psi_g))
    invPsi_gs[[g]] <- invPsi_g
  }

  # update indicator intercepts
  for(k in 1:nclust){
    if(N_ks[k]>0){
      lambda_k=Lambda_ks[[k]]
      tlambda_k=t(lambda_k)
      suminvPsi=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvPsi=matrix(0,1,nvar)
      for(g in 1:ngroup){
        if(N_gks[g,k]>0){
          invPsi_g=invPsi_gs[[g]]
          summeansminusalphaLambdainvPsi=summeansminusalphaLambdainvPsi+N_gks[g,k]*(mean_gs[g,]-(alpha_gks[[g,k]]+meanexpEta_gks[[g,k]])%*%tlambda_k)%*%invPsi_g
          suminvPsi=suminvPsi+N_gks[g,k]*invPsi_g
        }
      }
      tau_ks[k,]=t(solve(suminvPsi,t(summeansminusalphaLambdainvPsi)))
    }
  }

  # update factor means
  for(g in 1:ngroup){
    invPsi_g=invPsi_gs[[g]]
    for(k in 1:nclust){
      if(N_ks[k]>0){
        lambda_k=Lambda_ks[[k]]
        tlambda_k=t(lambda_k)
        tlambda_kinvPsi_g=tlambda_k%*%invPsi_g
        alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-meanexpEta_gks[[g,k]]%*%tlambda_k)%*%t(solve((tlambda_kinvPsi_g%*%lambda_k),(tlambda_kinvPsi_g)))
      }
    }
  }


  # update factor loadings
  if(EFA==1){
    for(k in 1:nclust){
      if(N_ks[k]>0){
        lambda_k=matrix(0,nvar,nfactors)
        sumSbetas=matrix(list(0),1,nvar)
        sumthetaalphas=matrix(list(0),1,nvar)
        summeansalphas=matrix(list(0),1,nvar)
        for(g in 1:ngroup){
          if(N_gks[g,k]>0){
            psi_g=Psi_gs[[g]]
            S_gk=S_gks[[g,k]]
            beta_gk=Beta_gks[[g,k]]
            theta_gk=Theta_gks[[g,k]]
            alpha_gk=alpha_gks[[g,k]]
            meanexpeta_gk=meanexpEta_gks[[g,k]]
            #talpha_gk=t(alpha_gk)
            N_gks_psi_g=N_gks[g,k]/diag(psi_g)
            S_gktbeta_gk=tcrossprod(S_gk,beta_gk)
            theta_gk_latmeanscp=theta_gk+crossprod(alpha_gk+meanexpeta_gk,alpha_gk)
            meansalpha=((mean_gs[g,]-tau_ks[k,])%*%alpha_gk)
            for(j in 1:nvar){
              sumSbetas[[j]]=sumSbetas[[j]]+N_gks_psi_g[j]*S_gktbeta_gk[j, ,drop=FALSE]
              sumthetaalphas[[j]]=sumthetaalphas[[j]]+N_gks_psi_g[j]*theta_gk_latmeanscp
              summeansalphas[[j]]=summeansalphas[[j]]+N_gks_psi_g[j]*meansalpha[j, ,drop=FALSE]
            }
          }
        }
        for(j in 1:nvar){
          lambda_k[j,]= t(solve(sumthetaalphas[[j]],t(sumSbetas[[j]]+summeansalphas[[j]])))
        }
        Lambda_ks[[k]]=lambda_k
      }
    }
  } else {
    # for(k in 1:nclust){
    #   if(N_ks[k]>0){
    #     lambda_k=matrix(0,nvar,nfactors)
    #     sumSbetas=matrix(list(0),1,nvar)
    #     sumthetaalphas=matrix(list(0),1,nvar)
    #     summeansalphas=matrix(list(0),1,nvar)
    #     for(g in 1:ngroup){
    #       if(N_gks[g,k]>0){
    #         psi_g=Psi_gs[[g]]
    #         S_gk=S_gks[[g,k]]
    #         beta_gk=Beta_gks[[g,k]]
    #         theta_gk=Theta_gks[[g,k]]
    #         alpha_gk=alpha_gks[[g,k]]
    #         meanexpeta_gk=meanexpEta_gks[[g,k]]
    #         #talpha_gk=t(alpha_gk)
    #         N_gks_psi_g=N_gks[g,k]/diag(psi_g)
    #         S_gktbeta_gk=S_gk%*%t(beta_gk)
    #         theta_gk_latmeanscp=(theta_gk+t(alpha_gk+meanexpeta_gk)%*%alpha_gk)
    #         meansalpha=((mean_gs[g,]-tau_ks[k,])%*%alpha_gk)
    #         for(j in 1:nvar){ # simultaneously for all variables with same non-zero loadings? Based on unique(design,MARGIN=1)?
    #           # nfactors_j=sum(design[j,])
    #           # indfactors_j=(design[j,]==1)
    #           # if(nfactors_j<nfactors){
    #           #   beta_gk=beta_gk[indfactors_j, ,drop=FALSE]
    #           # }
    #           # theta_gk=theta_gk[indfactors_j,indfactors_j]
    #           # alpha_gk=alpha_gk[indfactors_j]
    #           # meanexpeta_gk=meanexpeta_gk[indfactors_j]
    #           # if(g==1){
    #           #   sumSbetas[[j]]=matrix(0,1,nfactors)
    #           #   sumthetaalphas[[j]]=matrix(0,nfactors,nfactors)
    #           #   summeansalphas[[j]]=matrix(0,1,nfactors)
    #           # }
    #           sumSbetas[[j]]=sumSbetas[[j]]+N_gks_psi_g[j]*S_gktbeta_gk[j, ,drop=FALSE]
    #           sumthetaalphas[[j]]=sumthetaalphas[[j]]+N_gks_psi_g[j]*theta_gk_latmeanscp
    #           summeansalphas[[j]]=summeansalphas[[j]]+N_gks_psi_g[j]*meansalpha[j, ,drop=FALSE]
    #         }
    #
    #       }
    #     }
    #     for(j in 1:nvar){
    #       #indfactors_j=(design[j,]==1)
    #       lambda_k[j,]= t(solve(sumthetaalphas[[j]],t(sumSbetas[[j]]+summeansalphas[[j]])))
    #     }
    #     Lambda_ks[[k]]=lambda_k
    #   }
    # }
    for(k in 1:nclust){
      if(N_ks[k]>0){
        lambda_k=matrix(0,nvar,nfactors)
        for(j in 1:nvar){
          nfactors_j=sum(design[j,])
          indfactors_j=(design[j,]==1)
          sumSbeta=matrix(0,1,nfactors_j)
          sumthetaalpha=matrix(0,nfactors_j,nfactors_j)
          summeansalpha=matrix(0,1,nfactors_j)
          for(g in 1:ngroup){
            if(N_gks[g,k]>0){
              psi_g=diag(Psi_gs[[g]])
              S_gk=S_gks[[g,k]]
              beta_gk=Beta_gks[[g,k]]
              theta_gk=Theta_gks[[g,k]]
              alpha_gk=alpha_gks[[g,k]]
              meanexpeta_gk=meanexpEta_gks[[g,k]]
              if(nfactors_j<nfactors){
                beta_gk=beta_gk[indfactors_j, ,drop=FALSE]
              }
              theta_gk=theta_gk[indfactors_j,indfactors_j]
              alpha_gk=alpha_gk[indfactors_j]
              meanexpeta_gk=meanexpeta_gk[indfactors_j]
              sumSbeta=sumSbeta+(N_gks[g,k]/psi_g[j])*tcrossprod(S_gk[j,],beta_gk)
              sumthetaalpha=sumthetaalpha+(N_gks[g,k]/psi_g[j])*(theta_gk+tcrossprod(alpha_gk+meanexpeta_gk,alpha_gk))
              summeansalpha=summeansalpha+(N_gks[g,k]/psi_g[j])*((mean_gs[g,j]-tau_ks[k,j])%*%alpha_gk)
            }
          }
          if(nfactors_j==1){
            lambda_k[j,indfactors_j]= (sumSbeta+summeansalpha)/sumthetaalpha
          } else {
            lambda_k[j,indfactors_j]= t(solve(sumthetaalpha,t(sumSbeta+summeansalpha)))
          }
        }
        Lambda_ks[[k]]=lambda_k
      }
    }
  }



  # update unique variances
  nractivatedconstraints=0
  for(g in 1:ngroup){
    S_g=S_gs[[g]]
    sum2SbetaB_BthetaB=0
    for(k in 1:nclust){
      if(N_gks[g,k]>0){
        lambda_k=Lambda_ks[[k]]
        beta_gk=Beta_gks[[g,k]]
        theta_gk=Theta_gks[[g,k]]
	      S_gk=S_gks[[g,k]]
        #sum2SbetaB_BthetaB=sum2SbetaB_BthetaB+(N_gks[g,k]/N_gs[g])*(lambda_k%*%(2*beta_gk%*%S_gk-theta_gk%*%t(lambda_k))) # modelimplied reduced covariance matrix on sample level, based on old structure matrix and sigma_gk, weighting based on new z_gks
        sum2SbetaB_BthetaB=sum2SbetaB_BthetaB+(N_gks[g,k]/N_gs[g])*rowSums(lambda_k*(2*tcrossprod(S_gk,beta_gk)-lambda_k%*%theta_gk)) # faster way (result = vector), only works when no residual covariances!
      }
    }
    psi_g=diag(S_g)-sum2SbetaB_BthetaB # diag(diag(S_g-sum2SbetaB_BthetaB))
    d=psi_g #diag(psi_g)
    if (sum(d<.0001)>0){ # track "heywood" cases
      ind=d<.0001
      d[ind]=0.0001
      psi_g=d # diag(psi_g)=d
      nractivatedconstraints=nractivatedconstraints+sum(ind)
    }
    Psi_gs[[g]]=diag(psi_g) #psi_g
  }

  # update factor (co)variances
  for(g in 1:ngroup){
    for(k in 1:nclust){
      if(N_ks[k]>0){
        theta_gk=Theta_gks[[g,k]]
        phi_gk=theta_gk
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry to avoid accumulation of asymmetry over iterations
      }
    }
  }


  # update (inv)Sigma_gks
  Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
  invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
  for(k in 1:nclust){
    lambda_k=Lambda_ks[[k]]
    tlambda_k=t(lambda_k)
    for(g in 1:ngroup){
      psi_g=Psi_gs[[g]]
      invPsi_g=diag(1/diag(psi_g))
      invPsi_g_lambda_k=invPsi_g%*%lambda_k
      phi_gk=Phi_gks[[g,k]]
      invPhi_gk=solve(phi_gk)
      sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
      Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
      invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk+tlambda_k%*%invPsi_g_lambda_k)
      invPhi_gk_tLambdainvPsi_g_Lambda=(invPhi_gk_tLambdainvPsi_g_Lambda+t(invPhi_gk_tLambdainvPsi_g_Lambda))*(1/2)
      invSigma_gks[[g,k]]=invPsi_g-invPsi_g_lambda_k%*%solve(invPhi_gk_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda_k) # Woodbury identity
    }
  }

  output_list <- list(Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,tau_ks=tau_ks,alpha_gks=alpha_gks,Sigma_gks=Sigma_gks,invSigma_gks=invSigma_gks,nractivatedconstraints=nractivatedconstraints)

  return(output_list)
}

# computation of adjusted rand index
adjrandindex <- function(part1,part2){

  IM1=diag(max(part1))
  IM2=diag(max(part2))
  A=IM1[part1,]
  B=IM2[part2,]

  T = t(A)%*%B
  N = sum(T)
  Tc = colSums(T)
  Tr = rowSums(T)
  a = (sum(T^2) - N)/2
  b = (sum(Tr^2) - sum(T^2))/2
  c = (sum(Tc^2) - sum(T^2))/2
  d = (sum(T^2) + N^2 - sum(Tr^2) - sum(Tc^2))/2
  ARI = (choose(N,2)*(a + d) - ((a+b)*(a+c)+(c+d)*(b+d)))/(choose(N,2)^2 - ((a+b)*(a+c)+(c+d)*(b+d)))

  return(ARI)
}

# Procrustes rotation (orthogonal)
procr <- function(x,y){
  s <- svd(t(x)%*%y)
  U <- s$u # X = U D V'
  D <- s$d
  V <- s$v
  R <- U%*%t(V)
  yhat <- x%*%R # rotated x that approximates y

  return(yhat)
}

stirlingnr2 <- function(n,k){

  # The number of ways of partitioning a set of n elements into k nonempty sets (Stirling number of second kind).
  # For example, the set {1,2,3} can be partitioned into three subsets in one way: {{1},{2},{3}}; into two subsets in three ways: {{1,2},{3}}, {{1,3},{2}},
  # and {{1},{2,3}}; and into one subset in one way: {{1,2,3}}.


  S=0
  for(i in 0:k){
    S=S+(-1)^i*round(factorial(k),digits=15)/(round(factorial(i),digits=15)*round(factorial(k-i),digits=15))*(k-i)^n
  }
  S=S/round(factorial(k),digits=15)


  if(S<0){
    S=1
  }

  return(S)
}


