#' Mixture multigroup factor analysis for intercept (non-)variance
# ----------------------------------------------------------------
# Code written by Kim De Roover

#' @description
#' Finds clusters of groups based on their intercepts, building on metric invariance, given a user-specified number of clusters.
#' Loadings are invariant and residual variances are group-specific. Factor (co)variances are group-specific. Factor means are group- and cluster-specific. See referenced paper for details.

# INPUT:
#' @param data A list of "$covariances" (a vertically concatenated matrix or list of group-specific (co)variance matrices) and "$means" (a matrix with rows = group-specific means); or a matrix containing the vertically concatenated raw data for all groups.
#' @param N_gs Vector containing the sample size for each group (in the same order as they appear in the data).
#' @param nclust User-specified number of clusters.
#' @param nfactors User-specified number of factors.
#' @param maxiter Maximum number of iterations. Increase in case of non-convergence.
#' @param start Type of start (start = 1: pre-selected random starts, start = 2: start from a user-specified startpartition).
#' @param nruns Number of starts (based on pre-selected random partitions when start = 1).
#' @param preselect Percentage of best starts taken in pre-selection (increase to speed up startprocedure).
#' @param design For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).
#' @param startpartition Partition of groups (vector) to start from (use with start = 2 and nruns = 1).

# OUTPUT:
#' @return Output object (list) with:
#'
#' $z_gks = cluster memberships of groups (posterior classification probabilities)
#'
#' $pi_ks = mixing proportions (prior classification probabilities)
#'
#' $Lambda = invariant loadings
#'
#' $Psi_gs = group-specific unique variances, access unique variances of group g via Psi_gs[[g]]
#'
#' $Phi_gs = group-specific factor (co)variances, access (co)variances of group g via Phi_gs[[g]]
#'
#' $tau_ks = cluster-specific intercepts, access intercepts of cluster k via tau_ks[k,]
#'
#' $alpha_gks = group- and cluster-specific factor means, access factor means of group g in cluster k via alpha_gks[[g,k]]
#'
#' $bestloglik = final loglikelihood, loglikelihood of best start
#'
#' $logliks = loglikelihoods (first column) and number of activated constraints on unique variances (second column) for all starts
#'
#' $nrpars = number of free parameters, to be used for model selection in combination with bestloglik
#'
#' $convergence = 2 if best start converged on loglikelihood, 1 if converged on parameter changes, 0 if not converged
#'
#' $nractivatedconstraints = number of constraints on the unique variances (across groups, for best start) to avoid unique variances approaching zero

#' @references
#' De Roover, K. (2021). Finding clusters of groups with measurement invariance: Unraveling intercept non-invariance with mixture multigroup factor analysis. Structural Equation Modeling: A Multidisciplinary Journal, 28(5), 663-683.
#'



#' @export
mixmgfa_intercepts <- function(data,N_gs,nclust,nfactors=1,maxiter = 5000,start = 1,nruns = 25,design = 0,preselect = 10,startpartition){

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
      obsS_gs[[g]] <- CP_g-tcrossprod(mean_g,mean_g) #mean_g%*%t(mean_g) #(1/N_gs[g])*(t(X_gc)%*%X_gc)
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

          tau_ks <- matrix(0,nclust,nvar)
          for(k in 1:nclust){
            for(g in 1:ngroup){
              if(randpartvec[g]==k){
                tau_ks[k,]=tau_ks[k,]+(N_gks[g,k]/N_ks[k])*mean_gs[g,]
              }
            }
          }
          # # center data with initialized intercepts
          # Xsupcent <- matrix(0,N,nvar)
          # for(g in 1:ngroup){
          #   k=randpartvec[g]
          #   X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
          #   Xsupcent[Ncum[g,1]:Ncum[g,2],]=X_g-(matrix(tau_ks[k,],ncol=nvar,nrow=N_gs[g],byrow=TRUE))
          # }
          # S <- (1/N)*(t(Xsupcent)%*%Xsupcent)
          S=matrix(0,nvar,nvar)
          for(g in 1:ngroup){
            k=randpartvec[g]
            mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
            mu_gk=matrix(tau_ks[k,],nrow=1,ncol=nvar)
            S=S+N_gs[g]*(CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g))
          }
          S=(1/N)*S
          ed<-eigen(S, symmetric=TRUE, only.values = FALSE)
          val<-ed$values
          u<-ed$vectors
          totalerror=sum((val[-seq_len(nfactors)]))
          meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
          Uniq=rep(meanerror,nvar)
          Lambda=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
          if (EFA==0){
            Lambda <- procr(Lambda,design)
            Lambda=Lambda*design # non-zero loadings should be indicated with '1' for this to work properly
          }


          Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
          Phi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific factor covariances
          alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
          for(g in 1:ngroup){
            Psi_gs[[g]]=diag(Uniq)
            Phi_gs[[g]]=diag(nfactors)
            if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
              X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
              for(k in 1:nclust){
                X_c=X_g-t(matrix(tau_ks[k,],ncol=N_gs[g],nrow=nvar))
                udv=svd(X_c%*%Lambda) # this is a truncated svd
                U=udv$u
                V=udv$v
                Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
                Fvar=apply(Fscores,2,var)
                Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
                alpha_gks[[g,k]] <- colMeans(Fscores)
              }
            } else {
              tLambda=t(Lambda)
              S_g=obsS_gs[[g]]
              tLambdaS_g=tLambda%*%solve(S_g)
              for(k in 1:nclust){
                alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambdaS_g)%*%Lambda),(tLambdaS_g)))
              }
            }
          }

          Sigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
          invSigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
          tLambda=t(Lambda)
          for(g in 1:ngroup){
            psi_g=Psi_gs[[g]]
            invPsi_g=diag(1/diag(psi_g))
            invPsi_g_lambda=invPsi_g%*%Lambda
            phi_g=Phi_gs[[g]]
            invPhi_g=phi_g # phi_gk is still identity matrix
            sigma_g=Lambda %*% phi_g %*% tLambda + psi_g
            Sigma_gs[[g]]=(sigma_g+t(sigma_g))*(1/2) # avoid asymmetry due to rounding errors
            invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g+tLambda%*%invPsi_g_lambda)
            invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
            #invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda) # Woodbury identity
            invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda,tLambda)%*%invPsi_g # Woodbury identity
          }

          # if (is.matrix(data)==FALSE & is.data.frame(data)==FALSE){
          #   # cycle 1 update factor means
          #   for(g in 1:ngroup){
          #     invSigma_g=invSigma_gs[[g]]
          #     for(k in 1:nclust){
          #       if(N_ks[k]>0){
          #         alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambda%*%invSigma_g)%*%Lambda),(tLambda%*%invSigma_g)))
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
              mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
              S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
              S_gks[[g,k]] <- S_gk
              S_g=S_g+N_gks[g,k]*S_gk
            }
            S_gs[[g]] <- (1/N_gs[g])*S_g
          }

          # compute Beta_gs and theta_gks
          Beta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
          Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
          Theta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
          meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
          #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
          for(g in 1:ngroup){
            invsigma_g=invSigma_gs[[g]]
            phi_g=Phi_gs[[g]]
            lambda_phi_g=Lambda%*%phi_g
            beta_g=phi_g%*%crossprod(Lambda,invsigma_g) #t(Lambda)%*%invsigma_g
            tbeta_g=t(beta_g)
            Beta_gs[[g]]=beta_g;
            for(k in 1:nclust){
              S_gk=S_gks[[g,k]]
              theta_gk=phi_g-beta_g%*%(lambda_phi_g-S_gk%*%tbeta_g)
              Theta_gks[[g,k]]=theta_gk
              meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-tcrossprod(alpha_gks[[g,k]],Lambda))%*%tbeta_g
              #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
            }
            S_g=S_gs[[g]]
            Theta_gs[[g]]=phi_g-beta_g%*%(lambda_phi_g-S_g%*%tbeta_g)
          }

          Output_Mstep <- mixmgfa_intercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gs,Theta_gks,Theta_gs,meanexpEta_gks,Lambda,Psi_gs,Phi_gs,mean_gs,tau_ks,alpha_gks)
          Lambda=Output_Mstep$Lambda
          Psi_gs=Output_Mstep$Psi_gs
          Phi_gs=Output_Mstep$Phi_gs
          tau_ks=Output_Mstep$tau_ks
          alpha_gks=Output_Mstep$alpha_gks
          Sigma_gs=Output_Mstep$Sigma_gs
          invSigma_gs=Output_Mstep$invSigma_gs
          nractivatedconstraints=Output_Mstep$nractivatedconstraints

          S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
          S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
          tLambda=t(Lambda)
          for(g in 1:ngroup){
            S_g=matrix(0,nvar,nvar)
            mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
            for(k in 1:nclust){
              mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
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
            logdet_sigma_g=log(det(Sigma_gs[[g]]))
            invSigma_g=invSigma_gs[[g]]
            # X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
            # #tX_g=t(X_g)
            # for(k in 1:nclust){
            #   Xc_gk=X_g-t(matrix(tau_ks[k,]+alpha_gks[[g,k]]%*%t(Lambda),ncol=N_gs[g],nrow=nvar)) # centered data per group
            #   #tXc_gk=t(Xc_gk)
            #   loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g))
            #   for (n in 1:N_gs[g]){
            #     #loglik_gk=loglik_gk-(1/2)*(Xc_gk[n,]%*%invSigma_g%*%tXc_gk[,n])
            #     Xc_n=Xc_gk[n, ,drop=FALSE]
            #     loglik_gk=loglik_gk-(1/2)*(Xc_n%*%tcrossprod(invSigma_g,Xc_n))
            #   }
            #   # Xc_gk_list=split.data.frame(Xc_gk,seq_len(nrow(Xc_gk)),'list')
            #   # Xc_n_invSigma_g_tXc_n=lapply(Xc_gk_list,function(x){-(1/2)*(x%*%invSigma_g%*%t(x))})
            #   # loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g))+do.call(sum,Xc_n_invSigma_g_tXc_n)
            #   loglik_gks[g,k]=loglik_gk
            #   loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
            # }
            for(k in 1:nclust){
              loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
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
    S=matrix(0,nvar,nvar)
    for(g in 1:ngroup){
      k=randpartvec[g]
      mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
      mu_gk=matrix(tau_ks[k,],nrow=1,ncol=nvar)
      S=S+N_gs[g]*(CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g))
    }
    S=(1/N)*S
    ed<-eigen(S, symmetric=TRUE, only.values = FALSE)
    val<-ed$values
    u<-ed$vectors
    totalerror=sum((val[-seq_len(nfactors)]))
    meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
    Uniq=rep(meanerror,nvar)
    Lambda=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
    if (EFA==0){
      Lambda <- procr(Lambda,design)
      Lambda=Lambda*design # non-zero loadings should be indicated with '1' for this to work properly
    }


    Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
    Phi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific factor covariances
    alpha_gks <- matrix(list(NA), nrow=ngroup, ncol = nclust) # initialize group- and cluster-specific factor means
    for(g in 1:ngroup){
      Psi_gs[[g]]=diag(Uniq)
      Phi_gs[[g]]=diag(nfactors)
      if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){
        X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
        for(k in 1:nclust){
          X_c=X_g-t(matrix(tau_ks[k,],ncol=N_gs[g],nrow=nvar))
          udv=svd(X_c%*%Lambda) # this is a truncated svd
          U=udv$u
          V=udv$v
          Fscores=U[,seq_len(nfactors),drop=FALSE]%*%t(V) # compute component scores as initial estimates of factor scores
          Fvar=apply(Fscores,2,var)
          Fscores=scale(Fscores,center=FALSE,scale=sqrt(Fvar))
          alpha_gks[[g,k]] <- colMeans(Fscores)
        }
      } else {
        tLambda=t(Lambda)
        S_g=obsS_gs[[g]]
        tLambdaS_g=tLambda%*%solve(S_g)
        for(k in 1:nclust){
          alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambdaS_g)%*%Lambda),(tLambdaS_g)))
        }
      }
    }

    Sigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
    invSigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
    tLambda=t(Lambda)
    for(g in 1:ngroup){
      psi_g=Psi_gs[[g]]
      invPsi_g=diag(1/diag(psi_g))
      invPsi_g_lambda=invPsi_g%*%Lambda
      phi_g=Phi_gs[[g]]
      invPhi_g=solve(phi_g)
      sigma_g=Lambda %*% phi_g %*% tLambda + psi_g
      Sigma_gs[[g]]=(sigma_g+t(sigma_g))*(1/2) # avoid asymmetry due to rounding errors
      invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g+tLambda%*%invPsi_g_lambda)
      invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
      #invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda) # Woodbury identity
      invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda,tLambda)%*%invPsi_g # Woodbury identity
    }

    # if (is.matrix(data)==FALSE & is.data.frame(data)==FALSE){
    #   # cycle 1 update factor means
    #   for(g in 1:ngroup){
    #     invSigma_g=invSigma_gs[[g]]
    #     for(k in 1:nclust){
    #       if(N_ks[k]>0){
    #         alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambda%*%invSigma_g)%*%Lambda),(tLambda%*%invSigma_g)))
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
        mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
        S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
        S_gks[[g,k]] <- S_gk
        S_g=S_g+N_gks[g,k]*S_gk
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }

    # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
    for(g in 1:ngroup){
      logdet_sigma_g=log(det(Sigma_gs[[g]]))
      invSigma_g=invSigma_gs[[g]]
      for(k in 1:nclust){
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
        loglik_gks[g,k]=loglik_gk
      }
    }

    iter=0;
    conv1=1;
    conv2=1;
    ODLL=-Inf;
    pars=c(pi_ks,Lambda[design==1],unlist(lapply(Psi_gs,diag)),unlist(Phi_gs),as.vector(tau_ks),unlist(alpha_gks))
    while(min(conv1,conv2)>1e-4 && iter<100){
      prev_ODLL=ODLL
      prev_Lambda=Lambda
      prev_Psi_gs=Psi_gs
      prev_Phi_gs=Phi_gs
      prev_tau_ks=tau_ks
      prev_alpha_gks=alpha_gks
      prev_pars=pars
      iter=iter+1

      # ** CYCLE 1**
      # **E-step**: compute the posterior classification probabilities
      if(nclust>1 && nclust<ngroup){
        z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
        pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
      }

      N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
      N_ks=colSums(N_gks)


      # cycle 1 update indicator intercepts
      tLambda=t(Lambda)
      for(k in 1:nclust){
        if(N_ks[k]>0){
          suminvSigma=matrix(0,nvar,nvar)
          summeansminusalphaLambdainvSigma=matrix(0,1,nvar)
          for(g in 1:ngroup){
            if(N_gks[g,k]>0){
              invSigma_g=invSigma_gs[[g]]
              summeansminusalphaLambdainvSigma=summeansminusalphaLambdainvSigma+N_gks[g,k]*(mean_gs[g,]-alpha_gks[[g,k]]%*%tLambda)%*%invSigma_g
              suminvSigma=suminvSigma+N_gks[g,k]*invSigma_g
            }
          }
          tau_ks[k,]=t(solve(suminvSigma,t(summeansminusalphaLambdainvSigma)))
        }
      }

      # cycle 1 update factor means
      for(g in 1:ngroup){
        invSigma_g=invSigma_gs[[g]]
        tLambdainvSigma_g=tLambda%*%invSigma_g
        for(k in 1:nclust){
          if(N_ks[k]>0){
            alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambdainvSigma_g)%*%Lambda),(tLambdainvSigma_g)))
          }
        }
      }

      S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
        for(k in 1:nclust){
          mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
          S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
          S_gks[[g,k]] <- S_gk
          S_g=S_g+N_gks[g,k]*S_gk
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }

      # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
      loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
      for(g in 1:ngroup){
        logdet_sigma_g=log(det(Sigma_gs[[g]]))
        invSigma_g=invSigma_gs[[g]]
        for(k in 1:nclust){
          loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
          loglik_gks[g,k]=loglik_gk
        }
      }

      # ** CYCLE 2**
      # **E-step**: compute the posterior classification probabilities
      if(nclust>1 && nclust<ngroup){
        z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
        pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
      }

      N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
      N_ks=colSums(N_gks)

      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
	      for(k in 1:nclust){
            S_g=S_g+N_gks[g,k]*S_gks[[g,k]]
        }
        S_gs[[g]] <- (1/N_gs[g])*S_g
      }

      # compute Beta_gs and theta_gks
      Beta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
      Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      Theta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
      meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
      #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
      for(g in 1:ngroup){
        invsigma_g=invSigma_gs[[g]]
        phi_g=Phi_gs[[g]]
        lambda_phi_g=Lambda%*%phi_g
        beta_g=phi_g%*%crossprod(Lambda,invsigma_g) #t(Lambda)%*%invsigma_g
        tbeta_g=t(beta_g)
        Beta_gs[[g]]=beta_g;
        for(k in 1:nclust){
          S_gk=S_gks[[g,k]]
          theta_gk=phi_g-beta_g%*%(lambda_phi_g-S_gk%*%tbeta_g)
          Theta_gks[[g,k]]=theta_gk
          meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-tcrossprod(alpha_gks[[g,k]],Lambda))%*%tbeta_g
          #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
        }
        S_g=S_gs[[g]]
        Theta_gs[[g]]=phi_g-beta_g%*%(lambda_phi_g-S_g%*%tbeta_g)
      }

      Output_Mstep <- mixmgfa_intercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gs,Theta_gks,Theta_gs,meanexpEta_gks,Lambda,Psi_gs,Phi_gs,mean_gs,tau_ks,alpha_gks)
      Lambda=Output_Mstep$Lambda
      Psi_gs=Output_Mstep$Psi_gs
      Phi_gs=Output_Mstep$Phi_gs
      tau_ks=Output_Mstep$tau_ks
      alpha_gks=Output_Mstep$alpha_gks
      Sigma_gs=Output_Mstep$Sigma_gs
      invSigma_gs=Output_Mstep$invSigma_gs
      nractivatedconstraints=Output_Mstep$nractivatedconstraints

      S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      tLambda=t(Lambda)
      for(g in 1:ngroup){
        S_g=matrix(0,nvar,nvar)
        mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
        for(k in 1:nclust){
          mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
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
        logdet_sigma_g=log(det(Sigma_gs[[g]]))
        invSigma_g=invSigma_gs[[g]]
        for(k in 1:nclust){
          loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
          loglik_gks[g,k]=loglik_gk
          loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
        }
        m_i=max(loglik_gksw[g,]);
        for(k in 1:nclust){
          loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
        }
        ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
      }


      pars=c(pi_ks,Lambda[design==1],unlist(lapply(Psi_gs,diag)),unlist(Phi_gs),as.vector(tau_ks),unlist(alpha_gks))
      parsdiff <- mapply("-",pars,prev_pars)
      parsdiffdiv <- mapply("/",parsdiff,prev_pars)
      parsabsdiffdiv <- abs(parsdiffdiv)
      conv1 <- sum(parsabsdiffdiv)
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
      bestLambda=Lambda
      bestPsi_gs=Psi_gs
      bestPhi_gs=Phi_gs
      besttau_ks=tau_ks
      bestalpha_gks=alpha_gks
      bestSigma_gs=Sigma_gs
      bestinvSigma_gs=invSigma_gs
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
        bestLambda=Lambda
        bestPsi_gs=Psi_gs
        bestPhi_gs=Phi_gs
        besttau_ks=tau_ks
        bestalpha_gks=alpha_gks
        bestSigma_gs=Sigma_gs
        bestinvSigma_gs=invSigma_gs
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
  Lambda=bestLambda
  Psi_gs=bestPsi_gs
  Phi_gs=bestPhi_gs
  tau_ks=besttau_ks
  alpha_gks=bestalpha_gks
  Sigma_gs=bestSigma_gs
  invSigma_gs=bestinvSigma_gs
  ODLL=bestloglik
  loglik_gks=bestloglik_gks
  iter=bestiter
  conv1=bestconv1
  conv2=bestconv2

  pars=c(pi_ks,Lambda[design==1],unlist(lapply(Psi_gs,diag)),unlist(Phi_gs),as.vector(tau_ks),unlist(alpha_gks))
  while(min(conv1,conv2)>1e-6 && iter<maxiter+1){ # iterate till convergence for best start
    prev_ODLL=ODLL
    prev_Lambda=Lambda
    prev_Psi_gs=Psi_gs
    prev_Phi_gs=Phi_gs
    prev_tau_ks=tau_ks
    prev_alpha_gks=alpha_gks
    prev_pars=pars
    iter=iter+1

    # **CYCLE 1**
    # **E-step**: compute the posterior classification probabilities
    if(nclust>1 && nclust<ngroup){
      z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
      pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
    }

    N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=colSums(N_gks)


    # cycle 1 update indicator intercepts
    tLambda=t(Lambda)
    for(k in 1:nclust){
      if(N_ks[k]>0){
        suminvSigma=matrix(0,nvar,nvar)
        summeansminusalphaLambdainvSigma=matrix(0,1,nvar)
        for(g in 1:ngroup){
          if(N_gks[g,k]>0){
            invSigma_g=invSigma_gs[[g]]
            summeansminusalphaLambdainvSigma=summeansminusalphaLambdainvSigma+N_gks[g,k]*(mean_gs[g,]-alpha_gks[[g,k]]%*%tLambda)%*%invSigma_g
            suminvSigma=suminvSigma+N_gks[g,k]*invSigma_g
          }
        }
        tau_ks[k,]=t(solve(suminvSigma,t(summeansminusalphaLambdainvSigma)))
      }
    }

    # cycle 1 update factor means
    for(g in 1:ngroup){
      invSigma_g=invSigma_gs[[g]]
      tLambdainvSigma_g=tLambda%*%invSigma_g
      for(k in 1:nclust){
        if(N_ks[k]>0){
          alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,])%*%t(solve(((tLambdainvSigma_g)%*%Lambda),(tLambdainvSigma_g)))
        }
      }
    }

    S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
      for(k in 1:nclust){
        mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
        S_gk=CP_gs[[g]]-crossprod(mean_g,mu_gk)+crossprod(mu_gk,mu_gk-mean_g)
        S_gks[[g,k]] <- S_gk
        S_g=S_g+N_gks[g,k]*S_gk
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }

    # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
    for(g in 1:ngroup){
      logdet_sigma_g=log(det(Sigma_gs[[g]]))
      invSigma_g=invSigma_gs[[g]]
      for(k in 1:nclust){
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
        loglik_gks[g,k]=loglik_gk
      }
    }


    # ** CYCLE 2**
    # **E-step**: compute the posterior classification probabilities
    if(nclust>1 && nclust<ngroup){
      z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust)
      pi_ks=(1/ngroup)*colSums(z_gks) # update mixing proportions
    }

    N_gks=sweep(z_gks,N_gs,MARGIN=1,'*',check.margin = FALSE) #N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=colSums(N_gks)


    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      for(k in 1:nclust){
        S_g=S_g+N_gks[g,k]*S_gks[[g,k]]
      }
      S_gs[[g]] <- (1/N_gs[g])*S_g
    }

    # compute Beta_gs and theta_gks
    Beta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
    Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    Theta_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
    meanexpEta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # mean of expected eta's for each group-cluster combination
    #meanexpEta_gs <- matrix(list(0), nrow = ngroup, ncol = 1) # mean of expected eta's for each group
    for(g in 1:ngroup){
      invsigma_g=invSigma_gs[[g]]
      phi_g=Phi_gs[[g]]
      lambda_phi_g=Lambda%*%phi_g
      beta_g=phi_g%*%crossprod(Lambda,invsigma_g) #t(Lambda)%*%invsigma_g
      tbeta_g=t(beta_g)
      Beta_gs[[g]]=beta_g;
      for(k in 1:nclust){
        S_gk=S_gks[[g,k]]
        theta_gk=phi_g-beta_g%*%(lambda_phi_g-S_gk%*%tbeta_g)
        Theta_gks[[g,k]]=theta_gk
        meanexpEta_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-tcrossprod(alpha_gks[[g,k]],Lambda))%*%tbeta_g
        #meanexpEta_gs[[g]]=meanexpEta_gs[[g]]+(N_gks[g,k]/N_gs[g])*meanexpEta_gks[[g,k]]
      }
      S_g=S_gs[[g]]
      Theta_gs[[g]]=phi_g-beta_g%*%(lambda_phi_g-S_g%*%tbeta_g)
    }

    Output_Mstep <- mixmgfa_intercepts_Mstep(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gs,Theta_gks,Theta_gs,meanexpEta_gks,Lambda,Psi_gs,Phi_gs,mean_gs,tau_ks,alpha_gks)
    Lambda=Output_Mstep$Lambda
    Psi_gs=Output_Mstep$Psi_gs
    Phi_gs=Output_Mstep$Phi_gs
    tau_ks=Output_Mstep$tau_ks
    alpha_gks=Output_Mstep$alpha_gks
    Sigma_gs=Output_Mstep$Sigma_gs
    invSigma_gs=Output_Mstep$invSigma_gs
    nractivatedconstraints=Output_Mstep$nractivatedconstraints

    S_gks <- matrix(list(NA),nrow = ngroup, ncol = nclust)
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    tLambda=t(Lambda)
    for(g in 1:ngroup){
      S_g=matrix(0,nvar,nvar)
      mean_g=matrix(mean_gs[g,],nrow=1,ncol=nvar)
      for(k in 1:nclust){
        mu_gk=tau_ks[k,]+alpha_gks[[g,k]]%*%tLambda
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
      logdet_sigma_g=log(det(Sigma_gs[[g]]))
      invSigma_g=invSigma_gs[[g]]
      for(k in 1:nclust){
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_g+sum(S_gks[[g,k]]*invSigma_g)) # sum(S_gks[[g,k]]*invSigma_g)=sum(diag(S_gks[[g,k]]%*%invSigma_g))
        loglik_gks[g,k]=loglik_gk
        loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
      }
      m_i=max(loglik_gksw[g,]);
      for(k in 1:nclust){
        loglik_gksw[g,k]=exp(loglik_gksw[g,k]-m_i); # exp because we have to sum over clusters before we can take the log
      }
      ODLL=ODLL+log(sum(loglik_gksw[g,]))+m_i;
    }

    pars=c(pi_ks,Lambda[design==1],unlist(lapply(Psi_gs,diag)),unlist(Phi_gs),as.vector(tau_ks),unlist(alpha_gks))
    parsdiff <- mapply("-",pars,prev_pars)
    parsdiffdiv <- mapply("/",parsdiff,prev_pars)
    parsabsdiffdiv <- abs(parsdiffdiv)
    conv1 <- sum(parsabsdiffdiv)
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


  # set scale of factors across groups (and, in case of EFA, make them orthogonal)
  theta=matrix(0,nfactors,nfactors)
  for(g in 1:ngroup){
    theta_g=Theta_gs[[g]]
    theta=theta+(N_gs[g]/N)*theta_g;
  }
  if(nfactors>1){
    if(EFA==1){
      # find matrix square root via eigenvalue decomposition
      ed=eigen(theta)
      sqrtFscale=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
      invsqrtFscale=solve(sqrtFscale)
    }
    else {
      sqrtFscale=diag(diag(theta^(1/2)))
      invsqrtFscale=diag(diag((1/theta)^(1/2)))
    }
  }
  else {
    sqrtFscale=sqrt(theta)
    invsqrtFscale=1/sqrtFscale
  }

  for(g in 1:ngroup){
    phi_g=Phi_gs[[g]]
    phi_g=invsqrtFscale%*%phi_g%*%invsqrtFscale;
    Phi_gs[[g]]=((phi_g+t(phi_g))/2); # enforce perfect symmetry
  }
  Lambda=Lambda%*%sqrtFscale # compensate for (re)scaling of factors in the loadings
  for(g in 1:ngroup){ # transform the factor means accordingly
    for(k in 1:nclust){
      alpha_gks[[g,k]]=alpha_gks[[g,k]]%*%invsqrtFscale
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
  # Translation of factor means per cluster and cycle 1 update indicator intercepts
  for(k in 1:nclust){
    for(g in 1:ngroup){
      alpha_gks[[g,k]]=alpha_gks[[g,k]]-mean_alpha_ks[k,]
    }
    if(N_ks[k]>0){
      suminvSigma=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvSigma=matrix(0,1,nvar)
      for(g in 1:ngroup){
        invSigma_g=invSigma_gs[[g]]
        summeansminusalphaLambdainvSigma=summeansminusalphaLambdainvSigma+N_gks[g,k]*(mean_gs[g,]-tcrossprod(alpha_gks[[g,k]],Lambda))%*%invSigma_g
        suminvSigma=suminvSigma+N_gks[g,k]*invSigma_g
      }
      tau_ks[k,]=t(solve(suminvSigma,t(summeansminusalphaLambdainvSigma)))
    }
  }



  if(EFA==1){
    nrpars=nclust-1+(nvar*nfactors-(nfactors*(nfactors-1)*(1/2)))+(nfactors*(nfactors+1)/2)*(ngroup-1)+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }
  else {
    nrpars=nclust-1+sum(design)+(nfactors*(nfactors+1)/2)*ngroup-nfactors+nvar*nclust+nfactors*(ngroup-nclust)+nvar*ngroup-nractivatedconstraints;
  }


  output_list <- list(z_gks=z_gks,pi_ks=pi_ks,Lambda=Lambda,Psi_gs=Psi_gs,Phi_gs=Phi_gs,tau_ks=tau_ks,alpha_gks=alpha_gks,bestloglik=bestloglik,logliks=logliks,nrpars=nrpars,convergence=convergence,nractivatedconstraints=nractivatedconstraints)

  return(output_list)
} # end main function



# functions for posterior classification probabilities (E-step) and cycle 2 M-step


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



mixmgfa_intercepts_Mstep <- function(S_gs,S_gks,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gs,Theta_gks,Theta_gs,meanexpEta_gks,Lambda,Psi_gs,Phi_gs,mean_gs,tau_ks,alpha_gks){
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
  tLambda=t(Lambda)
  for(k in 1:nclust){
    if(N_ks[k]>0){
      suminvPsi=matrix(0,nvar,nvar)
      summeansminusalphaLambdainvPsi=matrix(0,1,nvar)
      for(g in 1:ngroup){
        if(N_gks[g,k]>0){
          invPsi_g=invPsi_gs[[g]]
          summeansminusalphaLambdainvPsi=summeansminusalphaLambdainvPsi+N_gks[g,k]*(mean_gs[g,]-(alpha_gks[[g,k]]+meanexpEta_gks[[g,k]])%*%tLambda)%*%invPsi_g
          suminvPsi=suminvPsi+N_gks[g,k]*invPsi_g
        }
      }
      tau_ks[k,]=t(solve(suminvPsi,t(summeansminusalphaLambdainvPsi)))
    }
  }

  # update factor means
  for(g in 1:ngroup){
    invPsi_g=invPsi_gs[[g]]
    tLambdainvPsi_g=tLambda%*%invPsi_g
    for(k in 1:nclust){
      if(N_ks[k]>0){
        alpha_gks[[g,k]]=(mean_gs[g,]-tau_ks[k,]-meanexpEta_gks[[g,k]]%*%tLambda)%*%t(solve(((tLambdainvPsi_g)%*%Lambda),(tLambdainvPsi_g)))
      }
    }
  }

  # update factor loadings
  for(j in 1:nvar){
    nfactors_j=sum(design[j,])
    d_j=design[j,]==1
    sumSbeta=matrix(0,1,nfactors_j)
    sumthetaalpha=matrix(0,nfactors_j,nfactors_j)
    summeansalpha=matrix(0,1,nfactors_j)
    for(g in 1:ngroup){
      psi_g=Psi_gs[[g]]
      S_g=S_gs[[g]]
      beta_g=Beta_gs[[g]]
      if(nfactors_j<nfactors){
        beta_g=beta_g[d_j, ,drop=FALSE]
      }
      sumSbeta=sumSbeta+(N_gs[g]/psi_g[j,j])*tcrossprod(S_g[j,],beta_g)
      for(k in 1:nclust){
        if(N_gks[g,k]>0){
          theta_gk=Theta_gks[[g,k]]
          alpha_gk=alpha_gks[[g,k]]
          meanexpeta_gk=meanexpEta_gks[[g,k]]
          theta_gk=theta_gk[d_j,d_j]
          alpha_gk=alpha_gk[d_j]
          meanexpeta_gk=meanexpeta_gk[d_j]
          #talpha_gk=t(alpha_gk)
          sumthetaalpha=sumthetaalpha+(N_gks[g,k]/psi_g[j,j])*(theta_gk+tcrossprod(alpha_gk+meanexpeta_gk,alpha_gk))
          summeansalpha=summeansalpha+(N_gks[g,k]/psi_g[j,j])*((mean_gs[g,j]-tau_ks[k,j])%*%alpha_gk)
        }
      }
    }
    if(nfactors_j==1){
      Lambda[j,d_j]= (sumSbeta+summeansalpha)/sumthetaalpha
    } else {
      Lambda[j,d_j]= t(solve(sumthetaalpha,t(sumSbeta+summeansalpha)))
    }
  }
  tLambda=t(Lambda)

  # update unique variances
  nractivatedconstraints=0
  for(g in 1:ngroup){
    S_g=S_gs[[g]]
    beta_g=Beta_gs[[g]]
    theta_g=Theta_gs[[g]]
    twoSbetaB_BthetaB=Lambda%*%(2*beta_g%*%S_g-theta_g%*%tLambda) # modelimplied reduced covariance matrix on sample level, based on old structure matrix and sigma_gk, weighting based on new z_gks
    psi_g=diag(diag(S_g-twoSbetaB_BthetaB))
    if (sum(diag(psi_g)<.0001)>0){ # track "heywood" cases
      ind=diag(psi_g)<.0001
      d=diag(psi_g);
      d[ind]=0.0001;
      psi_g=diag(d);
      nractivatedconstraints=nractivatedconstraints+sum(ind)
    }
    Psi_gs[[g]]=psi_g
  }

  # update factor (co)variances
  for(g in 1:ngroup){
    theta_g=Theta_gs[[g]]
    phi_g=theta_g
    Phi_gs[[g]]=((phi_g+t(phi_g))/2); # enforce perfect symmetry to avoid accummulation of asymmetry over iterations
  }


  # update (inv)Sigma_gs
  Sigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
  invSigma_gs <- matrix(list(NA), nrow = ngroup, ncol = 1)
  for(g in 1:ngroup){
    psi_g=Psi_gs[[g]]
    invPsi_g=diag(1/diag(psi_g))
    invPsi_g_lambda=invPsi_g%*%Lambda
    phi_g=Phi_gs[[g]]
    invPhi_g=solve(phi_g)
    sigma_g=Lambda %*% phi_g %*% tLambda + psi_g
    Sigma_gs[[g]]=(sigma_g+t(sigma_g))*(1/2) # avoid asymmetry due to rounding errors
    invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g+tLambda%*%invPsi_g_lambda)
    invPhi_g_tLambdainvPsi_g_Lambda=(invPhi_g_tLambdainvPsi_g_Lambda+t(invPhi_g_tLambdainvPsi_g_Lambda))*(1/2)
    #invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda)%*%t(invPsi_g_lambda) # Woodbury identity
    invSigma_gs[[g]]=invPsi_g-invPsi_g_lambda%*%solve(invPhi_g_tLambdainvPsi_g_Lambda,tLambda)%*%invPsi_g # Woodbury identity
  }

  output_list <- list(Lambda=Lambda,Psi_gs=Psi_gs,Phi_gs=Phi_gs,tau_ks=tau_ks,alpha_gks=alpha_gks,Sigma_gs=Sigma_gs,invSigma_gs=invSigma_gs,nractivatedconstraints=nractivatedconstraints)

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


