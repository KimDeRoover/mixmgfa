#' Mixture multigroup factor analysis for loading (non-)invariance
# ----------------------------------------------------------------
# Code written by Kim De Roover

#' @description
#' Find clusters of groups based on their factor loadings, given a user-specified number of clusters.
#' Intercepts and residual variances are kept group-specific. Factor (co)variances are group- and cluster-specific. Factor means are zero. See referenced paper for details.
#'
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
#'          (Using different design matrices for different clusters is currently not supported.)
#' @param startpartition Partition of groups (vector) to start from (use with start = 2 and nruns = 1).

# OUTPUT:
#' @return Output object (list) with:
#'
#' $z_gks = cluster memberships of groups (posterior classification probabilities)
#'
#' $pi_ks = mixing proportions (prior classification probabilities)
#'
#' $Lambda_ks = cluster-specific loadings, access loadings of cluster k via Lambda_ks[[k]]
#'
#' $Psi_gs = group-specific unique variances, access loadings of group g via Psi_gs[[g]]
#'
#' $Phi_gks = group- and cluster-specific factor (co)variances, access (co)variances of group g in cluster k via Phi_gks[[g,k]]
#'
#' $mu_gs = group-specific means, access means of group g via mu_gs[g,]
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
#'
#' @references
#' De Roover, K., Vermunt, J. K., & Ceulemans, E. (2020). Mixture multigroup factor analysis for unraveling factor loading noninvariance across many groups. Psychological Methods. Advance online publication.
#' https://psycnet.apa.org/doi/10.1037/met0000355

#' @export
mixmgfa_loadings <- function(data,N_gs,nclust,nfactors=1,maxiter = 5000,start = 1,nruns = 50,design=0,preselect = 10,startpartition){

  ngroup <- length(N_gs)
  if(nrow(N_gs)!=ngroup || is.null(nrow(N_gs))){ # make sure N_gs is a column vector
    N_gs_colvec=matrix(0,ngroup,1)
    for (g in 1:ngroup){
      N_gs_colvec[g,]=N_gs[g]
    }
    N_gs <- N_gs_colvec
  }
  N <- sum(N_gs)
  IM <- diag(nclust)

  if (nrow(data)==N){ # input is raw data matrix
    nvar=ncol(data)
    Xsup=as.matrix(data)
    Ncum <- matrix(0,ngroup,2)
    Ncum[1,1]=1
    Ncum[1,2]=N_gs[1]
    for(g in 2:ngroup){
      Ncum[g,1]=sum(N_gs[1:(g-1)])+1 # Ncum[g,1]: first row in Xsup for group g
      Ncum[g,2]=sum(N_gs[1:g]) # Ncum[g,2]: last row in Xsup for group g
    }
    # compute group-specific means, work with group-centered data from now on
    mu_gs=matrix(0,ngroup,nvar)
    Xsupcent <- matrix(0,N,nvar)
    for(g in 1:ngroup){
      X_g <- Xsup[Ncum[g,1]:Ncum[g,2],]
      mu_gs[g,] <- colMeans(X_g)
      Xsupcent[Ncum[g,1]:Ncum[g,2],]=scale(X_g,scale = FALSE)
    }
    Xsup <- Xsupcent

    # compute sample covariance matrices
    S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
    for(g in 1:ngroup){
      X_g=Xsup[Ncum[g,1]:Ncum[g,2],]
      S_gs[[g]] <- (1/N_gs[g])*(t(X_g)%*%X_g)
    }
  } else { # input is list or concatenation of covariance matrices
    if (is.list(data)){ # input is a list of covariance matrices
      nvar=ncol(data[[1]])
      S_gs=data
      mu_gs=matrix(0,ngroup,nvar)
    }
    else { # input is concatenation of covariance matrices that should be turned into a list
      nvar=ncol(data)
      S_gs <- matrix(list(NA),nrow = ngroup, ncol = 1)
      for(g in 1:ngroup){
        S_gs[[g]] <- data[nvar*(g-1)+1:nvar*g,]
      }
      mu_gs=matrix(0,ngroup,nvar)
    }
  }

  if (sum(design)==0){ # if design is unspecified, EFA is used
    design=matrix(1,nvar,nfactors)
    EFA=1
  } else {
    EFA=0
  }


  if(start==1){
    if(nclust>1){
      # pre-selection of random partitions
      nrtrialstarts=nruns*round(100/preselect) # generate 'nruns'*(100/preselect) different random partitions
      randpartvecs=matrix(0,nrtrialstarts,ngroup);
      for (trialstart in 1:nrtrialstarts){
        aris=1;
        iterrnd=0;
        while (sum(aris==1)>0 && iterrnd<5){
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
            iterrnd=iterrnd+1;
          }
        }
        randpartvecs[trialstart,]=randpartvec
      }
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
        N_gks=diag(N_gs[,1])%*%z_gks
        N_ks=colSums(N_gks)

        Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
        uniq_ks <- matrix(0,nclust,nvar)
        for(k in 1:nclust){
          # Xsup_k=matrix(0,N_ks[k],nvar)
          # for(g in 1:ngroup){
          #   if(randpartvec[g]==k){
          #     X=Xsup[Ncum[g,1]:Ncum[g,2],]
          #     Xsup_k[(sum(N_gks[1:g-1,k])+1):sum(N_gks[1:g,k]),]=X
          #   }
          # }
          # S_k <- (1/sum(N_gs[randpartvec==k]))*(t(Xsup_k)%*%Xsup_k)
          S_k=matrix(0,nvar,nvar)
          for(g in 1:ngroup){
            if(randpartvec[g]==k){
              S_k=S_k+N_gs[g]*S_gs[[g]]
            }
          }
          S_k=(1/sum(N_gs[randpartvec==k]))*S_k
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
        for(g in 1:ngroup){
          k=randpartvec[g]
          Psi_gs[[g]]=diag(uniq_ks[k,])
          for(k2 in 1:nclust){
            Phi_gks[[g,k2]]=diag(nfactors)
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
            phi_gk=Phi_gks[[g,k]]
            invPhi_gk=phi_gk # here phi_gk is still identity matrix
            sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
            Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
            #EV <- eigen((invPhi_gk+t(lambda_k)%*%invPsi_g%*%lambda_k), only.values=TRUE, symmetric = TRUE)
            #print(EV)
            invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk+tlambda_k%*%invPsi_g%*%lambda_k)
            invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk_tlambda_kinvPsi_g_lambda_k+t(invPhi_gk_tlambda_kinvPsi_g_lambda_k))*(1/2)
            invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_gk_tlambda_kinvPsi_g_lambda_k,tlambda_k)%*%invPsi_g; # Woodbury identity
          }
        }

        # compute Beta_gks and theta_gks
        Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
        for(g in 1:ngroup){
          S_g=S_gs[[g]]
          for(k in 1:nclust){
            lambda_k=Lambda_ks[[k]]
            invsigma_gk=invSigma_gks[[g,k]]
            phi_gk=Phi_gks[[g,k]]
            beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
            Beta_gks[[g,k]]=beta_gk;
            theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_g%*%t(beta_gk)
            Theta_gks[[g,k]]=theta_gk
          }
        }


        Output_Mstep <- mixmgfa_load_Mstep(S_gs,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,Lambda_ks,Psi_gs,Phi_gks)
        Lambda_ks=Output_Mstep$Lambda_ks
        Psi_gs=Output_Mstep$Psi_gs
        Phi_gks=Output_Mstep$Phi_gks
        Sigma_gks=Output_Mstep$Sigma_gks
        invSigma_gks=Output_Mstep$invSigma_gks
        nractivatedconstraints=Output_Mstep$nractivatedconstrains


        # compute observed-data log-likelihood for start
        ODLL_trialstart=0;
        loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
        loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
        for(g in 1:ngroup){
          # X_g=Xsup[Ncum[g,1]:Ncum[g,2],] # centered data per group
          # for(k in 1:nclust){
          #   logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
          #   invSigma_gk=invSigma_gks[[g,k]]
          #   loglik_gk=-(1/2)*(N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk))
          #   for (n in 1:N_gs[g]){
          #     X_n=X_g[n, ,drop=FALSE]
          #     loglik_gk=loglik_gk-(1/2)*(X_n%*%tcrossprod(invSigma_gk,X_n))
          #   }
          #   loglik_gks[g,k]=loglik_gk
          #   loglik_gksw[g,k]=log(pi_ks[k])+loglik_gk
          # }
          for(k in 1:nclust){
            logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
            loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gs[[g]]*invSigma_gks[[g,k]])) # sum(S_gs[[g]]*invSigma_gks[[g,k]])=sum(sum(diag(S_gs[[g]]%*%invSigma_gks[[g,k]]))
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
      sortedstarts <- sort(ODLLs_trialstarts,decreasing=TRUE,index.return=TRUE);
      index_beststarts=sortedstarts$ix[1:nruns]
      if (nrtrialstarts>nruns){
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
      nruns=1;
      randpartvecs=matrix(1,1,ngroup)
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
    Lambda_ks <- matrix(list(NA),nrow = 1, ncol=nclust)
    uniq_ks <- matrix(0,nclust,nvar)
    for(k in 1:nclust){
      S_k=matrix(0,nvar,nvar)
      for(g in 1:ngroup){
        if(randpartvec[g]==k){
          S_k=S_k+N_gs[g]*S_gs[[g]]
        }
      }
      S_k=(1/sum(N_gs[randpartvec==k]))*S_k
      ed<-eigen(S_k, symmetric=TRUE, only.values = FALSE)
      val<-ed$values
      u<-ed$vectors
      totalerror=sum((val[-seq_len(nfactors)]))
      meanerror=totalerror/(nvar-nfactors) # mean error variance: mean variance in discarded dimensions
      Uniq=rep(meanerror,nvar)
      lambda_k=u[,seq_len(nfactors),drop=FALSE] %*% sqrt(diag(val[seq_len(nfactors)]-Uniq[seq_len(nfactors)],nrow=nfactors,ncol=nfactors))
      Lambda_ks[[k]]=lambda_k
      uniq_ks[k,]=Uniq
    }
    Psi_gs <- matrix(list(NA), nrow = ngroup, ncol = 1) # initialize group-specific unique variances
    Phi_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust) # initialize group- and cluster-specific factor covariances
    for(g in 1:ngroup){
      k=randpartvec[g]
      Psi_gs[[g]]=diag(uniq_ks[k,])
      for(k2 in 1:nclust){
        Phi_gks[[g,k2]]=diag(nfactors)
      }
    }

    # initialize prior and posterior classification probabilities
    if(nclust>1){
      z_gks=IM[randpartvec,]
      pi_ks=(1/ngroup)*colSums(z_gks)
    }
    else {
      z_gks=t(randpartvec)
      pi_ks=1
    }


    Sigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    invSigma_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    for(k in 1:nclust){
      lambda_k=Lambda_ks[[k]]
      tlambda_k=t(lambda_k)
      for(g in 1:ngroup){
        psi_g=Psi_gs[[g]]
        invPsi_g=diag(1/diag(psi_g))
        phi_gk=Phi_gks[[g,k]]
        invPhi_gk=solve(phi_gk)
        sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
        Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
        invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk+tlambda_k%*%invPsi_g%*%lambda_k)
        invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk_tlambda_kinvPsi_g_lambda_k+t(invPhi_gk_tlambda_kinvPsi_g_lambda_k))*(1/2)
        invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_gk_tlambda_kinvPsi_g_lambda_k,tlambda_k)%*%invPsi_g; # Woodbury identity
      }
    }

    # compute the loglikelihood for each group-cluster combination, unweighted with mixing proportions, to be used for update of posterior classification probabilities
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust);
    for(g in 1:ngroup){
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gs[[g]]*invSigma_gks[[g,k]])) # sum(S_gs[[g]]*invSigma_gks[[g,k]])=sum(sum(diag(S_gs[[g]]%*%invSigma_gks[[g,k]]))
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
    pars=c(pars,unlist(lapply(Psi_gs,diag)),unlist(Phi_gks))
    while(min(conv1,conv2)>1e-4 && iter<101){
      prev_ODLL=ODLL;
      prev_Lambda_ks=Lambda_ks;
      prev_Psi_gs=Psi_gs;
      prev_Phi_gks=Phi_gks;
      prev_pars=pars;
      iter=iter+1

      # **E-step**: compute the posterior classification probabilities
      z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust, nfactors)


      N_gks=diag(N_gs[,1])%*%z_gks
      N_ks=colSums(N_gks)

      # update mixing proportions
      pi_ks=(1/ngroup)*colSums(z_gks)

      # compute Beta_gks and theta_gks
      Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
      for(g in 1:ngroup){
        S_g=S_gs[[g]]
        for(k in 1:nclust){
          lambda_k=Lambda_ks[[k]]
          invsigma_gk=invSigma_gks[[g,k]]
          phi_gk=Phi_gks[[g,k]]
          beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
          Beta_gks[[g,k]]=beta_gk;
          theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_g%*%t(beta_gk)
          Theta_gks[[g,k]]=theta_gk
        }
      }

      Output_Mstep <- mixmgfa_load_Mstep(S_gs,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,Lambda_ks,Psi_gs,Phi_gks)
      Lambda_ks=Output_Mstep$Lambda_ks
      Psi_gs=Output_Mstep$Psi_gs
      Phi_gks=Output_Mstep$Phi_gks
      Sigma_gks=Output_Mstep$Sigma_gks
      invSigma_gks=Output_Mstep$invSigma_gks
      nractivatedconstraints=Output_Mstep$nractivatedconstraints

      # check on change in observed-data log-likelihood
      ODLL=0;
      loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
      loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
      for(g in 1:ngroup){
        for(k in 1:nclust){
          logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
          loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gs[[g]]*invSigma_gks[[g,k]])) # sum(S_gs[[g]]*invSigma_gks[[g,k]])=sum(sum(diag(S_gs[[g]]%*%invSigma_gks[[g,k]]))
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
      pars=c(pars,unlist(lapply(Psi_gs,diag)),unlist(Phi_gks))
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
    # hier parameters beste start tot nu toe bijhouden
    logliks[run,]=c(ODLL,nractivatedconstraints);
    if (run==1) {
      bestz_gks=z_gks
      bestpi_ks=pi_ks
      bestLambda_ks=Lambda_ks
      bestPsi_gs=Psi_gs;
      bestPhi_gks=Phi_gks
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
        bestPsi_gs=Psi_gs;
        bestPhi_gks=Phi_gks
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
  pars=c(pars,unlist(lapply(Psi_gs,diag)),unlist(Phi_gks))
  while(min(conv1,conv2)>1e-6 && iter<maxiter+1){ # iterate till convergence for best start
    prev_ODLL=ODLL
    prev_Lambda_ks=Lambda_ks
    prev_Psi_gs=Psi_gs
    prev_Phi_gks=Phi_gks
    prev_pars=pars
    iter=iter+1

    # **E-step**: compute the posterior classification probabilities
    z_gks <- UpdPostProb(pi_ks, loglik_gks, ngroup, nclust, nfactors)


    N_gks=diag(N_gs[,1])%*%z_gks
    N_ks=colSums(N_gks)

    # update mixing proportions
    pi_ks=(1/ngroup)*colSums(z_gks)

    # compute Beta_gks and theta_gks
    Beta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    Theta_gks <- matrix(list(NA), nrow = ngroup, ncol = nclust)
    for(g in 1:ngroup){
      S_g=S_gs[[g]]
      for(k in 1:nclust){
        lambda_k=Lambda_ks[[k]]
        invsigma_gk=invSigma_gks[[g,k]]
        phi_gk=Phi_gks[[g,k]]
        beta_gk=phi_gk%*%t(lambda_k)%*%invsigma_gk
        Beta_gks[[g,k]]=beta_gk;
        theta_gk=phi_gk-beta_gk%*%lambda_k%*%phi_gk+beta_gk%*%S_g%*%t(beta_gk)
        Theta_gks[[g,k]]=theta_gk
      }
    }

    Output_Mstep <- mixmgfa_load_Mstep(S_gs,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,Lambda_ks,Psi_gs,Phi_gks)
    Lambda_ks=Output_Mstep$Lambda_ks
    Psi_gs=Output_Mstep$Psi_gs
    Phi_gks=Output_Mstep$Phi_gks
    Sigma_gks=Output_Mstep$Sigma_gks
    invSigma_gks=Output_Mstep$invSigma_gks
    nractivatedconstraints=Output_Mstep$nractivatedconstraints

    # check on change in observed-data log-likelihood
    ODLL=0;
    loglik_gks <- matrix(0, nrow = ngroup, ncol = nclust) # unweighted with mixing proportions, to be re-used for calculation posterior classification probabilities
    loglik_gksw <- matrix(0, nrow = ngroup, ncol = nclust) # weighted with mixing proportions
    for(g in 1:ngroup){
      for(k in 1:nclust){
        logdet_sigma_gk=log(det(Sigma_gks[[g,k]]))
        loglik_gk=-(1/2)*N_gs[g]*(nvar*log(2*pi)+logdet_sigma_gk+sum(S_gs[[g]]*invSigma_gks[[g,k]])) # sum(S_gs[[g]]*invSigma_gks[[g,k]])=sum(sum(diag(S_gs[[g]]%*%invSigma_gks[[g,k]]))
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
    pars=c(pars,unlist(lapply(Psi_gs,diag)),unlist(Phi_gks))
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


  # set scale of factors per cluster and, in case of EFA, make them orthogonal
  for(k in 1:nclust){
    if(N_ks[k]>1e-8){
      theta_k=matrix(0,nfactors,nfactors)
      for(g in 1:ngroup){
        if(N_gks[g,k]>0){
          theta_gk=Theta_gks[[g,k]]
          theta_k=theta_k+(N_gks[g,k]/N_ks[k])*theta_gk;
        }
      }
      if(nfactors>1){
        if(EFA==1){
          # find matrix square root via eigenvalue decomposition
          ed=eigen(theta_k)
          sqrtheta_k=ed$vectors%*%diag(ed$values)^(1/2)%*%solve(ed$vectors)
          invsqrtheta_k=solve(sqrtheta_k)
        }
        else {
          sqrtheta_k=diag(diag(theta_k^(1/2)))
          invsqrtheta_k=diag(diag((1/theta_k)^(1/2)))
        }
      }
      else {
        sqrtheta_k=sqrt(theta_k)
        invsqrtheta_k=1/sqrtheta_k
      }
      for(g in 1:ngroup){
        phi_gk=Phi_gks[[g,k]]
        phi_gk=invsqrtheta_k%*%phi_gk%*%invsqrtheta_k; # rescale to satisfy restriction on phi_k=sum((N_ik/N_i)*phi_ik)=I
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry
      }
      Lambda_ks[[k]]=Lambda_ks[[k]]%*%sqrtheta_k # compensate for (re)scaling of factors per clusters in the loadings
    }
  }


  if(EFA==1){
    nrpars=nclust-1+(nvar*nfactors-(nfactors*(nfactors-1)*(1/2)))*nclust+(nfactors*(nfactors+1)/2)*(ngroup-nclust)+nvar*ngroup*2-nractivatedconstraints;
  }
  else {
    nrpars=nclust-1+sum(design)*nclust+(nfactors*(nfactors+1)/2)*ngroup-nclust*nfactors+nvar*ngroup*2-nractivatedconstraints;
  }



  output_list <- list(z_gks=z_gks,pi_ks=pi_ks,Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,mu_gs=mu_gs,bestloglik=bestloglik,logliks=logliks,nrpars=nrpars,convergence=convergence,nractivatedconstraints=nractivatedconstraints)

  return(output_list)
} # end main function



# functions for E-step/posterior classification probabilities and M-step


# Update the cluster-membership probabilities z_gk
# Reuses loglik_gks to save time

UpdPostProb <- function(pi_ks, loglik_gks, ngroup, nclust, nfact){
  max_g <-rep(0,ngroup)
  z_gks <- matrix(NA,nrow=ngroup,ncol=nclust)

  for(g in 1:ngroup){
    for(k in 1:nclust){
      z_gks[g,k] <- log(pi_ks[k])+loglik_gks[g,k]
    }
    max_g[g] <- max(z_gks[g,]) # prevent arithmetic underflow
    z_gks[g,] <- exp(z_gks[g,]-rep(max_g[g],nclust))
  }

  # divide by the rowwise sum of the above calculated part
  z_gks <- diag(1/rowSums(z_gks))%*%z_gks
  z_gks <- round(z_gks,digits=16)
  z_gks <- diag(1/rowSums(z_gks))%*%z_gks

  return(z_gks)
}



mixmgfa_load_Mstep <- function(S_gs,N_gs,nvar,nclust,nfactors,design,N_gks,Beta_gks,Theta_gks,Lambda_ks,Psi_gs,Phi_gks){
  nractivatedconstraints <- 0
  ngroup <- length(N_gs)
  N_ks=colSums(N_gks)

  # update cluster-specific loadings
  for(k in 1:nclust){
    if(N_ks[k]>1e-6){
      lambda_k=matrix(0,nvar,nfactors)
      for(j in 1:nvar){
        nfactors_j=sum(design[j,])
        d_j=design[j,]==1
        sumSbeta_k=matrix(0,1,nfactors_j)
        sumtheta_k=matrix(0,nfactors_j,nfactors_j)
        for(g in 1:ngroup){
          if(N_gks[g,k]>0){
            psi_g=Psi_gs[[g]]
            S_g=S_gs[[g]]
            beta_gk=Beta_gks[[g,k]]
            beta_gk=beta_gk[d_j, ,drop=FALSE]
            theta_gk=Theta_gks[[g,k]]
            theta_gk=theta_gk[d_j,d_j]
            sumSbeta_k=sumSbeta_k+(N_gks[g,k]/psi_g[j,j])*S_g[j,]%*%t(beta_gk)
            sumtheta_k=sumtheta_k+(N_gks[g,k]/psi_g[j,j])*theta_gk
          }
        }
        lambda_k[j,d_j]= t(solve(sumtheta_k,t(sumSbeta_k)))
      }
      Lambda_ks[[k]]=lambda_k
    }
  }


  # update unique variances
  nractivatedconstraints=0
  for(g in 1:ngroup){
    S_g=S_gs[[g]]
    sum2SbetaB_BthetaB=0;
    for(k in 1:nclust){
      if(N_gks[g,k]>0){
        lambda_k=Lambda_ks[[k]]
        beta_gk=Beta_gks[[g,k]]
        theta_gk=Theta_gks[[g,k]]
        sum2SbetaB_BthetaB=sum2SbetaB_BthetaB+(N_gks[g,k]/N_gs[g])*(2*lambda_k%*%beta_gk%*%S_g-lambda_k%*%theta_gk%*%t(lambda_k)) # modelimplied reduced covariance matrix on sample level, based on old structure matrix and sigma_gk, weighting based on new z_gks
      }
    }
    psi_g=diag(diag(S_g-sum2SbetaB_BthetaB))
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
  for(k in 1:nclust){
    if(N_ks[k]>1e-8){
      # theta_k=matrix(0,nfactors,nfactors)
      # for(g in 1:ngroup){
      #   if(N_gks[g,k]>0){
      #     theta_gk=Theta_gks[[g,k]]
      #     theta_k=theta_k+(N_gks[g,k]/N_ks[k])*theta_gk;
      #   }
      # }
      #invsqrtheta_k=Re(solve(theta_k)^(1/2))
      for(g in 1:ngroup){
        theta_gk=Theta_gks[[g,k]]
        phi_gk=theta_gk#invsqrtheta_k%*%theta_gk%*%invsqrtheta_k; # rescale to satisfy restriction on phi_k=sum((N_ik/N_i)*phi_ik)=I
        Phi_gks[[g,k]]=((phi_gk+t(phi_gk))/2); # enforce perfect symmetry to avoid accummulation of asymmetry over iterations
      }
      # Lambda_ks[[k]]=Lambda_ks[[k]]%*%Re(theta_k^(1/2)) # compensate for (re)scaling of factors per clusters in the loadings
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
      phi_gk=Phi_gks[[g,k]]
      invPhi_gk=solve(phi_gk)
      sigma_gk=lambda_k %*% phi_gk %*% tlambda_k + psi_g
      Sigma_gks[[g,k]]=(sigma_gk+t(sigma_gk))*(1/2) # avoid asymmetry due to rounding errors
      invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk+tlambda_k%*%invPsi_g%*%lambda_k)
      invPhi_gk_tlambda_kinvPsi_g_lambda_k=(invPhi_gk_tlambda_kinvPsi_g_lambda_k+t(invPhi_gk_tlambda_kinvPsi_g_lambda_k))*(1/2)
      invSigma_gks[[g,k]]=invPsi_g-invPsi_g%*%lambda_k%*%solve(invPhi_gk_tlambda_kinvPsi_g_lambda_k,tlambda_k)%*%invPsi_g; # Woodbury identity
    }
  }

  output_list <- list(Lambda_ks=Lambda_ks,Psi_gs=Psi_gs,Phi_gks=Phi_gks,Sigma_gks=Sigma_gks,invSigma_gks=invSigma_gks,nractivatedconstraints=nractivatedconstraints)

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
