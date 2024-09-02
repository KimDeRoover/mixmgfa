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
#' @param cluster.spec Measurement parameters you want to cluster the groups on; "loadings", "intercepts", "residuals", c("loadings","intercepts"), c("intercepts","residuals"), c("loadings","residuals"), or c("loadings","intercepts","residuals").
#' Note: cluster.spec = "intercepts" and cluster.spec = c("intercepts","residuals") impose invariant loadings across all groups, cluster.spec = "residuals" also imposes invariant intercepts across all groups.
#' @param nsclust Vector of length one or two indicating the number(s) of clusters. In case of length two, the vector indicates the minimal and maximal number of clusters (it is recommended to set the minimal number to one).
#' @param maxiter Maximum number of iterations used in each MMG-FA analysis. Increase in case of non-convergence.
#' @param nruns Number of (preselected) random starts (important for avoiding local maxima in case of few groups and/or small groups).
#' @param design For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).
#'          (Using different design matrices for different clusters is currently not supported.)
#' @param rotation Rotation criterion to use in case of EFA; currently either "oblimin", "geomin", "varimax" or "target" (i.e., semi-specified oblique Procrustes rotation), whereas 0 = no rotation. (Note: The GPArotation package is loaded or installed for rotation.)
#' @param preselect Percentage of best starts taken in pre-selection of initial partitions (for huge datasets, increase to speed up multistart procedure).
#' @param targetT Target matrix to use when rotation = "target". Note that the same target is used for all clusters. For using cluster-specific target matrices, use rotation = 0 and rotate the cluster-specific factors afterwards with the GPFobl function of GPArotation.
#' @param targetW Weights to be used when rotation = "target". You can set an entire row to zero to make sure that the simple structure is not messed up by a 'complex variable' (i.e., with strong loadings for multiple factors). Set all weights equal to '1' if you prefer fully specified target rotation. When left unspecified while rotation = target, the default is a targetW where the zeros in the target get a weight of '1' and the non-zeros in the target get a weight of '0'. If this results in all zero weights or too many zeros for the rotation to be identified, a fully specified target rotation (all weights equal to '1') is used instead.
#' @param rescov Binary matrix to specify residual item covariances. Not yet supported in the current package version!
#' @param parcomp Specified as 0 or 1 to indicate if parallel computation should be performed to speed up the computation. If parcomp==1, analyses with different numbers of clusters are executed simultaneously by using multiple cores on your machine.
#' @param freecores The number of cores kept free for other tasks during parallel computation (if parcomp == 1).


# OUTPUT:
#' @return Output object (list) with:
#'
#' $overview = overview of fitted MMG-FA solutions with loglikelihood (loglik), number of parameters (nrpars), BIC_N (using total number of observations as sample size), BIC_G (using number of groups as sample size), AIC,
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
#' Leitgöb, H., Seddig, D., Asparouhov, T., Behr, D., Davidov, E., De Roover, K., ... & van de Schoot, R. (2023). Measurement invariance in the social sciences: Historical development, methodological challenges, state of the art, and future perspectives. Social Science Research, 110, 102805.

#' @export
mixmgfa <- function(data,N_gs=c(),nfactors=1, cluster.spec = c("loadings","intercepts","residuals"),nsclust = c(1,5),maxiter = 5000,nruns = 25,design=0,rotation=0,preselect = 10,targetT=0,targetW=0,rescov=0,parcomp=0,freecores=2){

  allowed=c("loadings","intercepts","residuals")
  if(sum(is.element(cluster.spec,allowed)==FALSE)>0){
    stop("The specification of cluster.spec seems to contain an error.")
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

  if(max(design)!=0){
    if(ncol(design)!=nfactors){
      stop("The specified design matrix does not match the specified number of factors.")
    }
  }

  if(parcomp==1){
    if(!require(doParallel)){ # test if package is loaded
      loadtest<-try(library(doParallel),silent=TRUE) # test if package can be loaded (if it is installed)
      if(class(loadtest)=="try-error"){
        install.packages("doParallel")
      }
    }
    if(!require(parallel)){ # test if package is loaded
      loadtest<-try(library(parallel),silent=TRUE) # test if package can be loaded (if it is installed)
      if(class(loadtest)=="try-error"){
        install.packages("parallel")
      }
    }
  }


  options(warn=-1)

  if (is.matrix(data)==TRUE | is.data.frame(data)==TRUE){ # input is raw data matrix
    if (is.null(N_gs) || !exists("N_gs")){ # create N_gs if not given as input
      T=table(data[,1]) #N_gs<-tabulate(data[,1])
      N_gs=as.numeric(T)
      grouplabels=names(T)
      N_gs=N_gs[N_gs!=0]
      N_gs=as.matrix(N_gs)
      checkfirstcolumn=0
      firstcolID=TRUE
      # Note: without N_gs and without a first column with group ID, the first variable is taken as the group ID so that the N_gs created above matches this column
    } else { # if N_gs is given as input, check whether first column of data contains group ID (below)
      checkfirstcolumn=1
      firstcolID=FALSE # will become TRUE below if tabulation of first column matches N_gs
    }
    ngroup <- length(N_gs)
    if(nrow(N_gs)!=ngroup || is.null(nrow(N_gs))){ # make sure N_gs is a column vector
      N_gs_colvec=matrix(0,ngroup,1)
      for (g in 1:ngroup){
        N_gs_colvec[g,]=N_gs[g]
      }
      N_gs <- N_gs_colvec
    }
    if(sum(N_gs)!=nrow(data)){ # N_gs does not match the data size
      stop("The analysis cannot be performed because the specified N_gs vector does not match the data size.")
    }
    if (checkfirstcolumn==1){
      if(sum(diff(data[,1])!=0)>ngroup){
        stop("The analysis cannot be performed because the data matrix is not ordered according to the group labels in column one. All rows of a group should be placed directly below one another.")
      }
      T<-table(data[,1])
      checkcol1=as.numeric(T)
      checkcol1=checkcol1[checkcol1!=0]
      checkcol1=matrix(checkcol1,nrow=length(checkcol1),ncol=1)
      firstcolID=all.equal(N_gs,checkcol1)
      firstcolID=(firstcolID[1]==TRUE)
    }

    if (firstcolID){
      nvar <- ncol(data)-1
      Xsup=as.matrix(data[,2:(nvar+1)])
      grouplabels=unique(data[,1])
      varlabels=colnames(data[,2:(nvar+1)])
    } else {
      cat("Note: The data should be ordered according to the group memberships (with all rows of a group placed directly below one another). Reorder the data and repeat the analysis if this is not the case.")
      cat("\n")
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
    # if(length(cluster.spec)==1 && cluster.spec=="loadings"){
    #   data2=obsS_gs
    # } else if (length(cluster.spec)==2 && is.element("loadings",cluster.spec) && is.element("residuals",cluster.spec)) {
    #   data2=obsS_gs
    # } else {
      data2=list(covariances=obsS_gs,means=mean_gs)
    # }
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
    # if(length(cluster.spec)==1 && cluster.spec[1]=="loadings"){
    #   data2=obsS_gs
    # } else if (length(cluster.spec)==2 && is.element("loadings",cluster.spec) && is.element("residuals",cluster.spec)) {
    #   data2=obsS_gs
    # } else {
      data2=list(covariances=obsS_gs,means=mean_gs)
    # }
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


  if(sum(rescov)>nvar){
    if(length(cluster.spec)==2 && is.element("loadings",cluster.spec) && is.element("residuals",cluster.spec)){
      # check specification of residual covariances (symmetry and diagonal)
      rescov=ceiling((rescov+t(rescov))/2)
      diag(rescov)=rep(1,nvar)
    } else {
      # residual covariance not (yet) supported
      rescov=0
    }
  }

  if(parcomp==1){
    # Set up the parallel processing:
    # Detect the amount of cores on the machine
    cores <- detectCores()

    # Create a cluster that links the different cores, to perform
    # tasks in parallel (in the different cores at the same time),
    # the cluster ensures the different tasks can communicate with each other.
    # at least 2 (the number specified by 'freecores') of your total cores are not used to ensure system stability and
    # being able to still do other processes on your machine
    cl <- makeCluster(cores - freecores)

    # Activate the cluster to be the backend of the parallel
    # operations of the foreach function
    registerDoParallel(cl)
  }

  MMGFAsolutions=matrix(list(NA),nrow = 1, ncol=nsclust[2]-nsclust[1]+1)
  if (length(cluster.spec)==1 && cluster.spec[1]=="loadings"){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_loadings(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]]<- do_mixmgfa_loadings(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (length(cluster.spec)==1 && cluster.spec=="intercepts"){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_intercepts(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_intercepts(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (length(cluster.spec)==1 && cluster.spec=="residuals"){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_residuals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_residuals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (is.element("loadings",cluster.spec) & is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec)==FALSE){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_loadingsintercepts(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_loadingsintercepts(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (is.element("loadings",cluster.spec) & is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec)){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_loadingsinterceptsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_loadingsinterceptsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (is.element("intercepts",cluster.spec) & is.element("residuals",cluster.spec) & is.element("loadings",cluster.spec)==FALSE){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_interceptsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,parcomp=parcomp)
      }
    }
  }
  if (is.element("loadings",cluster.spec) & is.element("residuals",cluster.spec) & is.element("intercepts",cluster.spec)==FALSE){
    if(parcomp==1){
      MMGFAsolutions <- foreach(nclust = nsclust[1]:nsclust[2], .packages=c("mixmgfa","GPArotation")) %dopar% {
        outputnclust2 <- do_mixmgfa_loadingsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,rescov=rescov,parcomp=parcomp)
        return(outputnclust2)
      }
      stopCluster(cl)
    } else {
      for(nclust in nsclust[1]:nsclust[2]){
        MMGFAsolutions[[nclust-nsclust[1]+1]] <- do_mixmgfa_loadingsresiduals(data2,N_gs,nclust,nfactors,varlabels,grouplabels,maxiter = maxiter,start = 1,nruns = nruns,design=design,rotation=rotation,preselect=preselect,targetT=targetT,targetW=targetW,rescov=rescov,parcomp=parcomp)
      }
    }
  }

  overview=matrix(0,nsclust[2]-nsclust[1]+1,8)
  for(nclust in nsclust[1]:nsclust[2]){
    output_nclust2=MMGFAsolutions[[nclust-nsclust[1]+1]]
    loglik=output_nclust2$loglik
    nrpars=output_nclust2$nrpars
    overview[nclust-nsclust[1]+1,]=cbind(nclust,loglik,nrpars,-2*loglik+nrpars*log(N),-2*loglik+nrpars*log(ngroup),-2*loglik+nrpars*2,output_nclust2$convergence,output_nclust2$nractivatedconstraints)
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
    prefix=seq(nsclust[1],nsclust[2])
    suffix="clusters"
    sollistnames<-c(paste(prefix,suffix,sep="."))
    names(MMGFAsolutions)<-noquote(sollistnames)
    #overview=as.data.frame(overview,row.names=FALSE)
    if(nrows>2){
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","AIC","screeratios","convergence","nr.activated.constraints")
    } else {
      colnames(overview)<-c("nr of clusters","loglik","nrpars","BIC_N","BIC_G","AIC","convergence","nr.activated.constraints")
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
    if(nrows>2){
      cat("Choose the best number of clusters ('K_best') based on AIC, BIC_G and the CHull scree ratios. Also look at the plots obtained by using 'plot(OutputObject$overview)'.")
      cat("\n")
      cat("Based on the BIC_G, look for the number of clusters that minimizes the BIC_G or that corresponds to an elbow point in the BIC_G plot (after which the decrease with additional clusters levels off).")
      cat("\n")
      cat("Based on the CHull (scree ratios AND plot), look for the number of clusters that has the maximal scree ratio AND check whether this corresponds to at least a mild elbow point in the lower plot.")
      cat("\n")
    } else {
      cat("Choose the best number of clusters ('K_best') based on the AIC and BIC_G.")
      cat("\n")
      cat("Based on the AIC/BIC_G, look for the number of clusters that minimizes the AIC/BIC_G.")
      cat("\n")
      cat("The CHull could not be performed, because too few numbers of clusters were considered.")
      cat("\n")
    }

    cat("Access the corresponding cluster memberships and parameter estimates by using OutputObject$MMGFAsolutions[[K_best]]$clustermemberships and, for example, OutputObject$MMGFAsolutions[[K_best]]$clusterspecific.loadings.")
    cat("\n")
    cat("The parameter sets are further subdivided in group- and/or cluster-specific parameter sets.")
    cat("\n")
    class(MMGFAsolutions)<-"mixmgfa"
    class(overview)<-"mixmgfaMS"
    output<-list(overview=overview,MMGFAsolutions=MMGFAsolutions)
    #class(output)<-"mixmgfa"

  } else {
    cat("No output is generated. Perhaps you made a typo in 'cluster.spec' (make sure you use lowercase letters) or requested a model that is not supported.")
    cat("\n")
    cat("Please try again.")
    cat("\n")
  }
}
