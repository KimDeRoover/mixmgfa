% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScaleRotateMixmgfa.R
\name{ScaleRotateMixmgfa}
\alias{ScaleRotateMixmgfa}
\title{Re-scaling and rotation function for output objects of mixture multigroup factor analysis}
\usage{
ScaleRotateMixmgfa(
  OutputObject,
  N_gs,
  cluster.spec,
  nsclust = c(),
  design = 0,
  rescale = 0,
  markers = 0,
  rotation = 0,
  targetT = 0,
  targetW = 0
)
}
\arguments{
\item{OutputObject}{Output object resulting from using the mixmgfa function, either the entire output object or the MMGFAsolutions part. This may be rotated or unrotated.}

\item{N_gs}{Vector with number of subjects (sample size) for each group (in the same order as they appear in the data).}

\item{cluster.spec}{Measurement parameters you clustered the groups on: "loadings", "intercepts", "residuals", c("loadings","intercepts"), c("intercepts","residuals"), c("loadings","residuals"), or c("loadings","intercepts","residuals").}

\item{nsclust}{Vector of length one or two indicating the number(s) of clusters you want to rotate or re-scale for. In case of length two, the vector indicates the minimal and maximal number of clusters. If left unspecified, nsclust is derived from OutputObject.}

\item{design}{For confirmatory factor analysis, matrix (with ncol = nfactors) indicating position of zero loadings with '0' and non-zero loadings with '1'. Leave unspecified for exploratory factor analysis (EFA).}

\item{rescale}{Equal to 1 when solutions need to be rescaled.}

\item{markers}{Matrix (with ncol = nfactors and nrow = number of items) indicating position of marker variable loadings with '1', wheres the rest of the matrix should be equal to '0'. Note that marker variable re-scaling does not make a lot of sense for an unrotated EFA solution.}

\item{rotation}{Rotation criterion to use in case of EFA; currently either "oblimin", "geomin", "varimax" or "target" (i.e., semi-specified oblique Procrustes rotation), whereas 0 = no rotation. (Note: The GPArotation package is loaded or installed for rotation.)}

\item{targetT}{Target matrix to use when rotation = "target". Note that targetT can be one target matrix (with nrow = number of items) or a vertical concatenation of cluster-specific target matrices (with nrow = number of items x number of clusters).
Note that, when nsclust contains >1 numbers of clusters, the same target matrix (i.e., the first one) is used for all clusters.}

\item{targetW}{Weights to be used when rotation = "target". You can set an entire row to zero to make sure that the simple structure is not messed up by a 'complex variable' (i.e., with strong loadings for multiple factors). Set all weights equal to '1' if you prefer fully specified target rotation.
When left unspecified while rotation = target, the default is a targetW where the zeros in the target get a weight of '1' and the non-zeros in the target get a weight of '0'. If this results in all zero weights or too many zeros for the rotation to be identified, a fully specified target rotation (all weights equal to '1') is used instead.}
}
\description{
Allows to rotate and/or re-scale (to marker variable scale) mixture multigroup factor analysis solutions. In this way, you don't need to re-estimate the solutions when you want to change the scale or rotations.
It also counter-rotates the factor (co)variances and factor means.
}
