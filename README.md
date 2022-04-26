## `mixmgfa` package
R-package for mixture multigroup factor analysis

# Installation
You can download the development version from GitHub as follows:

```javascript
install.packages("devtools"); library(devtools)

devtools::install_github("KimDeRoover/mixmgfa")
```


# How to use `mixmgfa`
For example, to cluster groups on loadings and intercepts:

```javascript
Output<-mixmgfa(data,N_gs,nfactors=1,cluster.spec=c("loadings","intercepts"),nsclust=c(1,6),maxiter=5000,nruns=25,design=design)

```
(where design = a matrix with as many rows as there are variables and as many columns as nfactors, containing zeros for zero loadings and ones for nonzero loadings.)


For selecting the number of clusters, inspect Output$overview and the plots obtained by:
```javascript
plot(Output$overview)

```

For accessing the corresponding solution, use, for example:
```javascript
Output$MMGFAsolutions$`3.clusters`

```
followed by a '$' sign to access the cluster memberships and parameter sets.

For printing a summary of the selected solution, use, for example:
```javascript
summary(Output$MMGFAsolutions,nclust=3)

```

# Contribution
This package is a work-in-progress, so please report bugs and check back for updates (including more rotation and scaling options and improvements of the multistart procedure).
