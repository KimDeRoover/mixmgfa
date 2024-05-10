## `mixmgfa` package
R-package for mixture multigroup factor analysis

# Installation
You can download the development version from GitHub as follows:

```javascript
install.packages("devtools"); library(devtools)

devtools::install_github("KimDeRoover/mixmgfa")
```


# How to use `mixmgfa`
For example, to cluster groups on loadings and intercepts, with the number of clusters ranging from 1 to 6:

```javascript
Output<-mixmgfa(data,N_gs,nfactors=1,cluster.spec=c("loadings","intercepts"),nsclust=c(1,6),maxiter=5000,nruns=25,design=design)

```
(where ```design``` should be a matrix with as many rows as there are variables and as many columns as nfactors, containing zeros for zero loadings and ones for nonzero loadings, used for CFA). It is possible to specify only one number of clusters with, for example, ```nsclust = 2```. When performing analyses with many different numbers of clusters, you may encounter longer computation times (especially for larger numbers of clusters and larger data sets). You can use the ```parcomp = 1``` option to mitigate this. The analyses are then performed in parallel instead of sequentially. By default, two cores on your machine are kept free for other tasks, but you can modify this by means of the ```freecores``` option.

**Please note that the package currently cannot deal with missing data.**


For selecting the number of clusters, inspect Output$overview and the plots obtained by:
```javascript
plot(Output$overview)

```
Note that you can choose to plot the AIC instead of the BIC_G by using the ```AIC = 1``` option. A CHull plot cannot be generated when you performed analyses with only one or two different numbers of clusters.


For accessing the corresponding solution, use, for example:
```javascript
Output$MMGFAsolutions$`3.clusters`

```
or
```javascript
Output$MMGFAsolutions[[3]]

```
followed by a '$' sign to access the cluster memberships and parameter sets.


For printing a summary of the selected solution, use, for example:
```javascript
summary(Output$MMGFAsolutions,nclust=3)

```


Note that it can happen that you want to repeat some analyses with more random starts to avoid local maxima (e.g., the ones with larger numbers of clusters), and/or that you want to add analyses with even more clusters. In that case, you can do so by modifying the ```nruns``` and/or ```nsclust``` options. Afterwards, you can vertically concatenate multiple overview tables obtained from the ```mixmgfa``` function and use the ```CHull_mixmgfa``` function to perform the CHull:
```javascript
NewOverviewTable<-CHull_mixmgfa(MergedOverviews)

```
This new overview table only contains the most optimal solution for each number of clusters, is sorted according to the numbers of clusters and includes the (updated) CHull scree ratios. You can get the CHull plot by using this overview table as an input to the ```plot``` function (as indicated above).


# Contribution
This package is a work-in-progress, so please report bugs and check back for updates.
