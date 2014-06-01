Quick Comparison of RUV and SVA for diagnostic purposes
========================================================


## Summary 

RUV and SVA are two methods for removing unknown/unwanted batch effects from RNA-seq data. In this analysis we show that: (1) SVA produces very similar estimates of unknown batch effects to both RUV with control probes and RUV with empirical control probes, (2) SVA produces very similar adjusted counts to the RUV approaches, and (3) SVA produces very similar DE results to the RUV approaches. This first analysis is based solely on the example data set in the RUVSeq vignette. 


### Install packages

This analysis depends on the current devel version of Bioconductor because the RUV package is only availabe in devel. See session information at the end of the comparison. This chunk of code is not set to be evaluated. If you want to rerun the code and ensure packages are installed, set `eval=TRUE` in the source code. 



```r
source("http://bioconductor.org/biocLite.R")
BiocInstaller::useDevel()
biocLite("sva")
biocLite("RUVSeq")
biocLite("zebrafishRNASeq",type="source")
install.packages("devtools")
library(devtools)
devtools::install_github('RSkittleBrewer', 'alyssafrazee')
install.packages("doRNG",type="source")
biocLite("ffpe")
```

Now load the libraries


```r
library(sva)
```

```
## Loading required package: corpcor
## Loading required package: mgcv
## Loading required package: nlme
## This is mgcv 1.7-29. For overview type 'help("mgcv-package")'.
```

```r
library(RUVSeq)
```

```
## Loading required package: EDASeq
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: ShortRead
## Loading required package: BiocParallel
## Loading required package: Biostrings
## Loading required package: S4Vectors
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## 
## The following object is masked from 'package:nlme':
## 
##     collapse
## 
## Loading required package: XVector
## Loading required package: Rsamtools
## Loading required package: GenomeInfoDb
## Loading required package: GenomicRanges
## Loading required package: GenomicAlignments
## Loading required package: edgeR
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
library(zebrafishRNASeq)
library(RSkittleBrewer)
library(ffpe)
```

```
## Loading required package: TTR
## Loading required package: xts
## Loading required package: zoo
## 
## Attaching package: 'zoo'
## 
## The following object is masked from 'package:Rsamtools':
## 
##     index
## 
## The following objects are masked from 'package:base':
## 
##     as.Date, as.Date.numeric
## 
## 
## Attaching package: 'xts'
## 
## The following objects are masked from 'package:GenomicAlignments':
## 
##     first, last
```

```
## Warning: replacing previous import by 'graphics::image' when loading
## 'methylumi'
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
```

### Load and process the zebrafish data

The first comparison will be on the zebrafish data like those in the RUVSeq vignette. For the vignette code see this link http://bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.R.


```r
data(zfGenes)
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
```


### Calculate RUV adjustment factors


```r
## RUV using the known spikeins
ruvCp <- RUVg(counts(set), spikes, k=1)

## RUV using control probes
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
ruvEmp <- RUVg(counts(set), empirical, k=1)
```

### Calculate sva adjustment factors


```r
dat = log(as.matrix(filtered)+1)
mod = model.matrix(~rep(c(0,1),each=3))
mod0 = cbind(mod[,1])
sv1 = sva(dat,mod,mod0,n.sv=1)
```

```
## Number of significant surrogate variables is:  1 
## Iteration (out of 5 ):1  2  3  4  5
```


### Compare sva and ruv adjustment factors


```r
# Get colors
trop = RSkittleBrewer('tropical')

par(mfrow=c(1,2))
plot(sv1$sv,ruvCp$W,col=trop[1],
     pch=19,cex=2,xlab="SVA",ylab="RUV Control Probes")
abline(c(0,1))

plot(sv1$sv,ruvEmp$W,col=trop[2],
     pch=19,cex=2,xlab="SVA",ylab="RUV Empirical Control")
abline(c(0,1))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 

```r
cor(sv1$sv,ruvCp$W)
```

```
##        W_1
## [1,] 0.974
```

```r
cor(sv1$sv,ruvEmp$W)
```

```
##        W_1
## [1,] 0.996
```


### Compare adjusted values


```r
n = dim(filtered)[2]
Id  = diag(n)
modsv = cbind(sv1$sv)
resid = dat %*% (Id - modsv %*% solve(t(modsv) %*% modsv) %*% t(modsv))
svaCounts = round(exp(resid))

corSvaCp = corSvaEmp = rep(NA,dim(dat)[1])
for(i in 1:dim(dat)[1]){corSvaCp[i] = cor(svaCounts[i,],ruvCp$normalizedCounts[i,])}
for(i in 1:dim(dat)[1]){corSvaEmp[i] = cor(svaCounts[i,],ruvEmp$normalizedCounts[i,])}

par(mfrow=c(1,2))
hist(corSvaCp,ylab="Gene Specific Correlation",main="SVA vs. RUV CP",col=trop[1],breaks=100)
hist(corSvaCp,ylab="Gene Specific Correlation",main="SVA vs. RUV CP",col=trop[2],breaks=100)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 

```r
quantile(corSvaCp,probs=seq(0,1,length=20))
```

```
##      0%  5.263%  10.53%  15.79%  21.05%  26.32%  31.58%  36.84%  42.11% 
## -0.6009  0.8287  0.8938  0.9249  0.9423  0.9548  0.9637  0.9704  0.9758 
##  47.37%  52.63%  57.89%  63.16%  68.42%  73.68%  78.95%  84.21%  89.47% 
##  0.9803  0.9843  0.9877  0.9905  0.9930  0.9949  0.9966  0.9980  0.9990 
##  94.74%    100% 
##  0.9997  1.0000
```

```r
quantile(corSvaEmp,probs=seq(0,1,length=20))
```

```
##       0%   5.263%   10.53%   15.79%   21.05%   26.32%   31.58%   36.84% 
## -0.01717  0.96958  0.98170  0.98727  0.99051  0.99259  0.99404  0.99514 
##   42.11%   47.37%   52.63%   57.89%   63.16%   68.42%   73.68%   78.95% 
##  0.99601  0.99673  0.99731  0.99781  0.99827  0.99866  0.99901  0.99931 
##   84.21%   89.47%   94.74%     100% 
##  0.99958  0.99978  0.99993  1.00000
```

### Calculate DE results for the three approaches


```r
###
## Calculate DE results for Control Probes RUV
###

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=ruvCp$normalizedCounts, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

cpResults <- lrt$table

###
## Calculate DE results for Empirical Control Probes RUV
###

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=ruvEmp$normalizedCounts, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

empResults <- lrt$table


###
## Calculate DE results for sva
###

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=svaCounts, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

svaResults = lrt$table


###
## Calculate DE results for no normalization
###

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=filtered, group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

noneResults = lrt$table
```

### Compare results for the three approaches

Make concordance at the top plots to compare DE results from the three approaches. 



```r
empLr = empResults$LR
names(empLr) = rownames(empResults)
cpLr = cpResults$LR
names(cpLr) = rownames(cpResults)
svaLr = svaResults$LR
names(svaLr) = rownames(svaResults)
noneLr = noneResults$LR
names(noneLr) = rownames(noneResults)

none_sva = CATplot(-noneLr,-svaLr,maxrank=1000,make.plot=F)
ruv_cp_sva = CATplot(-cpLr,-svaLr,maxrank=1000,make.plot=F)
ruv_emp_sva = CATplot(-empLr,-svaLr,maxrank=1000,make.plot=F)

plot(none_sva,ylim=c(0,1),col=trop[1],lwd=2,type="l")
lines(ruv_cp_sva,ylim=c(0,1),col=trop[2],lwd=2)
lines(ruv_emp_sva,ylim=c(0,1),col=trop[3],lwd=2)
legend(600,0.2,legend=c("None vs. SVA","RUV CP vs SVA","RUV Emp vs SVA"),col=trop[1:3],lwd=2)
```

![plot of chunk resultscomp](figure/resultscomp.png) 


