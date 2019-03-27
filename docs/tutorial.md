



## Running an epigenome-wide association study (EWAS) in `ewaff`

### Preparing R

We first load the `ewaff` library.

```r
library(ewaff)
```

We indicate how many processors are available for performing analyses.
If this is not set explicitly, then the default is to use only
a single processor.

```r
options(mc.cores=4)
```

### Load a dataset

We will in fact generate a random dataset for this demonstration.

```r
set.seed(20171031)
n <- 500 ## n samples
s <- 100  ## s features/cpg sites

## variable of interest and covariates
data <- data.frame(variable=c(rep("A",n/2), rep("B",n/2)),      ## variable of interest (two groups)
                   continuous=rnorm(n),                         ## continuous covariate
                   categorical=factor(sample(0:3,n,replace=T)), ## categorical covariate
                   batch=factor(sample(0:2, size=n, replace=T)))        ## 3 batches

## correlation of each cpg site with the variable of interest
r <- runif(s, min=-1, max=1)

## methylation matrix randomly generated
## with batch effects and associations with the variable
## of interest (rows=cpg sites, cols=samples)
methylation <- t(sapply(r, function(r) {
    v <- as.numeric(data$variable)
    b <- data$batch
    ## batch effect
    b <- rnorm(unique(b), mean=0, sd=0.2)[b]
    ## noise
    e <- rnorm(length(v), mean=0, sd=sqrt(1-r^2))
    ## mean
    m <- runif(1, min=0.2, max=0.8)
    ## signal with mean=0, sd=1
    y <- (r*scale(v) + b + e)
    ## signal with mean=m, sd such that signal is 0..1
    f <- runif(1, min=0, max=min(m,1-m)/max(abs(y)))
    y <- y*f + m
    ## return result
    y
}))
rownames(methylation) <- paste("s", 1:nrow(methylation), sep="")
colnames(methylation) <- paste("p", 1:ncol(methylation), sep="")

## specify cpg site locations
manifest <- data.frame(chr=c(rep(1,s/2), rep(2,s/2)), pos=sample(1:(150*s), s))
manifest <- manifest[order(manifest$chr, manifest$pos),]

## generate 10 likely outliers
outliers <- cbind(sample(1:nrow(methylation), 10, replace=T),
                  sample(1:ncol(methylation), 10, replace=T))
methylation[outliers] <- ifelse(rowMeans(methylation)[outliers[,1]] > 0.5, 0, 1)
```

### Handling outliers

We use the 'iqr' method for handling outliers.5B
The 'iqr' method sets methylation levels that are outside
3*IQR of a CpG site to NA.

```r
methylation <- ewaff.handle.outliers(methylation, method="iqr")[[1]]
```

We can check that most of the outliers added to the dataset in the simulation
were set to NA.

```r
## outliers added
outliers
```

```
##       [,1] [,2]
##  [1,]   83  221
##  [2,]   59  123
##  [3,]   15  405
##  [4,]   76  165
##  [5,]   63  398
##  [6,]   47   86
##  [7,]   14  420
##  [8,]   16  487
##  [9,]   50  281
## [10,]   14  390
```

```r
## outliers identified
which(is.na(methylation), arr.ind=T)
```

```
##     row col
## s47  47  86
## s76  76 165
## s50  50 281
## s14  14 390
## s63  63 398
## s15  15 405
## s14  14 420
## s16  16 487
```
Some 'outliers' were missed because they are not actually outliers.

### Running the EWAS

Now we run the EWAS.  Notice that in this example
we don't explicitly include `batch` as a covariate,
although that is possible.
Here we let surrogate variable analysis (SVA)
generate batch covariates.

```r
sites.ret <- ewaff.sites(methylation ~ variable + continuous + categorical,
                         variable.of.interest="variable",
                         methylation=methylation,
                         data=data,
                         generate.confounders="sva",
                         method="glm")
```

```
## Number of significant surrogate variables is:  1 
## Iteration (out of 5 ):1  2  3  4  5
```

We show the top 10 associations.

```r
top.idx <- order(sites.ret$table$p.value)[1:10]
sites.ret$table[top.idx,]
```

```
##        estimate           se          t       p.value   n      p.adjust
## s58 -0.31296886 0.0023881027 -131.05335  0.000000e+00 500  0.000000e+00
## s61  0.07901268 0.0013680848   57.75422 1.425255e-221 500 1.425255e-219
## s71 -0.26904696 0.0053172012  -50.59936 2.424591e-197 500 2.424591e-195
## s11  0.10080193 0.0022914169   43.99109 8.275566e-173 500 8.275566e-171
## s50  0.02731829 0.0006386444   42.77543 6.928180e-168 499 6.928180e-166
## s57 -0.04791004 0.0011253586  -42.57313 2.977016e-167 500 2.977016e-165
## s78 -0.19526830 0.0048316784  -40.41418 1.378311e-158 500 1.378311e-156
## s49  0.05169961 0.0013563136   38.11774 4.300964e-149 500 4.300964e-147
## s39  0.14158981 0.0038775013   36.51573 2.671361e-142 500 2.671361e-140
## s72  0.29505754 0.0081238012   36.32013 1.843770e-141 500 1.843770e-139
```

Just for interest, we see if SVA detected batch.

```r
fit <- lm(sites.ret$design[,"sv1"] ~ data$batch)
coef(summary(fit))
```

```
##                Estimate  Std. Error   t value     Pr(>|t|)
## (Intercept) -0.01568699 0.003258078 -4.814797 1.958675e-06
## data$batch1  0.03171436 0.004600994  6.892937 1.663669e-11
## data$batch2  0.01519735 0.004755811  3.195533 1.484469e-03
```
It does seem like it did.

### Generating an EWAS report


```r
sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$pos, methylation,
                         selected.cpg.sites="s58")
```

```
## [ewaff.summary] Tue Mar 27 16:17:10 2018 QQ plots 
## [ewaff.summary] Tue Mar 27 16:17:10 2018 Manhattan plots 
## [ewaff.summary] Tue Mar 27 16:17:11 2018 CpG site plots: 10 
## [FUN] Tue Mar 27 16:17:11 2018 Plotting s11
```

```
## Loading required package: betareg
```

```
## [FUN] Tue Mar 27 16:17:12 2018 Plotting s39 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s49 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s50 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s57 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s58 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s61 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s71 
## [FUN] Tue Mar 27 16:17:13 2018 Plotting s72 
## [FUN] Tue Mar 27 16:17:14 2018 Plotting s78 
## [ewaff.summary] Tue Mar 27 16:17:14 2018 Sample characteristics 
## [ewaff.sample.characteristics] Tue Mar 27 16:17:14 2018 summarizing variables 
## [summarize.variable] Tue Mar 27 16:17:14 2018 variableB 
## [summarize.variable] Tue Mar 27 16:17:14 2018 continuous 
## [summarize.variable] Tue Mar 27 16:17:14 2018 categorical1 
## [summarize.variable] Tue Mar 27 16:17:14 2018 categorical2 
## [summarize.variable] Tue Mar 27 16:17:14 2018 categorical3 
## [summarize.variable] Tue Mar 27 16:17:14 2018 sv1 
## [ewaff.covariate.associations] Tue Mar 27 16:17:14 2018 covariate associations 
## [FUN] Tue Mar 27 16:17:14 2018 continuous 
## [FUN] Tue Mar 27 16:17:14 2018 categorical1 
## [FUN] Tue Mar 27 16:17:14 2018 categorical2 
## [FUN] Tue Mar 27 16:17:14 2018 categorical3 
## [FUN] Tue Mar 27 16:17:14 2018 sv1
```

```r
ewaff.report(sum.ret, output.file="output/report.html",
             author="Dom Rand",
             study="Associations in my kind of data")
```

```
## [ewaff.report] Tue Mar 27 16:17:14 2018 Writing report as html file to output/report.html
```

```
## Loading required package: gridExtra
```

### Other kinds of EWAS

#### Methylation is not the outcome
Methylation does not have to be the outcome variable.
In fact, any valid GLM model is possible.
In the following example, we make our binary variable the outcome.

```r
log.ret <- ewaff.sites(variable ~ methylation + continuous + categorical,
                       variable.of.interest="variable",
                       methylation=methylation,
                       data=data,
                       family="binomial",
                       generate.confounders="sva",
                       method="glm")                     
```

```
## Number of significant surrogate variables is:  1 
## Iteration (out of 5 ):1  2  3  4  5
```

Associations are identical to
those identified when methylation was the outcome.

```r
table(sites.ret$table$p.adjust < 0.05, log.ret$table$p.adjust < 0.05)
```

```
##        
##         FALSE TRUE
##   FALSE    21    0
##   TRUE      3   76
```


#### Variable of interest is complex

The variable of interest may be categorical with more than two categories.

```r
cats.ret <- ewaff.sites(methylation ~ categorical + variable + continuous,
                        variable.of.interest="categorical",
                        methylation=methylation,
                        data=data,
                        generate.confounders="sva",
                        method="limma")
```

```
## Number of significant surrogate variables is:  1 
## Iteration (out of 5 ):1  2  3  4  5
```

In this case, an f-statistic and p-value is calculated for the variable
along with statistics each binary 'dummy' variable.

```r
cats.ret$table[1:2,]
```

```
##           f   p.value categorical1.estimate categorical1.se categorical1.t
## s1 1.688244 0.1686009         -0.0008108143      0.01315824    -0.06162026
## s2 1.582263 0.1927593         -0.0099781474      0.01116591    -0.89362607
##    categorical1.p.value categorical2.estimate categorical2.se
## s1            0.9508902           0.009613312      0.01381623
## s2            0.3719571           0.004912092      0.01172427
##    categorical2.t categorical2.p.value categorical3.estimate
## s1      0.6957987            0.4868821            0.02637100
## s2      0.4189680            0.6754217           -0.01875292
##    categorical3.se categorical3.t categorical3.p.value   n p.adjust
## s1      0.01379154       1.912115           0.05643962 500        1
## s2      0.01170331      -1.602360           0.10971522 500        1
```

The variable of interest may actuually be multiple variables.

```r
vars.ret <- ewaff.sites(methylation ~ categorical + variable + continuous,
                        variable.of.interest=c("continuous","variable"),
                        methylation=methylation,
                        data=data,
                        generate.confounders="sva",
                        method="limma")
```

```
## Number of significant surrogate variables is:  1 
## Iteration (out of 5 ):1  2  3  4  5
```

Here again, an f-statistic is calcualted along with individual statistics
for each variable.

```r
vars.ret$table[1:5,]
```

```
##             f      p.value continuous.estimate continuous.se continuous.t
## s1  0.5152476 5.976733e-01        -0.003676184   0.004756168   -0.7729298
## s2 40.7276598 4.243540e-17        -0.004569485   0.004026465   -1.1348628
## s3 22.3026060 5.330194e-10         0.002470783   0.001706407    1.4479450
## s4 31.7877296 1.031762e-13        -0.005395120   0.005480739   -0.9843783
## s5 27.2812407 5.777345e-12         0.001329901   0.001180262    1.1267844
##    continuous.p.value variableB.estimate variableB.se variableB.t
## s1          0.4399336        0.006244405  0.009694717   0.6441039
## s2          0.2569832       -0.073640700  0.008207329  -8.9725538
## s3          0.1482668       -0.022583659  0.003478248  -6.4928256
## s4          0.3254114       -0.088577937  0.011171643  -7.9288194
## s5          0.2603807        0.017608508  0.002405782   7.3192443
##    variableB.p.value   n     p.adjust
## s1      5.198070e-01 500 1.000000e+00
## s2      6.040519e-18 500 4.243540e-15
## s3      2.055153e-10 500 5.330194e-08
## s4      1.486003e-14 500 1.031762e-11
## s5      1.017656e-12 500 5.777345e-10
```
