



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
    v <- sign(data$variable=="A")
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
##  [1,]   74  173
##  [2,]   98  452
##  [3,]   53  425
##  [4,]    8  409
##  [5,]   62   96
##  [6,]   82  149
##  [7,]   81  121
##  [8,]   38  151
##  [9,]   48  148
## [10,]   95  135
```

```r
## outliers identified
which(is.na(methylation), arr.ind=T)
```

```
##     row col
## s81  81 121
## s95  95 135
## s48  48 148
## s82  82 149
## s38  38 151
## s74  74 173
## s8    8 409
## s98  98 452
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
                         random.subset=1,
                         n.confounders=1,
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
## s7   0.00213271 1.805607e-05  118.11594  0.000000e+00 500  0.000000e+00
## s50 -0.07868680 7.695914e-04 -102.24490  0.000000e+00 500  0.000000e+00
## s75  0.09976573 7.585788e-04  131.51662  0.000000e+00 500  0.000000e+00
## s33 -0.16459231 2.516830e-03  -65.39668 4.140361e-245 500 4.140361e-243
## s1  -0.31685327 4.872176e-03  -65.03322 4.859680e-244 500 4.859680e-242
## s98 -0.00739083 1.193550e-04  -61.92311 2.071756e-234 499 2.071756e-232
## s25 -0.16856506 3.189139e-03  -52.85598 3.119356e-205 500 3.119356e-203
## s69  0.17176635 3.529924e-03   48.66007 2.283911e-190 500 2.283911e-188
## s27  0.06192653 1.437160e-03   43.08952 2.742106e-169 500 2.742106e-167
## s59 -0.08562959 2.003567e-03  -42.73858 6.607039e-168 500 6.607039e-166
```

Just for interest, we see if SVA detected batch.

```r
fit <- lm(sites.ret$design[,"sv1"] ~ data$batch)
coef(summary(fit))
```

```
##                Estimate  Std. Error   t value     Pr(>|t|)
## (Intercept) -0.01096824 0.003530798 -3.106448 2.001749e-03
## data$batch1  0.02012731 0.004862450  4.139334 4.091502e-05
## data$batch2  0.01165863 0.004889418  2.384462 1.747823e-02
```
It does seem like it did.

### Generating an EWAS report


```r
sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$pos, methylation,
                         selected.cpg.sites="s58")
```

```
## [ewaff.summary] Tue Jun 25 12:28:24 2024 QQ plots 
## [ewaff.summary] Tue Jun 25 12:28:24 2024 Manhattan plots 
## [ewaff.summary] Tue Jun 25 12:28:24 2024 CpG site plots: 11 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s1 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s7 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s25 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s27 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s33 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s50 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s59 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s69 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s75 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s98 
## [FUN] Tue Jun 25 12:28:24 2024 Plotting s58 
## [ewaff.summary] Tue Jun 25 12:28:24 2024 Sample characteristics 
## [ewaff.sample.characteristics] Tue Jun 25 12:28:24 2024 summarizing variables 
## [summarize.variable] Tue Jun 25 12:28:24 2024 variableB 
## [summarize.variable] Tue Jun 25 12:28:24 2024 continuous 
## [summarize.variable] Tue Jun 25 12:28:24 2024 categorical1 
## [summarize.variable] Tue Jun 25 12:28:24 2024 categorical2 
## [summarize.variable] Tue Jun 25 12:28:24 2024 categorical3 
## [summarize.variable] Tue Jun 25 12:28:24 2024 sv1 
## [ewaff.covariate.associations] Tue Jun 25 12:28:24 2024 covariate associations 
## [FUN] Tue Jun 25 12:28:24 2024 continuous 
## [FUN] Tue Jun 25 12:28:24 2024 categorical1 
## [FUN] Tue Jun 25 12:28:24 2024 categorical2 
## [FUN] Tue Jun 25 12:28:24 2024 categorical3 
## [FUN] Tue Jun 25 12:28:24 2024 sv1
```

```r
ewaff.report(sum.ret, output.file="output/report.html",
             author="Dom Rand",
             study="Associations in my kind of data")
```

```
## [ewaff.report] Tue Jun 25 12:28:24 2024 Writing report as html file to output/report.html
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
data$variable01 = sign(data$variable=="A")
log.ret <- ewaff.sites(variable01 ~ methylation + continuous + categorical,
                       variable.of.interest="variable01",
                       methylation=methylation,
                       data=data,
                       family="binomial",
                       generate.confounders="sva",
                       random.subset=1,
                       n.confounders=1,
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
##   FALSE    16    0
##   TRUE      7   77
```


#### Variable of interest is complex

The variable of interest may be categorical with more than two categories.

```r
cats.ret <- ewaff.sites(methylation ~ categorical + variable + continuous,
                        variable.of.interest="categorical",
                        methylation=methylation,
                        data=data,
                        generate.confounders="sva",
                        random.subset=1,
                        n.confounders=1,
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
##            f   p.value categorical1.estimate categorical1.se categorical1.t
## s1 0.1767259 0.9121571          -0.001994886     0.007176413     -0.2779781
## s2 1.5418963 0.2027868           0.001482874     0.001636584      0.9060787
##    categorical1.p.value categorical2.estimate categorical2.se categorical2.t
## s1            0.7811455           0.001984101     0.006987262      0.2839597
## s2            0.3653359          -0.001797667     0.001593449     -1.1281615
##    categorical2.p.value categorical3.estimate categorical3.se categorical3.t
## s1            0.7765602          0.0024931244     0.007042619      0.3540053
## s2            0.2597993          0.0002118757     0.001606073      0.1319216
##    categorical3.p.value   n p.adjust
## s1            0.7234859 500        1
## s2            0.8951000 500        1
```

The variable of interest may actuually be multiple variables.

```r
vars.ret <- ewaff.sites(methylation ~ categorical + variable + continuous,
                        variable.of.interest=c("continuous","variable"),
                        methylation=methylation,
                        data=data,
                        generate.confounders="sva",
                        random.subset=1,
                        n.confounders=1,
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
##              f       p.value continuous.estimate continuous.se continuous.t
## s1 2119.232961 4.050498e-243       -9.330774e-04  0.0024238458   -0.3849574
## s2    1.240829  2.900440e-01        3.746043e-04  0.0005523450    0.6782072
## s3   88.799345  1.133802e-33       -6.556652e-04  0.0040375141   -0.1623933
## s4  152.001930  3.594024e-52        8.924994e-05  0.0007562391    0.1180181
## s5  121.030019  1.671741e-43       -7.244696e-04  0.0014309037   -0.5063021
##    continuous.p.value variableB.estimate variableB.se variableB.t
## s1          0.7004347        -0.31685291  0.004867744  -65.092359
## s2          0.4979579         0.00159511  0.001109259    1.437996
## s3          0.8710626        -0.10804985  0.008108430  -13.325619
## s4          0.9061013         0.02647609  0.001518734   17.432992
## s5          0.6128704         0.04463645  0.002873645   15.533044
##    variableB.p.value   n      p.adjust
## s1     1.656168e-244 500 4.050498e-241
## s2      1.510684e-01 500  1.000000e+00
## s3      7.944998e-35 500  1.133802e-31
## s4      2.144654e-53 500  3.594024e-50
## s5      1.357139e-44 500  1.671741e-41
```
