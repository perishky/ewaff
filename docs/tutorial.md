



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
##         estimate           se          t       p.value   n      p.adjust
## s7  -0.002167527 1.804655e-05 -120.10757  0.000000e+00 500  0.000000e+00
## s50  0.081268506 8.003392e-04  101.54258  0.000000e+00 500  0.000000e+00
## s75 -0.093580406 7.122508e-04 -131.38687  0.000000e+00 500  0.000000e+00
## s33  0.159574185 2.350343e-03   67.89399 2.474930e-252 500 2.474930e-250
## s1   0.314147596 4.889542e-03   64.24888 1.025477e-241 500 1.025477e-239
## s98  0.007285512 1.194706e-04   60.98165 1.634132e-231 499 1.634132e-229
## s25  0.162049835 3.051587e-03   53.10346 4.395821e-206 500 4.395821e-204
## s69 -0.177208720 3.392930e-03  -52.22882 4.605431e-203 500 4.605431e-201
## s17  0.151232349 3.349666e-03   45.14848 2.891114e-177 500 2.891114e-175
## s59  0.085073236 2.045896e-03   41.58239 2.626082e-163 500 2.626082e-161
```

Just for interest, we see if SVA detected batch.

```r
fit <- lm(sites.ret$design[,"sv1"] ~ data$batch)
coef(summary(fit))
```

```
##                Estimate  Std. Error   t value     Pr(>|t|)
## (Intercept) -0.01082079 0.003532752 -3.062990 2.310223e-03
## data$batch1  0.01979530 0.004865142  4.068802 5.495376e-05
## data$batch2  0.01156477 0.004892125  2.363956 1.846466e-02
```
It does seem like it did.

### Generating an EWAS report


```r
sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$pos, methylation,
                         selected.cpg.sites="s58")
```

```
## [ewaff.summary] Thu Jan 27 01:24:53 2022 QQ plots 
## [ewaff.summary] Thu Jan 27 01:24:53 2022 Manhattan plots 
## [ewaff.summary] Thu Jan 27 01:24:53 2022 CpG site plots: 11 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s1 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s7 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s17 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s25 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s33 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s50 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s59 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s69 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s75 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s98 
## [FUN] Thu Jan 27 01:24:53 2022 Plotting s58 
## [ewaff.summary] Thu Jan 27 01:24:53 2022 Sample characteristics 
## [ewaff.sample.characteristics] Thu Jan 27 01:24:53 2022 summarizing variables 
## [summarize.variable] Thu Jan 27 01:24:53 2022 variableB 
## [summarize.variable] Thu Jan 27 01:24:53 2022 continuous 
## [summarize.variable] Thu Jan 27 01:24:53 2022 categorical1 
## [summarize.variable] Thu Jan 27 01:24:53 2022 categorical2 
## [summarize.variable] Thu Jan 27 01:24:53 2022 categorical3 
## [summarize.variable] Thu Jan 27 01:24:53 2022 sv1 
## [ewaff.covariate.associations] Thu Jan 27 01:24:53 2022 covariate associations 
## [FUN] Thu Jan 27 01:24:53 2022 continuous 
## [FUN] Thu Jan 27 01:24:53 2022 categorical1 
## [FUN] Thu Jan 27 01:24:53 2022 categorical2 
## [FUN] Thu Jan 27 01:24:53 2022 categorical3 
## [FUN] Thu Jan 27 01:24:53 2022 sv1
```

```r
ewaff.report(sum.ret, output.file="output/report.html",
             author="Dom Rand",
             study="Associations in my kind of data")
```

```
## [ewaff.report] Thu Jan 27 01:24:53 2022 Writing report as html file to output/report.html
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
##   FALSE    17    0
##   TRUE      8   75
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
## s1 0.1563443 0.9255858          -0.001642581     0.007209047      -0.227850
## s2 1.5711490 0.1954726           0.001777555     0.001739497       1.021879
##    categorical1.p.value categorical2.estimate categorical2.se categorical2.t
## s1            0.8198571           0.002263977     0.007007019      0.3231013
## s2            0.3073384          -0.001775241     0.001690749     -1.0499733
##    categorical2.p.value categorical3.estimate categorical3.se categorical3.t
## s1            0.7467552          0.0023652227     0.007058216     0.33510207
## s2            0.2942437          0.0001363788     0.001703102     0.08007671
##    categorical3.p.value   n p.adjust
## s1            0.7376903 500        1
## s2            0.9362087 500        1
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
##             f       p.value continuous.estimate continuous.se continuous.t
## s1 2069.57711 7.670475e-241       -0.0009363271  0.0024327999   -0.3848763
## s2    5.99221  2.683384e-03        0.0003991080  0.0005883500    0.6783514
## s3   71.33490  6.079467e-28       -0.0008383210  0.0051567687   -0.1625671
## s4  118.62958  8.418467e-43        0.0001080568  0.0009157912    0.1179928
## s5  114.03147  1.917764e-41       -0.0007817620  0.0015434795   -0.5064933
##    continuous.p.value variableB.estimate variableB.se variableB.t
## s1          0.7004947        0.314148004  0.004885031   64.308296
## s2          0.4978666       -0.003990602  0.001181399   -3.377861
## s3          0.8709258        0.123596290  0.010354725   11.936221
## s4          0.9061213       -0.028311180  0.001838897  -15.395741
## s5          0.6127364       -0.046802303  0.003099287  -15.100991
##    variableB.p.value   n      p.adjust
## s1     3.526973e-242 500 7.670475e-239
## s2      7.882165e-04 500  2.683384e-01
## s3      4.941532e-29 500  6.079467e-26
## s4      5.710129e-44 500  8.418467e-41
## s5      1.228769e-42 500  1.917764e-39
```
