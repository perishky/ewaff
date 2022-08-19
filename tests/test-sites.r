library(ewaff)
options(mc.cores=4)

source("simulation-functions.r")

###################################
## construct a random dataset
set.seed(20171031)

n <- 500 ## n samples
s <- 10   ## s features

## variable of interest and covariates
data <- data.frame(variable=c(rep("A",n/2), rep("B",n/2)), ## variable of interest (two groups)
                   continuous=rnorm(n),                ## continuous covariate
                   categorical=factor(sample(0:3,n,replace=T))) ## categorical covariate

## correlation of each cpg site with the variable of interest
r <- runif(s, min=-1, max=1) 

## methylation matrix (rows=cpg sites, cols=samples)
methylation <- t(sapply(r, function(r) rcor(as.numeric(as.factor(data$variable)), r))) 

## batch variable, 3 batches
batch <- sample(0:2,size=nrow(data),replace=T)

## add a systematic effect for each batch
for (b in unique(batch)) {
    idx <- which(batch == b)
    methylation[,idx] <- methylation[,idx] + rnorm(nrow(methylation))
}

###################################
## perform EWAS
ret <- ewaff.sites(methylation ~ variable + .,
                   variable.of.interest="variable",
                   methylation=methylation,
                   data=data,
                   generate.confounders="sva",
                   random.subset=0.9,
                   method="glm")

## list summary statistics
ret$table
##        estimate         se            t       p.value   n      p.adjust
## 1   0.316122665 0.04848556   6.51993471  1.743581e-10 500  1.743581e-09
## 2   0.943354480 0.04403419  21.42322717  1.864826e-72 500  1.864826e-71
## 3  -0.003178783 0.06596018  -0.04819246  9.615824e-01 500  1.000000e+00
## 4  -0.570421539 0.04596814 -12.40906262  5.829200e-31 500  5.829200e-30
## 5  -0.355428046 0.04985035  -7.12990075  3.596122e-12 500  3.596122e-11
## 6  -0.426801680 0.03842641 -11.10698743  9.836881e-26 500  9.836881e-25
## 7   0.027996441 0.01710605   1.63664008  1.023439e-01 500  1.000000e+00
## 8  -0.966882957 0.01283460 -75.33408678 1.217838e-272 500 1.217838e-271
## 9  -0.085198308 0.04106466  -2.07473541  3.852928e-02 500  3.852928e-01
## 10  0.005233260 0.04818505   0.10860753  9.135580e-01 500  1.000000e+00

## check association statistics of the first cpg site manually
coef(summary(glm(methylation[1,] ~ ., data=as.data.frame(ret$design))))["variableB",]
##     Estimate   Std. Error      t value     Pr(>|t|) 
## 3.161227e-01 4.848556e-02 6.519935e+00 1.743581e-10 


ret2 <- ewaff.sites(methylation ~ variable + .,
                    variable.of.interest="variable",
                    methylation=methylation,
                    data=data,
                    generate.confounders="sva",
                    random.subset=0.9,
                 method="limma")
ret2$table
##         estimate         se             t       p.value   n      p.adjust
## 1   3.025911e-01 0.04804993  6.297429e+00  6.670168e-10 500  6.670168e-09
## 2   9.253519e-01 0.04132022  2.239465e+01  2.716168e-77 500  2.716168e-76
## 3  -7.479063e-03 0.06619187 -1.129907e-01  9.100836e-01 500  1.000000e+00
## 4  -5.655820e-01 0.04595159 -1.230821e+01  1.436577e-30 500  1.436577e-29
## 5  -3.581953e-01 0.04961421 -7.219612e+00  1.969358e-12 500  1.969358e-11
## 6  -4.275574e-01 0.03836230 -1.114525e+01  6.708131e-26 500  6.708131e-25
## 7   5.189992e-02 0.02759045  1.881083e+00  6.054579e-02 500  6.054579e-01
## 8  -9.672749e-01 0.01301092 -7.434333e+01 3.440480e-271 500 3.440480e-270
## 9  -6.922289e-02 0.03324151 -2.082423e+00  3.781534e-02 500  3.781534e-01
## 10 -3.605528e-05 0.04524844 -7.968291e-04  9.993645e-01 500  1.000000e+00



ret3 <- ewaff.sites(methylation ~ variable + .,
                    variable.of.interest="variable",
                    methylation=methylation,
                    data=data,
                    generate.confounders="sva",
                    random.subset=0.9,
                    method="rlm")

ret3$table
##       estimate         se           z      p.value   n     p.adjust
## 1   0.23484112 0.05444909   4.3130404 1.610247e-05 500 1.610247e-04
## 2   0.88196621 0.05409215  16.3048837 9.111219e-60 500 9.111219e-59
## 3  -0.01903942 0.06900372  -0.2759188 7.826104e-01 500 1.000000e+00
## 4  -0.57120497 0.04837506 -11.8078390 3.555833e-32 500 3.555833e-31
## 5  -0.36746076 0.05227654  -7.0291715 2.077634e-12 500 2.077634e-11
## 6  -0.42516487 0.03961289 -10.7329942 7.124941e-27 500 7.124941e-26
## 7   0.12457382 0.06722023   1.8532191 6.385096e-02 500 6.385096e-01
## 8  -0.96682115 0.01350887 -71.5693829 0.000000e+00 500 0.000000e+00
## 9  -0.01569902 0.02092514  -0.7502470 4.531059e-01 500 1.000000e+00
## 10 -0.02050618 0.04484952  -0.4572217 6.475117e-01 500 1.000000e+00

###############################################
## perform another EWAS with the variable of interest as dependent variable
## so use logistic regression
ret4 <- ewaff.sites(variable ~ methylation + .,
                    variable.of.interest="variable",
                    methylation=methylation,
                    data=data,
                    family="binomial",
                    generate.confounders="sva",
                    random.subset=0.9,
                    method="glm")

## list summary statistics
ret4$table
##         estimate           se             t      p.value   n     p.adjust
## 1   1.076818e+00 1.834087e-01  5.871139e+00 4.328116e-09 500 4.328116e-08
## 2   4.540911e+00 4.104924e-01  1.106211e+01 1.915441e-28 500 1.915441e-27
## 3  -1.409149e-02 1.240645e-01 -1.135819e-01 9.095692e-01 500 1.000000e+00
## 4  -2.134386e+00 2.207007e-01 -9.670950e+00 4.006423e-22 500 4.006423e-21
## 5  -1.203186e+00 1.825627e-01 -6.590537e+00 4.382391e-11 500 4.382391e-10
## 6  -2.397286e+00 2.655627e-01 -9.027196e+00 1.761260e-19 500 1.761260e-18
## 7   5.662776e-01 3.004990e-01  1.884457e+00 5.950313e-02 500 5.950313e-01
## 8  -1.067435e+02 3.935895e+04 -2.712051e-03 9.978361e-01 500 1.000000e+00
## 9  -5.201814e-01 2.500793e-01 -2.080066e+00 3.751951e-02 500 3.751951e-01
## 10 -8.269307e-06 1.815050e-01 -4.555965e-05 9.999636e-01 500 1.000000e+00

## manually verify statistics for the first cpg site 
coef(summary(glm(variable ~ methylation[1,] + ., data=as.data.frame(ret4$design), family="binomial")))["methylation[1, ]",]
##     Estimate   Std. Error      z value     Pr(>|z|) 
## 1.076818e+00 1.834087e-01 5.871139e+00 4.328116e-09 



############################
## dataset with multiple variables of interest
## here we are asking if adding two or more variables
## to the null model (without them) significantly
## improves the fit of the model.

data$variable1 <- data$variable
data$variable2 <- sample(c("X","Y","Z"), nrow(data), replace=T)
data$variable <- NULL

ret5 <- ewaff.sites(methylation ~ .,
                    variable.of.interest=c("variable1","variable2"),
                    methylation=methylation,
                    data=data,
                    family=gaussian, 
                    generate.confounders="sva",
                    random.subset=0.9,
                    method="glm")

ret5$table
##               f       p.value variable1B.estimate variable1B.se variable1B.t
## 1     3.6672782  1.232638e-02         0.264651694    0.09424123   2.80823675
## 2    20.8025047  1.053210e-12         0.859356867    0.11747963   7.31494381
## 3     0.2786005  8.408527e-01        -0.010650354    0.06999358  -0.15216187
## 4    36.3309211  3.186729e-21        -0.551625446    0.05412765 -10.19119482
## 5    16.9086941  1.814526e-10        -0.367657487    0.05290070  -6.94995481
## 6    44.4668425  1.926539e-25        -0.419124244    0.03861021 -10.85526886
## 7     3.0396792  2.870219e-02         0.134786950    0.13792171   0.97727144
## 8  1810.5543660 2.488993e-265        -0.969651850    0.01320512 -73.43000217
## 9     2.6889428  4.581528e-02        -0.006664799    0.10896246  -0.06116601
## 10    1.6451558  1.780631e-01        -0.029797742    0.06202383  -0.48042410
##    variable1B.p.value variable2Y.estimate variable2Y.se variable2Y.t
## 1        5.179367e-03          0.14986113    0.11449236    1.3089181
## 2        1.053045e-12          0.33742627    0.14272436    2.3641813
## 3        8.791216e-01         -0.03925514    0.08503423   -0.4616393
## 4        2.943524e-22         -0.08650937    0.06575893   -1.3155532
## 5        1.165145e-11          0.07713968    0.06426833    1.2002753
## 6        9.331653e-25         -0.13148624    0.04690701   -2.8031255
## 7        3.289151e-01         -0.47801472    0.16755917   -2.8528115
## 8       2.878655e-267          0.01627414    0.01604271    1.0144258
## 9        9.512518e-01         -0.36194879    0.13237697   -2.7342278
## 10       6.311394e-01          0.16071595    0.07535189    2.1328723
##    variable2Y.p.value variable2Z.estimate variable2Z.se variable2Z.t
## 1         0.191173107         0.170232804    0.11325665    1.5030712
## 2         0.018457443         0.186469088    0.14118394    1.3207528
## 3         0.644544021         0.037363305    0.08411646    0.4441854
## 4         0.188936752        -0.109278939    0.06504919   -1.6799430
## 5         0.230610033         0.096108307    0.06357468    1.5117388
## 6         0.005260834        -0.030834002    0.04640074   -0.6645153
## 7         0.004515921        -0.314554593    0.16575070   -1.8977572
## 8         0.310878215         0.007046226    0.01586957    0.4440088
## 9         0.006478146        -0.254062915    0.13094823   -1.9401784
## 10        0.033429526         0.039183771    0.07453861    0.5256842
##    variable2Z.p.value   n      p.adjust
## 1          0.13346226 500  1.232638e-01
## 2          0.18719778 500  1.053210e-11
## 3          0.65710382 500  1.000000e+00
## 4          0.09360320 500  3.186729e-20
## 5          0.13124234 500  1.814526e-09
## 6          0.50667191 500  1.926539e-24
## 7          0.05831316 500  2.870219e-01
## 8          0.65723145 500 2.488993e-264
## 9          0.05292903 500  4.581528e-01
## 10         0.59934449 500  1.000000e+00

ret6 <- ewaff.sites(methylation ~ .,
                    variable.of.interest=c("variable1","variable2"),
                    methylation=methylation,
                    data=data,
                    family=gaussian, 
                    generate.confounders="sva",
                    random.subset=0.9,
                    method="limma")

ret6$table
##               f       p.value variable1B.estimate variable1B.se variable1B.t
## 1     3.6776608  1.215208e-02         0.264651694    0.09410811   2.81220919
## 2    20.8680671  9.606266e-13         0.859356867    0.11729494   7.32646185
## 3     0.2791858  8.404314e-01        -0.010650354    0.06992017  -0.15232163
## 4    36.3675636  2.996835e-21        -0.551625446    0.05410038 -10.19633283
## 5    16.9235941  1.771534e-10        -0.367657487    0.05287741  -6.95301629
## 6    44.3954622  2.045137e-25        -0.419124244    0.03864124 -10.84655266
## 7     3.0497422  2.831450e-02         0.134786950    0.13769398   0.97888775
## 8  1737.9018862 6.153800e-262        -0.969651850    0.01347831 -71.94164819
## 9     2.6971647  4.531366e-02        -0.006664799    0.10879625  -0.06125945
## 10    1.6478805  1.774451e-01        -0.029797742    0.06197253  -0.48082178
##    variable1B.p.value variable2Y.estimate variable2Y.se variable2Y.t
## 1        5.116105e-03          0.14986113    0.11433063    1.3107697
## 2        9.698286e-13          0.33742627    0.14249998    2.3679039
## 3        8.789955e-01         -0.03925514    0.08494504   -0.4621240
## 4        2.772782e-22         -0.08650937    0.06572579   -1.3162165
## 5        1.137717e-11          0.07713968    0.06424003    1.2008041
## 6        9.878939e-25         -0.13148624    0.04694470   -2.8008747
## 7        3.281147e-01         -0.47801472    0.16728250   -2.8575298
## 8       7.116790e-264          0.01627414    0.01637461    0.9938644
## 9        9.511774e-01         -0.36194879    0.13217505   -2.7384048
## 10       6.308561e-01          0.16071595    0.07528956    2.1346378
##    variable2Y.p.value variable2Z.estimate variable2Z.se variable2Z.t
## 1         0.190544812         0.170232804    0.11309666    1.5051974
## 2         0.018273467         0.186469088    0.14096198    1.3228325
## 3         0.644195893         0.037363305    0.08402823    0.4446518
## 4         0.188711994        -0.109278939    0.06501642   -1.6807900
## 5         0.230402865         0.096108307    0.06354669    1.5124047
## 6         0.005296323        -0.030834002    0.04643803   -0.6639817
## 7         0.004449740        -0.314554593    0.16547702   -1.9008959
## 8         0.320775616         0.007046226    0.01619788    0.4350092
## 9         0.006396830        -0.254062915    0.13074849   -1.9431423
## 10        0.033282314         0.039183771    0.07447696    0.5261193
##    variable2Z.p.value   n      p.adjust
## 1          0.13291264 500  1.215208e-01
## 2          0.18650331 500  9.606266e-12
## 3          0.65676620 500  1.000000e+00
## 4          0.09343605 500  2.996835e-20
## 5          0.13107060 500  1.771534e-09
## 6          0.50701191 500  2.045137e-24
## 7          0.05789757 500  2.831450e-01
## 8          0.66374582 500 6.153800e-261
## 9          0.05256683 500  4.531366e-01
## 10         0.59904147 500  1.000000e+00


cbind(ret5$table$variable1B.t, ret6$table$variable1B.t)
##               [,1]         [,2]
##  [1,]   2.80823675   2.81220919
##  [2,]   7.31494381   7.32646185
##  [3,]  -0.15216187  -0.15232163
##  [4,] -10.19119482 -10.19633283
##  [5,]  -6.94995481  -6.95301629
##  [6,] -10.85526886 -10.84655266
##  [7,]   0.97727144   0.97888775
##  [8,] -73.43000217 -71.94164819
##  [9,]  -0.06116601  -0.06125945
## [10,]  -0.48042410  -0.48082178

cbind(ret5$table$f, ret6$table$f)
##               [,1]         [,2]
##  [1,]    3.6672782    3.6776608
##  [2,]   20.8025047   20.8680671
##  [3,]    0.2786005    0.2791858
##  [4,]   36.3309211   36.3675636
##  [5,]   16.9086941   16.9235941
##  [6,]   44.4668425   44.3954622
##  [7,]    3.0396792    3.0497422
##  [8,] 1810.5543660 1737.9018862
##  [9,]    2.6889428    2.6971647
## [10,]    1.6451558    1.6478805


ret.int1 <- ewaff.sites(
    methylation ~ variable1 * continuous + categorical ,
    variable.of.interest="variable1:continuous",
    methylation=methylation,
    data=data,
    generate.confounders=NULL,
    random.subset=0.9,
    method="glm")

ret.int2 <- ewaff.sites(
    methylation ~ variable1 * continuous + categorical ,
    variable.of.interest="variable1*continuous",
    methylation=methylation,
    data=data,
    generate.confounders=NULL,
    random.subset=0.9,
    method="glm")

all(abs(ret.int1$table$t - ret.int2$table$t) < 2e-16)
## [1] TRUE

fit <- lm(methylation[4,] ~ variable1 * continuous + categorical, data=data)
abs(coef(summary(fit))["variable1B:continuous","t value"] - ret.int1$table$t[4]) < 2e-16
## [1] TRUE
