library(ewaff)
options(mc.cores=4)

source("simulation-functions.r")

###################################
## construct a random dataset
set.seed(20180220)

n <- 100 ## n samples
s <- 1000  ## s features

## variable of interest and covariates
data <- data.frame(variable=c(rep("A",n/2), rep("B",n/2)), ## variable of interest (two groups)
                   continuous=rnorm(n),                ## continuous covariate
                   categorical=factor(sample(0:3,n,replace=T))) ## categorical covariate

manifest <- data.frame(chr=1, pos=sample(1:(200*s), s))
manifest <- manifest[order(manifest$chr, manifest$pos),]

methylation <- generate.methylation(n, manifest$pos)

## check out methylation correlation structure
r <- sapply(2:nrow(methylation), function(i) cor(methylation[i,], methylation[i-1,]))
d <- tail(manifest$pos,-1) - head(manifest$pos,-1)
d <- 1/exp(d/200)
coef(summary(lm(r~d)))["d",]
##     Estimate    Std. Error       t value      Pr(>|t|) 
## 9.554976e-01  2.020437e-02  4.729164e+01 5.754340e-257 

median(cor(methylation))
## [1] 0.8328628


###################################
## test associations at each CpG site
ret <- ewaff.sites(methylation ~ variable + .,
                 variable.of.interest="variable",
                 methylation=methylation,
                 data=data,
                 generate.confounders="sva",
                 random.subset=0.9,
                 method="glm")

##############################
## test for bumps (there should be none)
bumps.ret <- ewaff.bumps(ret$table$estimate, ret$table$se, ret$table$p.value,
                         methylation,
                         manifest$chr, manifest$pos)
bumps.ret[which(bumps.ret$p.adjust < 0.05),]

##############################
## generate a variable that has a bump in the data
var <- generate.true.bump.var(methylation, manifest$chr, manifest$pos, cluster.sites=20, bump.sites=10, cluster.position=0.5, r=0.6, maxgap=500)
data$bump <- var$var

#############
## test associations at each cpg site
ret <- ewaff.sites(methylation ~ .,
                 variable.of.interest="bump",
                 methylation=methylation,
                 data=data,
                 generate.confounders="sva",
                 random.subset=0.9,
                 method="glm")

ret$table[var$bump.idx,]
##        estimate           se        t      p.value   n     p.adjust
## 230 0.007111196 0.0025068394 2.836718 6.071581e-03 100 1.0000000000
## 231 0.006132350 0.0021617766 2.836718 6.071581e-03 100 1.0000000000
## 232 0.006565303 0.0014801050 4.435701 3.622295e-05 100 0.0362229542
## 233 0.005553250 0.0009832278 5.647980 3.863106e-07 100 0.0003863106
## 234 0.003759809 0.0006413745 5.862112 1.669338e-07 100 0.0001669338
## 235 0.002853838 0.0007081232 4.030143 1.486717e-04 100 0.1486717339
## 236 0.003031969 0.0010774327 2.814068 6.465347e-03 100 1.0000000000
## 237 0.003482685 0.0017010264 2.047402 4.466431e-02 100 1.0000000000
## 238 0.002158379 0.0012959301 1.665506 1.006248e-01 100 1.0000000000
## 239 0.009816485 0.0028214413 3.479245 9.027849e-04 100 0.9027848870

ret$table[ret$table$p.adjust < 0.05,]
##        estimate           se        t      p.value   n     p.adjust
## 232 0.006565303 0.0014801050 4.435701 3.622295e-05 100 0.0362229542
## 233 0.005553250 0.0009832278 5.647980 3.863106e-07 100 0.0003863106
## 234 0.003759809 0.0006413745 5.862112 1.669338e-07 100 0.0001669338

###################
## test for bumps (there should be one)
bumps.ret <- ewaff.bumps(ret$table$estimate, ret$table$se, ret$table$p.value,
                         methylation,
                         manifest$chr, manifest$pos)
bumps.ret[which(bumps.ret$p.adjust < 0.05),]
##   chr start   end n start.idx end.idx start.orig end.orig       p.orig
## 1   1 48753 49292 4       234     237        227      237 1.823888e-08
## 2   1 47978 48723 7       227     233        227      237 1.823888e-08
##        p.value     p.adjust
## 1 7.379003e-10 8.448959e-07
## 2 1.334402e-09 1.527890e-06


