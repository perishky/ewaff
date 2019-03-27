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

median(cor(methylation))


###################################
## test associations at each CpG site
ret <- ewaff.sites(methylation ~ variable + .,
                 variable.of.interest="variable",
                 methylation=methylation,
                 data=data,
                 generate.confounders="sva",
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
                 method="glm")

ret$table[var$bump.idx,]

ret$table[ret$table$p.adjust < 0.05,]

###################
## test for bumps (there should be one)
bumps.ret <- ewaff.bumps(ret$table$estimate, ret$table$se, ret$table$p.value,
                         methylation,
                         manifest$chr, manifest$pos)
bumps.ret[which(bumps.ret$p.adjust < 0.05),]


