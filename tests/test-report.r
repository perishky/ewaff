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

manifest <- data.frame(chr=c(rep(1,s/2), rep(2,s/2)), pos=sample(1:(150*s), s))
manifest <- manifest[order(manifest$chr, manifest$pos),]

methylation <- generate.methylation(n, manifest$pos)
rownames(methylation) <- paste("cpg", 1:nrow(methylation), sep="")

########################################
## construct variable with associations
var <- generate.true.bump.var(methylation, manifest$chr, manifest$pos, cluster.sites=20,
                              bump.sites=10, cluster.position=0.5, r=0.5, maxgap=500)
data$bump <- var$var


############################
## add some missing values

methylation[sample(1:length(methylation), 5)] <- NA

data$bump[sample(1:nrow(data), 2)] <- NA
data$variable[sample(1:nrow(data), 3)] <- NA

###################################
## test associations at each CpG site
sites.ret <- ewaff.sites(methylation ~ .,
                   variable.of.interest="bump",
                   methylation=methylation,
                   data=data,
                   generate.confounders="sva",
                   random.subset=0.9,
                   method="glm")

##########################
## generate report
sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$pos, methylation,
                         selected.cpg.sites=paste("cpg", var$bump.idx, sep=""))

ewaff.report(sum.ret, output.file="output/report.html",
             author="Robin Banks",
             study="Mining bitcoin")
