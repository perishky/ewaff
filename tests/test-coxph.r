library(ewaff)
options(mc.cores=4)

source("simulation-functions.r")
library(simsurv)

###################################
## construct a random dataset with
## a single CpG site association
set.seed(20180220)

n <- 100 ## n samples
s <- 1000  ## s features

## generate survival data correlated with random methylation pattern
data <- data.frame(binvar=sample(c("A","B"), n, replace=T),
                   convar=rnorm(n),
                   methylation=runif(n))
data <- cbind(data,
              simsurv(lambdas=0.1,
                      gammas=1.5,
                      betas=c(methylation = 2),
                      x=data)[,-1])

manifest <- data.frame(chr=c(rep(1,s/2), rep(2,s/2)), pos=sample(1:(150*s), s))
manifest <- manifest[order(manifest$chr, manifest$pos),]

methylation <- generate.methylation(n, manifest$pos)
rownames(methylation) <- paste("cpg", 1:nrow(methylation), sep="")

## insert the survival correlated methylation pattern
## at the location of the random CpG site most strongly
## correlated with the methylation pattern
r <- abs(cor(t(methylation), data$methylation))
methylation[which.max(r),] <- data$methylation

data$methylation <- NULL

############################
## add some missing values

methylation[sample(1:length(methylation), 5)] <- NA

data$binvar[sample(1:nrow(data), 3)] <- NA

###################################
## test associations at each CpG site
sites.ret <- ewaff.sites(Surv(eventtime, status) ~ methylation + binvar + convar,
                         methylation=methylation,
                         data=data,
                         generate.confounders="pca",
                         n.confounders=5,
                         random.subset=0.5,
                         method="coxph")

##########################
## generate report
sum.ret <- ewaff.summary(sites.ret, manifest$chr, manifest$pos, methylation,
                         selected.cpg.sites=paste("cpg", which.max(r)+ -5:5, sep=""))

ewaff.report(sum.ret, output.file="output/coxph.html",
             author="Gloria Gaynor",
             study="Surviving stuff")



