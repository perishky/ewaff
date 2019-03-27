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
methylation <- t(sapply(r, function(r) rcor(as.numeric(data$variable), r))) 

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
                   method="glm")

## list summary statistics
ret$table

## check association statistics of the first cpg site manually
coef(summary(glm(methylation[1,] ~ ., data=as.data.frame(ret$design))))["variableB",]


ret2 <- ewaff.sites(methylation ~ variable + .,
                 variable.of.interest="variable",
                 methylation=methylation,
                 data=data,
                 generate.confounders="sva",
                 method="limma")
ret2$table



ret3 <- ewaff.sites(methylation ~ variable + .,
                 variable.of.interest="variable",
                 methylation=methylation,
                 data=data,
                 generate.confounders="sva",
                 method="rlm")

ret3$table

###############################################
## perform another EWAS with the variable of interest as dependent variable
## so use logistic regression
ret4 <- ewaff.sites(variable ~ methylation + .,
                  variable.of.interest="variable",
                  methylation=methylation,
                  data=data,
                  family="binomial",
                  generate.confounders="sva",
                  method="glm")

## list summary statistics
ret4$table

## manually verify statistics for the first cpg site 
coef(summary(glm(variable ~ methylation[1,] + ., data=as.data.frame(ret4$design), family="binomial")))["methylation[1, ]",]



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
                  method="glm")

ret5$table

ret6 <- ewaff.sites(methylation ~ .,
                  variable.of.interest=c("variable1","variable2"),
                  methylation=methylation,
                  data=data,
                  family=gaussian, 
                  generate.confounders="sva",
                  method="limma")

ret6$table


cbind(ret5$table$variable1B.t, ret6$table$variable1B.t)

cbind(ret5$table$f, ret6$table$f)


