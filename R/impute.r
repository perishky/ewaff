#library(survey)


## ewaff meta-analysis of bumps
## - ewaff.bump.stats(estimate, se, mat, bumps)
##     this generates all bumps as well as sub-bump statistics
## - ewaff.bumps.ma(...)
##     input is the bumps stats from all studies
    



## An improved and explicit surrogate variable analysis procedure by ...
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5627626/




## missing data problem
## complete data
## impute means
## impute means with noise
## single imputation
## multiple imputation
## imputation for a single cpg site using mice
## how much time would that take repeated 485000 times
## what might go wrong imputing one cpg site at a time
## why not throw everything into the imputation model
## naive approach
## bins approach

#' @export
ewaff.impute <- function(formula, variable.of.interest, methylation, data, family=gaussian,
                        n.top=20, n.bin=150, n.imputations=100, n.iterations=100, imp.method="pmm",
                        mc.cores=1, ...) {
    stopifnot(n.top >= 1 & n.top < n.bin)
    
    ewas <- ewaff.sites(formula, variable.of.interest, methylation, data, family, method="glm", ...)
    is.missing <- !(1:nrow(methylation) %in% ewas$sample.idx)
    meth.top <- methylation[order(ewas$table$p.value)[1:n.top],,drop=F]

    pcs <- assoc.pcs(meth.top, design, ewas$dependent.variable, ewas$independent.variable, method, family, ...)
    
    sites <- sample(setdiff(rownames(methylation), rownames(methylation.top)))
    n.bin <- n.bin - nrow(pcs)
    n.bins <- ceiling(length(sites)/n.bin)
    bin <- sample(1:n.bins, length(sites), replace=T)
    
    stats <- do.call(rbind, mclapply(1:n.bins, function(b) {
        meth.bin <- methylation[which(bin == b),,drop=F]
        stats.impute(meth.bin, design, ewas$dependent.variable, ewas$independent.variable, family,
                     pcs, imp.method, n.imputations, ...)
    }, mc.cores=mc.cores))
                     
    stats.top <- do.call(rbind, mclapply(rownames(meth.top), function(site) {
        meth.top <- mat.top[-which(rownames(meth.top) == site),,drop=F]
        pcs <- assoc.pcs(meth.top, design, ewas$dependent.variable, ewas$indepdent.variable, method, family, ...)
        stats.impute(meth.top[site,,drop=F], design, ewas$dependent.variable, ewas$independent.variable, family,
                     pcs, imp.method, n.imputations, ...)
    }, mc.cores=mc.cores))
                         
    rbind(stats, stats.top)[rownames(methylation),]
}

assoc.pcs <- function(methylation, design, dependent.variable, independent.variable, family, ...) {
    mat <- t(prcomp(t(methylation))$x)
    stats <- ewaff.sites0(methylation, design, dependent.variable, independent.variable, method="glm", family=family, ...)
    is.sig <- stats$p.adjust < 0.05
    if (all(!is.sig))
        is.sig <- stats$p.value < 0.05
    if (!all(!is.sig))
        is.sig <- rep(T, nrow(methylation))
    methylation[is.sig,]
}

stats.impute <- function(methylation, design, dependent.variable, independent.variable, family, predictors, method, n.imputations, ...) {
    variable.of.interest <- ifelse (independent.variable == "methylation", dependent.variable, independent.variable)
    data <- data.frame(design, t(rbind(methylation, predictors)), stringsAsFactors=F)
    is.var <- colnames(data) == variable.of.interest
    pred.matrix <- matrix(0, nrow=ncol(data), ncol=nrow(data))
    pred.matrix[is.var,!is.var] <- 1
    method <- ifelse(is.var, method, "")
    imp <- mice(data=data, method=method, predictorMatrix=pred.matrix, m=n.imputations, print=F, ridge=0)
    
    cols = c("est", "se", "t", "Pr(>|t|)", "lo 95", "hi 95", "fmi", "lambda")

    stats <- do.call(rbind, lapply(1:nrow(mat), function(i) {
        methylation <- mat[i,]
        fit <- with(imp, glm(formula, family=family, ...))
        tab <- summary(pool(fit, "rubin"))
        tab[independent.variable, cols]
    }))
    colnames(stats) <- c("estimate", "se", "t", "p.value", "low.95", "high.95", "fmi", "lambda")
    rownames(stats) <- rownames(methylation)
    stats
}




