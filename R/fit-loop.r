fit.loop <- function(methylation, design, dependent.variable, independent.variable, fit.function, stats.function, ...) {
    ret <- rep(NA,4)
    dep.idx <- which(colnames(design) == dependent.variable)

    do.call(rbind, mclapply(1:nrow(methylation), function(i) {
        design[,"methylation"] <- methylation[i,]   ## assign CpG site methylation
        idx <- which(!is.na(design[,"methylation"]))
        try({
            fit <- fit.function(x=design[idx,-dep.idx,drop=F], ## independent variables
                                y=design[idx,dep.idx],         ## dependent variable
                                ...)
            ret <- stats.function(fit, independent.variable) 
        }, silent=T)
        ret
    }))
}
