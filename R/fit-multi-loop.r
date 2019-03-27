fit.multi.loop <- function(methylation, design, dependent.variable, independent.variables, fit.function, stats.function, ...) {
    ret <- rep(NA,2 + length(independent.variables)*4)
    dep.idx <- which(colnames(design) == dependent.variable)
    ind.idx <- which(colnames(design) %in% independent.variables)

    do.call(rbind, mclapply(1:nrow(methylation), function(i) {
        design[,"methylation"] <- methylation[i,]   ## assign CpG site methylation
        idx <- which(!is.na(design[,"methylation"]))
        try({
            fit <- fit.function(x=design[idx,-dep.idx,drop=F], ## independent variables
                                y=design[idx,dep.idx],         ## dependent variable
                                ...)
            fit0 <- fit.function(x=design[idx,-c(dep.idx,ind.idx)],
                                 y=design[idx,dep.idx],
                                 ...)
            
            delta.dev <- abs(fit0$deviance-fit$deviance)
            delta.df <- abs(fit0$df.residual-fit$df.residual)
            dispersion <- fit$deviance/fit$df.residual
            f.stat <- (delta.dev/delta.df)/dispersion
            p.value <- pf(f.stat, delta.df, fit$df.residual, lower.tail=FALSE)
            ## derived from stats:::anova.glm, stats:::stat.anova, gaussian()$dev.resids, summary.glm
            ret <- c(f=f.stat, p.value=p.value,
                     unlist(sapply(independent.variables, function(independent.variable) {
                         stats.function(fit, independent.variable)
                     }, simplify=F)))
        }, silent=T)
        ret
    }))
}

