sites.coxph <- function(methylation, design, dependent.variable, ...) {
    formula <- as.formula(paste(dependent.variable, "~ ."))
    design <- as.data.frame(design)
    ret <- rep(NA, 5)    
    stats <- do.call(rbind, mclapply(1:nrow(methylation), function(i) {
        design[,"methylation"] <- methylation[i,]
        try({
            fit <- coxph(formula, design, ...)
            ret <- coef(summary(fit))["methylation",]
            names(ret) <- c("estimate","exp.estimate","se","z","p.value")
            ret[["n"]] <- fit$n
            ret
        }, silent=T)
        ret
    }))
    stats <- as.data.frame(stats)
    stats$n <- floor(stats$n)
    stats
}

