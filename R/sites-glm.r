sites.glm <- function(methylation, design, dependent.variable, independent.variable, method, family, ...) {
    if (method == "rlm") {
        stats <- fit.loop(methylation, design, dependent.variable, independent.variable,
                          fit.rlm, stats.rlm, maxit=200, ...)
    }
    else if (method == "glm") {
        if (length(independent.variable) == 1)
            stats <- fit.loop(methylation, design, dependent.variable, independent.variable,
                              fit.glm, stats.glm, family=family, ...)
        else
            stats <- fit.multi.loop(methylation, design, dependent.variable, independent.variable,
                                    fit.glm, stats.glm, family=family, ...)
    }
    else {
        if (length(independent.variable) == 1)
            stats <- fit.limma(methylation, design, dependent.variable, independent.variable, ...)
        else
            stats <- fit.limma.multi(methylation, design, dependent.variable, independent.variable, ...)
    }
    stats <- as.data.frame(stats)
    stats$n <- rowSums(!is.na(methylation))
    stats
}



fit.rlm <- function(x, y, ...) {
    fit <- rlm(x=x, y=y, ...)
}

stats.rlm <- function(fit, var) {
    ret <- coeftest(fit, vcov=vcovHC(fit, type="HC0"))[var,]
    names(ret) <- c("estimate","se","z","p.value")
    ret
}

fit.glm <- function(x, y, ...) {
    fit <- glm.fit(x=x, y=y, ...)
}

stats.glm <- function(fit, var) {
    ret <- coef(summary.glm(fit))[var,]
    names(ret) <- c("estimate","se","t","p.value")
    ret
}



