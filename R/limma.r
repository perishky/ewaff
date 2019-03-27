fit.limma <- function(methylation, design, dependent.variable, independent.variable, ...) {
    dep.idx <- which(colnames(design) == dependent.variable)
    fit <- lmFit(methylation, design=design[,-dep.idx,drop=F])
    fit <- eBayes(fit, ...)

    se <- (sqrt(fit$s2.post) * fit$stdev.unscaled[,independent.variable])
    ##alpha <- 0.975
    ##margin.error <- (se * qt(alpha, df=fit$df.total))
    ##estimate.ci <- fit$coefficient[,independent.variable] +/- margin.error
    data.frame(estimate=fit$coef[,independent.variable],
               se=se,
               t=fit$t[,independent.variable],
               p.value=fit$p.value[,independent.variable])
}

fit.limma.multi <- function(methylation,design,dependent.variable,independent.variables,...) {
    dep.idx <- which(colnames(design) == dependent.variable)
    fit <- lmFit(methylation, design=design[,-dep.idx,drop=F])
    fit <- eBayes(fit, ...)

    ## eBayes uses the is.fullrank() function to decide
    ## whether or not to calculate an F-statistic based on a somewhat arbitrary
    ## threshold. There is no warning and omitting the F-statistic
    ## causes topTable() below to fail. Here we compute it if it is missing.
    if (!("F" %in% names(fit))) {
        f.stat <- classifyTestsF(fit, fstat.only = TRUE)
        fit$F <- as.vector(f.stat)
        df1 <- attr(f.stat, "df1")
        df2 <- attr(f.stat, "df2")
        if (df2[1] > 1e+06) 
            fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
        else
            fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
    }
    stats <- topTable(fit,
                      coef=independent.variables,
                      number=nrow(methylation),
                      sort.by="none")

    data.frame(f=stats$F,
               p.value=stats$P.Value,
               sapply(independent.variables, function(independent.variable) {
                   se <- (sqrt(fit$s2.post) * fit$stdev.unscaled[,independent.variable])
                   data.frame(estimate=fit$coef[,independent.variable],
                              se=se,
                              t=fit$t[,independent.variable],
                              p.value=fit$p.value[,independent.variable])
               }, simplify=F))             
}


