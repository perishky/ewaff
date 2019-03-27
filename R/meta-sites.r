#' ewaff.meta.sites
#'
#' Meta-analyse EWAS summary statistics.
#'
#' Note: \code{\link{mclapply}} is used to fit meta-analysis models
#' using multiple cores.  
#' 
#' @param estimate Matrix of regression estimates (rows=features, columns=samples).
#' @param se Matrix of standard errors of the estimates (rows=features, columns=samples).
#' @param ... Arguments to \code{\link{metafor::rma.uni}}.
#' @return Data frame containing meta-analyzed statistics and heterogeneity measures.
#' 
#' @export
ewaff.meta.sites <- function(estimate, se, ..., verbose=T) {
    stopifnot(is.matrix(estimate) && is.matrix(se))
    stopifnot(nrow(estimate) == nrow(se))
    stopifnot(ncol(estimate) == ncol(se))
    
    ret <- do.call(rbind, mclapply(1:nrow(estimate), function(i) {
        msg(i, verbose=(i %% 10000 == 1 && verbose))
        fit <- rma.uni(yi=estimate[i,], sei=se[i,], ...)
        c(estimate=fit$b["intrcpt",1],
          se=fit$se,
          z=fit$zval,
          p.value=fit$pval,
          ci.ub=fit$ci.ub,
          ci.lb=fit$ci.lb,
          Q=fit$QE,
          Q.p=fit$QEp,
          I2=fit$I2,
          H2=fit$H2,
          tau2=fit$tau2)
    }))
    rownames(ret) <- rownames(estimate)
    colnames(ret)[1] <- "estimate"
    as.data.frame(ret)
}
