#' Test associations at CpG sites
#'
#' Fit generalized linear model (GLM) to methylation levels each CpG site.
#'
#' Note: \code{\link{mclapply}} is used to fit regression models
#' using multiple cores.  
#' 
#' @param formula An object of class \code{\link{formula}}:
#' a symbolic description of the model to be fitted.
#' DNA methylation is referred to by the variable name 'methylation'.
#' @param methylation DNA methylation matrix, one row per CpG site,
#' one column per sample.
#' @param variable.of.interest Name of variable(s) in the model formula
#' for which to save summary statistics.  If it is the dependent variable
#' in the formula, then it must be numeric or binary; otherwise, it may
#' be any type of variable or even a vector of variables.
#' The value is ignored if \code{method == "coxph"}.
#' @param data Data frame of variables to include in the model.
#' @param family See description for \code{\link{glm}}.
#' @param method Method for regressions: "glm", "rlm", "limma" or "coxph" (Default: "glm").
#' @param generate.confounders Generate variables from the methylation data
#' to adjust for unknown confounders. May be \code{NULL} for none (default),
#' "sva" or "smartsva" to generate surrogate variables, or "pca" to generate
#' prinicipal components. If \code{method == "coxph"}, then \code{generate.confounders}
#' must be \code{NULL} or \code{"pca"}.
#' @param n.confounders Number of unknown confounders to generate.  A value of
#' \code{NULL} allows \code{\link{sva}} to select the number automatically
#' or if principal components the maximum number (Default: NULL).
#' @param most.variable Generate confounders from the
#' #' given most variable CpG sites rather than the whole matrix (Default: NULL).
#' @param random.subset Generate surrogate variables from the
#' given percentage of randomly selected CpG sites rather than the whole matrix (Default: 0.05, i.e. 5 percent).
#' @param ... Arguments to \code{\link{glm}}, \code{\link{rlm}}, \code{\link{limma::eBayes}},
#' or \code{\link{survival::coxph}}.
#' @return List containing a table of association statistics ('table') and the
#' model design matrix ('design').
#'
#' @export
ewaff.sites <- function(formula,
                        variable.of.interest,
                        methylation,
                        data,
                        family=gaussian,
                        method="glm",
                        generate.confounders=NULL,
                        n.confounders=NULL,
                        most.variable=NULL,
                        random.subset=0.05,
                        ...) {

    stopifnot(method %in% c("glm","rlm","limma","coxph"))
    stopifnot(is.null(generate.confounders)
              || generate.confounders == "sva"
              || generate.confounders == "smartsva"
              || generate.confounders == "pca"
              && !is.null(n.confounders) && n.confounders > 0 && n.confounders < ncol(methylation))
    stopifnot(method != "coxph" || is.null(generate.confounders) || generate.confounders=="pca")

    if (method == "coxph")
        variable.of.interest <- as.character(formula)[2]
    
    design <- build.design.matrix(formula=formula,
                                  variable.of.interest=variable.of.interest,
                                  methylation=methylation,
                                  data=data,
                                  family=family,
                                  method=method,
                                  generate.confounders=generate.confounders,
                                  n.confounders=n.confounders,
                                  most.variable=most.variable,
                                  random.subset=random.subset,
                                  ...)

    ## get names of the dependent variable
    ## and the independent variable of interest
    dependent.variable <- design$dependent.variable
    independent.variable <- design$independent.variable
    family <- design$family
    sample.idx <- design$sample.idx
    design <- design$matrix
    
    methylation <- methylation[,sample.idx,drop=F]

    stopifnot(method %in% c("glm","coxph") || dependent.variable == "methylation")

    if (length(independent.variable) == 0 || length(dependent.variable) == 0)
        stop("The variable.of.interest is missing")

    if (length(dependent.variable) > 1)
        stop("When methylation is an independent variable, variable.of.interest must name a single numeric or binary variable")

    if (length(independent.variable) > 1 && !(method %in% c("glm","limma")))
        stop("When method is not 'limma' or 'glm', variable.of.interest must be a single numeric or binary variable")

    if (method != "coxph")
        stats <- sites.glm(methylation, design, dependent.variable, independent.variable, method, family, ...)
    else
        stats <- sites.coxph(methylation, design, dependent.variable, ...)

    stats$p.adjust <- p.adjust(stats$p.value, "bonferroni")
    rownames(stats) <- rownames(methylation)

    list(class="sites",
         formula=paste(as.character(formula)[c(2,1,3)], collapse=" "),
         variable.of.interest=variable.of.interest,
         data=data,
         family=family,
         method=method,
         generate.confounders=generate.confounders,
         n.confounders=n.confounders,
         dependent.variable=dependent.variable,
         independent.variable=independent.variable,
         design=design[,-which(colnames(design) == "methylation"),drop=F],
         sample.idx=sample.idx,
         table=stats)
}

is.sites.object <- function(object)
    is.list(object) && "class" %in% names(object) && object$class == "sites"

