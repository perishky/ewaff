#' Generate surrogate variables using the ewaff.sites syntax
#'
#' Provide a function to manually compute surrogate variables, using  
#' the same syntax for submission of a ewaff.sites call
#' 
#' Note: surrogate variables are computed using \code{\link{sva}}
#' 
#' @param formula An object of class \code{\link{formula}}:
#' a symbolic description of the model to be fitted.
#' DNA methylation is referred to by the variable name 'methylation'.
#' @param methylation DNA methylation matrix, one row per CpG site,
#' one column per sample.
#' @param variable.of.interest Name of variable in the GLM model formula
#' for which to save summary statistics.
#' @param data Data frame providing variables referenced in model formula.
#' @param family (Default: gaussian).
#' @param n.sv Number of surrogate variables to compute (Default: 10).
#' @param most.variable Generate surrogate variables from the
#' given most variable CpG sites rather than the whole matrix (Default: NULL).
#' @param random.subset Generate surrogate variables from the
#' given percentage of randomly selected CpG sites rather than the whole matrix (Default: 0.05, i.e. 5 percent).
#' @param algorithm Algorithm to use, either "sva" or "smartsva" (Default: "sva").
#' @return Matrix containing surrogate variables computed
#'
#' @export
ewaff.sva <- function(formula,
                      variable.of.interest,
                      methylation,
                      data,
                      family=gaussian,                      
                      n.sv=10,
                      most.variable=NULL,
                      random.subset=0.05,
                      algorithm="sva",
                      ...) {
    design <- build.design.matrix(formula=formula,
                                  variable.of.interest=variable.of.interest,
                                  methylation=methylation,
                                  data=data,
                                  family=family,
                                  method="glm",
                                  generate.confounders=algorithm,
                                  n.confounders=n.sv,
                                  most.variable=most.variable,
                                  random.subset=random.subset,
                                  ...)$matrix
    design[,grep("^sv[0-9]+", colnames(design)),drop=F]
}
