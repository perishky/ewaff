#' ewaff.meta.bumps
#'
#' Identify differentially methylated regions from meta-analyzed summary statistics.
#'
#' @param meta.estimate Vector of meta-analysed effect estimates.
#' @param meta.se Vector of meta-analysed errors of the coefficients.
#' @param estimate Matrix of effect estimates (rows=features, cols=analyses).
#' @param se Matrix of standard errors of the coefficients (rows=features, cols=analyses).
#' @param methylation Reference methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome.
#' @param pos Feature chromosome position.
#' @param p.cutoff Unadjusted p-value cutoff for membership in a candidate hump (Default: 0.05).
#' @param maxgap Maximum distance between consecutive features (Default: 500bp).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing all candidate bumps and their meta-analysed statistics.
#'
#' @examples
#' estimates <- ...   ## rows=features, columns=analyses
#' ses <- ...         ## rows=features, columns=analyses
#' methylation <- ... ## rows=features, columns=samples
#' chr <- ...         ## feature chromosome
#' pos <- ...         ## feature position
#' 
#' ma.ret <- ewaff.meta.sites(estimates, ses, method="FE")
#'
#' bumps <- ewaff.meta.bumps(ma.ret$estimate, ma.ret$se, estimates, ses, methylation, chr, pos)
#'
#' bumps[which(bumps$p.adjust < 0.05),
#'       c("chr","start","end","n", "z", "p.value","p.adjust")]

#' @export
ewaff.meta.bumps <- function(meta.estimate, meta.p.value,
                             estimate, se, methylation, chr, pos, maxgap=500, p.cutoff=0.05, verbose=T) {

    stopifnot(is.vector(meta.estimate))
    stopifnot(is.vector(meta.p.value))
    stopifnot(is.matrix(estimate))
    stopifnot(is.matrix(se))
    stopifnot(is.matrix(methylation))
    stopifnot(length(meta.estimate) == length(meta.p.value))
    stopifnot(length(meta.estimate) == nrow(methylation))
    stopifnot(length(meta.estimate) == nrow(estimate))
    stopifnot(length(meta.estimate) == nrow(se))
    stopifnot(ncol(estimate) == ncol(se))
    
    candidates <- ewaff.bump.candidates(estimate=meta.estimate,
                                        p.value=meta.p.value,
                                        chr=chr, 
                                        pos=pos,
                                        maxgap=maxgap,
                                        p.cutoff=p.cutoff,
                                        verbose=verbose)

    stats <- shrink.candidates(candidates$start.idx, candidates$end.idx,
                               function(start.idx,end.idx) {
                                   idx <- start.idx:end.idx
                                   rho <- ivwfe.rho(methylation[idx,,drop=F])
                                   z <- sapply(1:ncol(estimate), function(j) {
                                       ivwfe.getz(estimate[idx,j], se[idx,j], rho=rho)
                                   })                                   
                                   sum(z, na.rm=T)/sqrt(ncol(estimate))
                               })

    collate.bump.stats(stats, chr, pos)
}
