#' ewaff.bumps
#'
#' Identify differentially methylated regions from EWAS summary statistics.
#'
#' @param estimate Vector of EWAS effect estimates (corresponds to rows of \code{methylation}).
#' @param se Vector of standard errors of the coefficients.
#' @param p.value Vector of p-values.
#' @param methylation Methylation matrix (rows=features, columns=samples).
#' @param chr Feature chromosome (corresponds to rows of \code{methylation}).
#' @param pos Feature chromosome position.
#' @param p.cutoff Unadjusted p-value cutoff for membership in a candidate hump (Default: 0.05).
#' @param maxgap Maximum distance between consecutive features (Default: 500bp).
#' @param verbose If \code{TRUE} (default), then output status messages.
#' @return A data frame listing all candidate bumps and their summary statistics.
#' 
#' @examples
#' chr <- ... ## chromosome of each CpG site
#' pos <- ... ## position of each CpG site
#' data <- ... ## data frame containing variable of interest and covariates
#' methylation <- ... ## methylation matrix
#' ewas.ret <- ewaff.sites(methylation ~ variable + .,
#'                       variable.of.interest="variable",
#'                       methylation=methylation,
#'                       data=data)
#'
#' bumps <- ewaff.bumps(estimate=ewas.ret$table$estimate, se=ewas.ret$table$se, p.value=ewas.ret$table$p.value,
#'                      methylation=methylation,
#'                      chr=chr, pos=pos)
#'                      
#' bumps[which(bumps$p.adjust < 0.05),
#'       c("chr","start","end","n", "z", "p.value","p.adjust")]
#'
#' @export
ewaff.bumps <- function(estimate, se, p.value, methylation, chr, pos, maxgap=500, p.cutoff=0.05, verbose=T) {
    stopifnot(is.vector(estimate))
    stopifnot(is.vector(se))
    stopifnot(is.vector(p.value))
    stopifnot(is.matrix(methylation))
    stopifnot(length(estimate) == length(se))
    stopifnot(length(estimate) == nrow(methylation))
    stopifnot(length(estimate) == length(p.value))
    stopifnot(length(estimate) == length(chr))
    stopifnot(length(estimate) == length(pos))
    
    candidates <- ewaff.bump.candidates(estimate=estimate,
                                        p.value=p.value,
                                        chr=chr, 
                                        pos=pos,
                                        maxgap=maxgap,
                                        p.cutoff=p.cutoff,
                                        verbose=verbose)

 
    stats <- shrink.candidates(candidates$start.idx, candidates$end.idx,
                               function(start.idx,end.idx) {
                                   idx <- start.idx:end.idx
                                   ivwfe.getz(estimate[idx], se[idx], methylation[idx,,drop=F])
                               })

    collate.bump.stats(stats, chr, pos)
}


collate.bump.stats <- function(stats, chr, pos) {   
    stats <- with(stats, data.frame(chr=chr[start.idx],
                                    start=pos[start.idx],
                                    end=pos[end.idx],
                                    n=end.idx-start.idx+1,
                                    start.idx=start.idx,
                                    end.idx=end.idx,
                                    start.orig=start.orig,
                                    end.orig=end.orig,
                                    p.orig=2*pnorm(-abs(z.orig), lower.tail=T),
                                    p.value=2*pnorm(-abs(z), lower.tail=T)))
    number.tests <- length(chr) + calculate.number.shrink.tests(stats)
    stats$p.adjust <- pmin(1, stats$p.value * number.tests)
    stats
}
