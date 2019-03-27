#' Summarize CpG distributions
#'
#' Summarize CpG site quantiles (0, 10, 15, 50, 75, 90, 100), the 
#' mean, and SD. 
#'
#' @param methylation DNA methylation matrix, one row per CpG site,
#' one column per sample.
#' @param cpgs.are.rows logical argument to allow CpG sites to be in 
#' columns instead of rows. 
#' @param ... Arguments to \code{\link{glm}}.
#' @return matrix with columns of each requested summary statistics 
#' and row for each CpG in \code{methylation}
#'
#' @export
ewaff.cpg.summary <- function(methylation, cpgs.are.rows=T){
    if (!cpgs.are.rows)
    	methylation <- t(methylation)

    t(apply(methylation, 1, function(object){
	    nas <- is.na(object)
	    object <- object[!nas]
	    qq <- quantile(object, c(0.00, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00))
	    qq <- signif(c(qq[1L:3L], mean(object), qq[4L], sd(object), qq[5L:7L]), 
	    			max(3, getOption("digits")-3))
	    names(qq)[c(1,4:6,9)] <- c("Min.", "Mean", "Median", "SD", "Max.")
		c(qq, `NA's` = sum(nas), `N` = length(object) - sum(nas) )
    }))
}


lambda <- function(p.values){
	qunif(median(p.values, na.rm=T), lower.tail = F)/ 
  		qunif(0.5, 1) 
}

lambda.chi <- function(p.values){
	qchisq(median(p.values, na.rm=T), df = 1, lower.tail = F)/ 
 		qchisq(0.5, 1) 
}

 