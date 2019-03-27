#' Handle outliers
#'
#' @param methylation DNA methylation matrix, one row per CpG site,
#' one column per sample.
#' @param method The method for handling outliers.  Options are 'iqr' (default) and 'winsorize'.
#' The 'iqr' option sets observations beyond a user-defined mutliple of the inter-quartile
#' (\code{iqr.limit}) to \code{NA}.  The 'winsorize' option sets observations beyond
#' the \code{winsorize.pct} or \code{1-winsorize.pct} percentiles to the corresponding
#' percentile.
#' @param iqr.limit The multiple of IQR beyond which observations
#' should be set to NA. Default is 3.
#' @param winsorize.pct The percentile used for winsorizing observations. Default is 0.05.
#' @return List containing the methylation matrix with outliers modified
#' depending on the method selected, and a log of how many observations were
#' identified for each probe.
#' 
#' @export
ewaff.handle.outliers <- function(methylation, method=c("iqr","winsorize"), iqr.limit=3, winsorize.pct=0.05) {
    
    method <- match.arg(method)
    
    stopifnot(is.matrix(methylation))
    if(nrow(methylation) < ncol(methylation))
        warning("expecting CpG methylation as rows (long dataset)") 

    if (method == "iqr") {
        iqr.trim(methylation, iqr.limit)
    }
    else {
        winsorize(methylation, winsorize.pct)
    }
    
}

iqr.trim<-function(methylation, iqr.limit=3) { 
    ## find outlying observations more extreme that the iqr.limit
    quantiles <- matrixStats::rowQuantiles(methylation, probs = c(0.25, 0.75), na.rm = T)
    iqr <- quantiles[,2] - quantiles[,1]
    low <- quantiles[,1]
    upper <- quantiles[,2]

    initial.missing <- rowSums(is.na(methylation))

    methylation[methylation < low - iqr.limit * iqr] <- NA
    outliers.lower <- rowSums(is.na(methylation)) - initial.missing

    methylation[methylation > upper + iqr.limit * iqr] <- NA
    outliers.upper <- rowSums(is.na(methylation)) - initial.missing - outliers.lower

    n <- ncol(methylation) - initial.missing - outliers.upper - outliers.lower
    log <- data.frame(outliers.lower, outliers.upper, n) 
    
    return(list(methylation=methylation, log=log)) 
} 


winsorize <- function(methylation,pct=0.05) {
    quantiles <- matrixStats::rowQuantiles(methylation, probs=c(pct,1-pct), na.rm=T)
    low <- quantiles[,1]
    upper <- quantiles[,2]

    outliers.lower <- rowSums(methylation < low, na.rm=T)
    outliers.upper <- rowSums(methylation > upper, na.rm=T)
    
    idx <- which(methylation < low, arr.ind=T)
    methylation[idx] <- low[idx[,1]]
    
    idx <- which(methylation > upper, arr.ind=T)
    methylation[idx] <- upper[idx[,1]]

    n <- rowSums(!is.na(methylation))
    log <- data.frame(outliers.lower, outliers.upper, n)
    
    return(list(methylation=methylation, log=log))
}
