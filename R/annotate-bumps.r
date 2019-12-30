#' ewaff.annotate.bumps
#'
#' Annotate bumps.
#'
#' @param bumps Output from \code{\link{ewaff.bumps}()}.
#' @param annotation A vector annotating the features
#' in the input of \link{ewaff.bumps}.
#' @return A vector providing the annotation for each bump.
#'
#' methylation <- ... ## methylation matrix
#' annotation <- ... ## data frame with gene annotations for each gene
#' bumps <- ewaff.bumps(..., methylation=methylation, ...)
#' idx <- match(rownames(methylation), annotation$cpg)
#' bumps$gene <- ewaff.annotate.bumps(bumps, annotation$gene[idx])
#' 
#' @export
ewaff.annotate.bumps <- function(bumps, annotation) {
    sapply(1:nrow(bumps), function(i) {
        start.idx <- bumps$start.idx[i]
        end.idx <- bumps$end.idx[i]
        values <- annotation[start.idx:end.idx]
        unique.values <- unique(values)
        unique.values[which.max(tabulate(match(values, unique.values)))]
    })
}
