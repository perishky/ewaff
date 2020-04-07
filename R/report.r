#' Generate CpG associations report
#'
#' Generate HTML file that summarises CpG site associations. 
#'
#' @param  object Output from \code{\link{ewaff.summary}}.
#' @param  output.file Default = "ewas-report.html".
#' If the file extension is not .htm, .html, .HTM or .HTML then
#' output will be in markdown format.
#' @param  author Default = "Analyst". Author name to be specified on report.
#' @param  study Default = "Illumina methylation data". Study name to be specified on report.
#' @param  ... Arguments to be passed to \code{\link{knitr::knit}}
#' @export
#' @return NULL
#' 
#' @export
ewaff.report <- function(object,
                         output.file = "report.html",
                         author = "Analyst",
                         study = "Illumina methylation data",
                         ...) {
    stopifnot(is.summary.object(object))
    
    ewaff:::msg("Writing report as html file to", output.file)
    report.path <- system.file("report", package="ewaff")
    require(knitr)
    require(Cairo)
    require(gridExtra)

    options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
    
    opts <- opts_chunk$get()
    on.exit(opts_chunk$set(opts))
    opts_chunk$set(warning=FALSE, echo=FALSE, message=FALSE, results="asis",
                   fig.width=6, fig.height=6, dev="CairoPNG")
    knit.report(file.path(report.path, "report.rmd"),output.file, ...)
}

#' Summarize results.
#'
#' Generates variable and covariate summary tables, QQ plots,
#' Manhattan plots, a list of associations, plots of the strongest
#' associations and plots of selected CpG sites.
#'
#' @param object Object returned by \code{\link{ewaff.sites}}.
#' @param chr Chromosome of each CpG site in \code{methylation}.
#' @param pos Position of each CpG site.
#' @param methylation Matrix of methylation levels used in the analysis.
#' @param selected.cpg.sites Vector of CpG site names to plot (Default: character(0)).
#' @param additional.variables A data frame of variables to show associations with the variable of interest.
#' @param parameters Default = ewaff.report.parameters(). List of parameter values.
#' See \code{\link{ewaff.report.parameters}()}.
#' @export
#' @return List
#' 
#' @export
ewaff.summary <- function(object, chr, pos, methylation,
                          selected.cpg.sites=character(0),
                          additional.variables=NULL,
                          parameters=ewaff.report.parameters(),
                          verbose=T) {
    
    stopifnot(is.sites.object(object))

    stopifnot(is.null(additional.variables) || nrow(additional.variables) == ncol(methylation))
    
    if ("f" %in% colnames(object$table)) {
        stop("Summaries and reports for complex models are not yet supported")
    }
    
    if (ncol(methylation) > length(object$sample.idx)) {
        methylation <- methylation[,object$sample.idx]
        if (!is.null(additional.variables))
            additional.variables <- additional.variables[object$sample.idx,,drop=F]
    }
    
    p.values <- object$table$p.value
    p.adjusted <- object$table$p.adjust
    estimates <- object$table$estimate
    
    stopifnot(length(p.values) == length(chr))
    stopifnot(length(p.values) == length(pos))
    stopifnot(length(p.values) == nrow(methylation))
    stopifnot(parameters$max.plots < length(p.values))    
    stopifnot(all(selected.cpg.sites %in% rownames(methylation)))
    
    parameters$practical.threshold <- p.values[order(p.values)[parameters$max.plots+1]]

    if (is.na(parameters$sig.threshold)) 
        parameters$sig.threshold <- 0.05/nrow(methylation)

    sig.idx <- which(p.values < parameters$sig.threshold)
    practical.idx <- which(p.values < parameters$practical.threshold)
    selected.idx <- match(selected.cpg.sites, rownames(methylation))    

    cpg.idx <- union(practical.idx, selected.idx)
    cpg.idx <- cpg.idx[order(chr[cpg.idx], pos[cpg.idx])]

    cpg.stats <- data.frame(chromosome=chr[cpg.idx],
                            position=pos[cpg.idx],
                            estimate=estimates[cpg.idx],
                            p.value=p.values[cpg.idx],
                            p.adjust=p.adjusted[cpg.idx])

    rownames(cpg.stats) <- rownames(methylation)[cpg.idx]
    
    msg("QQ plots", verbose=verbose)
    qq.plot <- ewaff.qq.plot(p.values=p.values,
                              sig.threshold=parameters$sig.threshold,
                              lambda.method=parameters$qq.inflation.method)

    msg("Manhattan plots", verbose=verbose)
    manhattan.plot <- ewaff.manhattan.plot(p.values=p.values,
                                            estimates=estimates,
                                            chr=chr,
                                            pos=pos,
                                            sig.threshold=parameters$sig.threshold)

    plot.sites <- rownames(methylation)[union(practical.idx, selected.idx)]
    msg("CpG site plots:", length(plot.sites), verbose=verbose)
    variable.of.interest <- ifelse(object$independent.variable == "methylation",
                                   object$dependent.variable, object$independent.variable)
    cpg.plots <- sapply(plot.sites, function(cpg) {
        msg("Plotting", cpg, verbose=verbose)
        ewaff.cpg.plot(variable.of.interest, as.data.frame(object$design), methylation[cpg,], cpg)
    }, simplify=F)
    
    sample.characteristics <- NULL
    covariate.associations <- NULL
    additional.associations <- NULL
    if (object$method != "coxph") {
        msg("Sample characteristics", verbose=verbose)
        data <- as.data.frame(object$design[,-1])
        sample.characteristics <- ewaff.sample.characteristics(variable.of.interest, data)
        covariate.associations <- ewaff.covariate.associations(variable.of.interest, data)
        if (!is.null(additional.variables)) {
            additional.variables <- cbind(data[,variable.of.interest,drop=F], additional.variables)
            additional.associations <- ewaff.covariate.associations(variable.of.interest, additional.variables)
        }
    }
    ##parameters$winsorize.pct <- object$winsorize.pct
    
    list(class="ewaff.summary",
         parameters=parameters,
         qq.plot=qq.plot,
         manhattan.plot=manhattan.plot,
         cpg.stats=cpg.stats,
         cpg.plots=cpg.plots,
         practical.sites=rownames(methylation)[practical.idx],
         significant.sites=rownames(methylation)[sig.idx],
         selected.sites=rownames(methylation)[selected.idx],         
         sample.characteristics=sample.characteristics,
         covariate.associations=covariate.associations,
         additional.associations=additional.associations)
}

is.summary.object <- function(object)
    is.list(object) && "class" %in% names(object) && object$class == "ewaff.summary"

#' Specify parameters for report
#'
#' @param sig.threshold P-value threshold for significance (Default: NA).
#' If NA, then threshold used will be 0.05 divided by the number of tests/probes.
#' @param max.plots Maximum number of plots to generate (Default: 10).
#' @param qq.inflation.method Method for calculating genomic inflation lambda.
#' Valid values are "median", "regression" or "robust" (Default: "median").
#' @return List of parameter values
#'
#' @export
ewaff.report.parameters <- function(sig.threshold=NA,
                                    max.plots=10,
                                    qq.inflation.method="median") {
    list(sig.threshold=sig.threshold,
         max.plots=max.plots,
         qq.inflation.method=qq.inflation.method)
}
