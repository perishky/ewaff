build.design.matrix  <- function(formula,
                                 variable.of.interest,
                                 methylation,
                                 data,
                                 family=gaussian,
                                 method="glm",
                                 generate.confounders=NULL,
                                 n.confounders=NULL,
                                 most.variable=NULL,
                                 ...) {
    stopifnot(is.null(generate.confounders)
              || generate.confounders == "sva"
              || generate.confounders == "smartsva"
              || generate.confounders == "pca"
              && !is.null(n.confounders) && n.confounders > 0 && n.confounders < ncol(methylation))

    stopifnot(is.null(most.variable)
              || (is.numeric(most.variable)
                  && most.variable > 1 
                  && most.variable <= nrow(methylation)))
    
    stopifnot(method != "coxph" || is.null(generate.confounders) || generate.confounders=="pca")

    if (!is.data.frame(data)) {
        stopifnot(is.matrix(data))
        data <- as.data.frame(data, stringsAsFactors=F)
    }
    
    ## fit model for one CpG site to get a fit object
    ## from which we can extract the design matrix corresponding
    ## to the formula
    data$methylation <- runif(nrow(data))
    if (method != "coxph") {
        fit <- glm(formula, data=data, family=family, x=TRUE, y=TRUE, ...)
        family <- fit$family

        ## get names of the dependent variable
        ## and the independent variable of interest    
        if ("methylation" %in% colnames(fit$x)) {
            dependent.variable <- variable.of.interest
            independent.variable <- "methylation"
        } else {
            dependent.variable <- "methylation"
            variable.of.interest <- unlist(lapply(variable.of.interest, function(varname) {               
                if (!is.numeric(data[[varname]])) {
                    if (!is.factor(data[[varname]]))
                        data[[varname]] <- as.factor(data[[varname]])
                    varname <- paste(varname, levels(data[[varname]]), sep="")
                }
                intersect(varname, colnames(fit$x))
            }))
            independent.variable <- variable.of.interest
        }

        if (length(independent.variable) == 0)
            stop(variable.of.interest, " is missing")
        
        if (length(dependent.variable) > 1)
            stop("When methylation is an independent variable, variable.of.interest must be a single numeric or binary variable")

            ## extract the design matrix from glm fit object
        design <- cbind(fit$x, fit$y)
        colnames(design) <- c(colnames(fit$x), dependent.variable)
    } else {
        fit <- coxph(formula, data=data, x=TRUE, y=TRUE, ...)
        family <- NA
        independent.variable <- "methylation"
        dependent.variable <- variable.of.interest
        design <- cbind(fit$x, fit$y)
        colnames(design) <- c(colnames(fit$x), extract.survival.variables(formula))
    }
    
    observation.names <- rownames(fit$x)                
    sample.idx <- match(rownames(fit$x), rownames(data))
        
    if (length(sample.idx) < 2)
        stop("Too many missing values in 'data'")
    
    methylation <- methylation[,sample.idx,drop=F]
    data <- data[sample.idx,,drop=F]
        
    if (!is.null(generate.confounders)) {
        ## replace missing values with CpG methylation mean 
        na.idx <- which(is.na(methylation), arr.ind=TRUE)
        if (nrow(na.idx) > 0) {
            mean.values <- rowMeans(methylation[na.idx[,1],,drop=F], na.rm=TRUE)
            nan.idx <- which(is.nan(mean.values))
            if (length(nan.idx) > 0) ## if cpg site is entirely missing
                mean.values[nan.idx] <- 0.5
            methylation[na.idx] <- mean.values
        }
        
        ## retain the most variable CpG sites
        if (!is.null(most.variable)) {
            var.idx <- order(rowVars(methylation, na.rm=T), decreasing=T)[1:most.variable]
            methylation <- methylation[var.idx,,drop=F]
        }
        
        ## perform surrogate variable analysis and add SVs to the design matrix
        if (generate.confounders %in% c("sva","smartsva")) {
            ## prepare complete and null models for SVA
            mod <- design[,-which(colnames(design) == "methylation"),drop=F]
            mod0 <- mod[,-which(colnames(mod) %in% variable.of.interest),drop=F]
                        
            ## perform SVA
            if (generate.confounders == "smartsva") {
                if (is.null(n.confounders)) {
                    res <- t(resid(lm(t(methylation) ~ ., data=as.data.frame(mod))))
                    n.confounders <- EstDimRMT(res, FALSE)$dim + 1
                }               
                sva.ret <- smartsva.cpp(methylation, mod=mod, mod0=mod0, n.sv=n.confounders)
            }
            else
                sva.ret <- sva(methylation, mod=mod, mod0=mod0, n.sv=n.confounders)
                        
            if (sva.ret$n.sv > 0) {
                ## add SVs to the design matrix
                svs <- sva.ret$sv
                if (sva.ret$n.sv == 1) svs <- matrix(svs, ncol=1)
                colnames(svs) <- paste("sv", 1:ncol(svs), sep="")
                design <- cbind(design, svs)
            }
        } else if (generate.confounders == "pca") {
            pcs <- prcomp(t(methylation), scale=T, center=T)$x[,1:n.confounders,drop=F]
            p.values <- sapply(1:ncol(pcs), function(i) {
                data$methylation <- pcs[,i]
                if (method != "coxph") {
                    fit <- glm(formula, data=data, family=family, ...)
                    coef(summary(fit))[independent.variable,"Pr(>|t|)"]
                } else {
                    fit <- coxph(formula, data=data, ...)
                    coef(summary(fit))[independent.variable,"Pr(>|z|)"]
                }
            })
            pcs <- pcs[,which(p.values > 0.05),drop=F]
            if (nrow(pcs) > 0) {
                design <- cbind(design, pcs)
            }
        }
    }

    rownames(design) <- observation.names
    
    list(dependent.variable=dependent.variable,
         independent.variable=independent.variable,
         family=family,
         matrix=design,
         sample.idx=sample.idx)
}


