extract.survival.variables <- function(formula) {
    if (class(formula) == "formula")
        formula <- as.character(formula)[2]
    strsplit(sub("Surv\\(([^)]+)\\).*", "\\1", formula), "[ ]*,[ ]*")[[1]]
}

is.survival.expression <- function(formula) {
    if (class(formula) == "formula")
        formula <- as.character(formula)[2]
    length(grep("Surv\\(.+\\)", formula)) > 0
}
