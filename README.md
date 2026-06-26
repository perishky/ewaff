# ewaff

An R package for running efficient and flexible epigenome-wide association studies (EWAS) of DNA methylation.

## Key features

- **Association testing** (`ewaff.sites`): Fit regression models at each CpG site using GLM, robust linear models (RLM), limma, or Cox proportional hazards. Supports any variable type as the variable of interest, including continuous, binary, and categorical.
- **Confounder adjustment**: Automatically generate surrogate variables (SVA/SmartSVA) or principal components (PCA) from the methylation matrix to correct for unknown batch effects and confounders.
- **Outlier handling** (`ewaff.handle.outliers`): Remove or winsorize per-site methylation outliers prior to analysis.
- **Meta-analysis** (`ewaff.meta.sites`): Meta-analyse EWAS summary statistics across multiple cohorts using random-effects models (via `metafor`).
- **Reporting** (`ewaff.summary`, `ewaff.report`): Generate HTML reports including QQ plots, Manhattan plots, and per-CpG scatter plots.
- **Parallel processing**: All computationally intensive steps use `mclapply` for multi-core execution.

## Installation

```r
library(remotes)
install_github("perishky/ewaff")
```

## Example analysis

See the [tutorial](http://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html) for a worked example using a simulated, reproducible dataset that demonstrates outlier handling, EWAS with SVA-based confounder generation, and report generation.


