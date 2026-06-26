# ewaff

An R package for running efficient and flexible epigenome-wide association studies (EWAS) of DNA methylation.

## Key features

- **Association testing** (`ewaff.sites`): Fit regression models at each CpG site using GLM, robust linear models (RLM), limma, or Cox proportional hazards. Supports any variable type as the variable of interest, including continuous, binary, and categorical.
- **Confounder adjustment**: Automatically generate surrogate variables (SVA/SmartSVA) or principal components (PCA) from the methylation matrix to correct for unknown batch effects and confounders.
- **Outlier handling** (`ewaff.handle.outliers`): Remove or winsorize per-site methylation outliers prior to analysis.
- **Meta-analysis** (`ewaff.meta.sites`): Meta-analyse EWAS summary statistics across multiple cohorts using random-effects models (via `metafor`).
- **Reporting** (`ewaff.summary`, `ewaff.report`): Generate HTML reports including QQ plots, Manhattan plots, and per-CpG scatter plots.
- **Parallel processing**: All computationally intensive steps use `mclapply` for multi-core execution.

## Comparison with meffil

[meffil](https://github.com/perishky/meffil) and ewaff are complementary packages from the same group, targeting different stages of the DNA methylation analysis workflow.

| | meffil | ewaff |
|---|---|---|
| **Primary focus** | Raw data processing & normalization | Statistical association testing |
| **Starting point** | Raw IDAT files from Illumina arrays | A pre-normalized methylation matrix |
| **Normalization** | Functional (quantile) normalization with control-probe PCA | Not applicable |
| **Array support** | 450k, EPIC v1/v2, MSA | Array-agnostic (works with any beta matrix) |
| **QC** | Sample and probe QC with outlier detection and reporting | Per-site outlier handling (IQR or winsorize) |
| **Cell counts** | Estimates cell-type composition from IDAT data | Not included |
| **EWAS** | Basic linear regression via `meffil.ewas` | GLM, robust LM, limma, or Cox PH via `ewaff.sites` |
| **Confounder generation** | PCs from control probes during normalization | SVA, SmartSVA, or PCA at the EWAS step |
| **Meta-analysis** | Not included | Random-effects meta-analysis across cohorts |
| **Memory efficiency** | Supports GDS files (matrix never loaded into RAM) | Requires in-memory matrix |

In practice the two packages are often used together: meffil produces the normalized beta matrix, and ewaff then performs the association analysis. ewaff can also be used independently with any normalized methylation matrix produced by other tools (e.g. minfi, SeSAMe).

## Installation

```r
library(remotes)
install_github("perishky/ewaff")
```

## Example analysis

See the [tutorial](http://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html) for a worked example using a simulated, reproducible dataset that demonstrates outlier handling, EWAS with SVA-based confounder generation, and report generation.


