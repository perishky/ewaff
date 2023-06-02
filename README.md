# ewaff

Efficient and flexible algorithms for EWAS

## Installation 

```
library(remotes)
install_github("perishky/ewaff")
```

## Example analysis

[tutorial](http://htmlpreview.github.io/?https://github.com/perishky/ewaff/blob/master/docs/tutorial.html).

## Features to add

- multi-core optional for windows or other debilitated users
- use tableone for output (per model) in EWAS object and report
- EWAS comparison function
- expand EWAS report for multi-variable ewas
- support adjustment for case-control
- flexible batch control
- correlation threshold for association between variable of interest and SVAs
- **a function called ewaff.pace() that runs the PACE workflow**
- Another way to potentially speed up SVA: https://cran.r-project.org/package=irlba
