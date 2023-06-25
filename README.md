
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tjbal

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**R** package for implementing trajectory balancing, a kernel-based
reweighting method for causal inference with panel data.

**Repo:** [GitHub](https://github.com/xuyiqing/tjbal)

**Examples:** R code used in the
[tutorial](https://yiqingxu.org/packages/tjbal/articles/tutorial.html)
can be downloaded from [here](tjbal_examples.R).

**Reference**: Hazlett, Chad and Yiqing Xu, 2018. “Trajectory Balancing:
A General Reweighting Approach to Causal Inference with Time-Series
Cross-Sectional Data.” Working Paper, UCLA and Stanford. Available at
SSRN: <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3214231>.

## Installation

<!---
You can install **fect** directly from CRAN by typing the following command in the **R** console: 
&#10;
```r
install.packages('fect', type = 'source')
```
--->

You can install the development version of the package from Github by
typing the following commands:

``` r
install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
# devtools::install_github('chadhazlett/kbal')
devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel")
devtools::install_github('xuyiqing/tjbal')
```

Note that installing **kbal** (from Github) is required. **tjbal** also
depends on the following packages, which will be installed automatically
when **tjbal** is installed. You can also install them manually:

``` r
## for plotting
require(ggplot2)  
## for parallel computing 
require(foreach)  
require(doParallel) 
require(parallel)
## for data manipulation 
require(plyr)
```

**panelView** for panel data visualization is also highly recommended:

``` r
devtools::install_github('xuyiqing/panelView')
```

### Notes on installation failures

1.  Mac users who have updated to MacOS BigSur or Monterey will likely
    encounter compilation problems. See
    [here](http://yiqingxu.org/public/BigSurError.pdf) for a potential
    solution.
2.  Windows users please consider upgrading R to 4.0.0 or higher and
    installing the [latest
    Rtools](https://cran.r-project.org/bin/windows/Rtools/) to avoid
    C++17 complier errors when installing fastplm.
3.  For Rcpp, RcppArmadillo and MacOS “-lgfortran” and “-lquadmath”
    error, click
    [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)
    for details.
4.  Installation failure related to OpenMP on MacOS, click
    [here](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/)
    for a solution.
5.  To fix these issues, try installing gfortran from
    [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS%20clang4%20R%20Binaries%20from%20https://github.com/coatless/r-macos-clang).

## Report bugs

Please report bugs to **yiqingxu \[at\] stanford.edu** with your sample
code and data file. Much appreciated!
