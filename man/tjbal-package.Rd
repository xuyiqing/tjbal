\name{tjbal-package}
\alias{tjbal-package}
\docType{package}
\title{Trajectory Balancing}
\description{Implements mean balancing and kernel balancing with time-series cross-sectional data.}
\details{
  Provides a general reweighting approach to causal inference with time-series cross-sectional (TSCS) data. It includes two estimators, mean balancing and kernel balancing. The former reweights control units such that the averages of the pre-treatment outcomes and covariates are approximately equal between the treatment and (reweighted) control groups. The latter relaxes the linearity assumption and seeks approximate balance on a kernel-based feature expansion of the pre-treatment outcomes and covariates. 

  See \code{\link{tjbal}} for details.
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu (Maintainer) <yiqingxu@ucsd.edu>, UCSD
}
\references{
  Hazlett, Chad and Yiqing Xu, 2018. ``Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data.'' Working Paper, UCLA and UCSD.
}
