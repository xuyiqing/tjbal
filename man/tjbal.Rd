\name{tjbal}
\alias{tjbal}
\title{Trajectory Balancing}
\description{Implements mean balancing and kernel balancing algorithms with time-series cross-sectional data.}
\usage{tjbal(formula=NULL, data, Y, D, X = NULL, X.avg.time = NULL, 
  index, index, trim.npre = 0, Y.match.time = NULL, Y.match.npre = NULL, 
  demean = TRUE, estimator = "meanfirst", sigma=NULL, 
  print.baltable = TRUE, vce = "jackknife", conf.lvl = 0.95, 
  nsims = NULL, parallel = TRUE, cores = 4, seed = 1234)
       }
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted. The first variable on the right-hand-side is a dichotomous treatment indicator (D); the rest of the right-hand-side variables are time-invariant controls (X). If X are time-varying, 
  they will be averaged within each unit based on values specified by \code{X.avg.time}.}
\item{data}{a data frame (must be a balanced panel).}
\item{Y}{outcome.}
\item{D}{treatment.}
\item{X}{covariates. If a covariate is time-varying, it will be averaged based on \code{X.avg.time} before balancing.}
\item{X.avg.time}{a list of time periods over which covariates are being averaged. Ignored if the treatment starts at different times.}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{trim.npre}{a numeric value indicating the smallest number of pre-treatment periods for a treated unit to be preserved. The default is 0.}    
\item{Y.match.time}{a set of pre-treatment time periods in which the outcome variable is being balanced on.}
\item{Y.match.npre}{a numeric value indicating the number of pre-treatment outcome periods to be balanced on. If \code{Y.match.npre = 0}, no pre-treatment outcome will be part of the balancing scheme.}
\item{demean}{a logical flag indicating whether a demeaning procedure will be performed to take out the average of pre-treatment outcomes for each unit.}
\item{estimator}{a string specifying the balancing approach: "mean" for mean-balancing, "kernel" for kernel-balancing, and "meanfirst" (default) for kernel balancing with mean balancing constraints. "meanfirst" will prioritize balancing on covariate means over higher-order terms and interactions.}
\item{sigma}{a numeric value specifying the bandwidth of the kernel transformation (will be multiplied by the number of covariates to be balanced on). The default is 2. Ignored if \code{estimator = "mean"}.}
\item{print.baltable}{a logical flag that controls whether to print out the balance table after the algorithm is run. Ignored if treatment timing is different.}
\item{vce}{a string specifying the variance estimator. \code{vce = "none"}: no uncertainty estimates; \code{vce = "fixed.weights"}: treating balancing weights as fixed; \code{vce = "bootstrap"}: conducting non-parametric bootstrapping by reshuffle both the treated and control units; \code{vce = "jackknife"}: conducting jackknife by omitting one treated unit at a time. When the treatment timing is different, \code{vce = "jackknife"} is the only available choice for producing uncertainty estimates.}
\item{conf.lvl}{a positive number in the range of 0 and 1 specifying confidence levels for uncertainty estimates. The default value is 0.95.}
\item{nsims}{an integer specifying the number of bootstrap or jackknife runs. Ignored if \code{vce = "none"} or there are too few treated units.}
\item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrap/jackknife simulations. .}
\item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
\item{seed}{an integer that sets the seed in random number generation. Ignored if \code{vce = "none"} or there are too few treated units.}
}
\details{
  \code{tjbal} provides a general reweighting approach to causal inference with time-series cross-sectional (TSCS) data. It includes two estimators, mean balancing and kernel balancing. The former reweights control units such that the averages of the pre-treatment outcomes and covariates are approximately equal between the treatment and (reweighted) control groups. The latter relaxes the linearity assumption and seeks approximate balance on a kernel-based feature expansion of the pre-treatment outcomes and covariates. The resulting approach inherits the ability of synthetic control and latent factor models to tolerate time-varying confounders, but (1) improves feasibility and stability with reduced user discretion; (2) accommodates both short and long pre-treatment time periods with many or few treated units; and (3) balances on the high-order ``trajectory" of pre-treatment outcomes rather than their period-wise average. We illustrate this method with simulations and two empirical examples. 
}
\value{
  \item{data.wide}{a matrix storing data in wide form -- each row represents one unit.}
  \item{id.tr}{a vector of row numbers for the treatment units.}
  \item{id.co}{a vector of row numbers for the control units.}
  \item{Y.tr}{data of the treated unit outcome.}
  \item{Y.co}{data of the control unit outcome.}
  \item{Y.var}{a vector of outcome variable names.}
  \item{matchvar.list}{a list of covariates to be balanced on for each subgroup.}
  \item{Ttot}{the set of all time periods.}
  \item{N}{the total number of units.}
  \item{Ntr}{the number of treated units.}
  \item{Nco}{the number of control units.}
  \item{T0}{a vector that stores the unique values of the timing of the treatment.}
  \item{T0.all}{a vector that stores the timing of the treatment for each unit.}
  \item{T0.tr}{a vector that stores the timing of the treatment for each treated unit.}
  \item{weights.co}{weights of the control unit; added up to 1.}
  \item{Y.bar}{average values of treated outcomes, counterfactual outcomes, and control outcomes over time.}
  \item{att}{average treatment effect on the treated over time (it is realigned and based on the timing of the treatment if it is different for each unit).}
  \item{att.avg}{average treatment effect on the treated (averaged both across unit and over time).}  
  \item{est.att}{inference for \code{att}.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{ntreated}{a vector of numbers of treated units for all subgroups.}
  \item{sameT0}{TRUE if the timing of the treatment is the same.}  
  \item{Ttot}{the set of all time periods.}  
  \item{Tpre}{the set of pre-treatment periods; same timing.}
  \item{Tpst}{the set of post-treatment periods; same timing.}
  \item{bal.type}{the type of balancing scheme being used: "mbal" for mean balancing and "kbal" for kernel balancing; same timing.}
  \item{ndims}{the number of dimensions being balanced; same timing.}
  \item{b}{tunning parameter for kbal; same timing.}
  \item{kbal.out}{output from kbal; same timing.}
  \item{success}{a logical flag indicating whether convergence is achieve with balancing.}
  \item{bias.ratio}{the ratio of L1 measure of distance after balancing over L1 measure of distance before balancing. A smaller number indicates more improvement brought by balancing.}
  \item{bal.table}{balance table; same timing.}
  \item{att.sims}{jackknife/bootstrap results for \code{att}; same timing.}
  \item{att.avg.sims}{jackknife/bootstrap results for \code{att.avg}; same timing.}
  \item{att.sub.att}{jackknife results for \code{att} for all subgroups; different timing.}
  \item{att.sub.att.avg}{jackknife results for \code{att.avg} for all subgroups; different timing.}
  \item{sub.weights.co}{a matrix of weights for the control units for all subgroups; different timing.}
  \item{sub.Ytr.avg}{a matrix of average outcomes for the treated units for all subgroups; different timing.}
  \item{sub.Yct.avg}{a matrix of average Y(0)'s for the treated units for all subgroups; different timing.}
  \item{sub.att}{a matrix of ATT for all subgroups; different timing.}
  \item{sub.ntr}{a matrix of numbers of treated units over time for all subgroups; different timing.}
  \item{sub.att.adj}{a matrix of realigned ATT for all subgroups; different timing.}
  \item{group.stats}{statistics for all subgroups, including T0, the number of treated units, sigma, convenience, and bias ratio; different timing.}
  \item{bal.table.list}{a list of the balance table for all subgroups; different timing.}
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu (Maintainer) <yiqingxu@stanford.edu>, Stanford
}
\references{
  Hazlett, Chad and Yiqing Xu, 2020. ``Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data.'' Working Paper, UCLA and UCSD.
}
\seealso{
  \code{\link{plot.tjbal}} and \code{\link{print.tjbal}}
}
\examples{
library(tjbal)
data(tjbal)
out <- tjbal(roa ~ treat + so_portion + rev2007, data = npc, 
       index = c("gvkey","fyear"), estimator = "mean") 
print(out) 
}


