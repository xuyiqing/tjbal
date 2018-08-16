\name{tjbal}
\alias{tjbal}
\title{Trajectory Balancing}
\description{Implements mean balancing and kernel balancing with time-series cross-sectional data.}
\usage{tjbal(formula=NULL, data, Y, D, X = NULL, X.avg.time = NULL,
       index, Y.match.periods = NULL, demean = TRUE, kernel = FALSE, 
       sigma=NULL, maxnumdims = NULL, method="ebal", 
       whiten=FALSE, test = FALSE, nsigma = 16, print.baltable = TRUE, 
       bootstrap = FALSE, conf.lvl = 0.95, nboots = 500, parallel = TRUE, cores = 4, seed = 1234)
       }
\arguments{
\item{formula}{an object of class "formula": a symbolic description of
  the model to be fitted. The first variable on the right-hand-side is a dichotomous treatment indicator (D); 
  the rest of right-hand-side variables are time-invariant controls (X). If X are time-varying, 
  they will be averaged within each unit based on values specified by \code{X.avg.time}.}
\item{data}{a data frame (must be a balanced panel).}
\item{Y}{outcome.}
\item{D}{treatment.}
\item{X}{time-invariant covariates.}
\item{X.avg.time}{}
\item{index}{a two-element string vector specifying the unit (group)
    and time indicators. Must be of length 2.}
\item{Y.match.periods}{}
\item{demean}{a logical flag indicating whether a demeaning procedure will be performed to take out the average of pre-treatment outcomes for each unit.}
\item{kernel}{a logical flag indicating whether the kernel transformation will be used; if \code{kernel = FALSE}, the mean balancing estimator will be used.}
\item{sigma}{a numeric value specifying the bandwidth of the kernel transformation. Ignored if \code{kernel = FALSE}.}
\item{maxnumdims}{a positive integer specifying the maximum number of dimensions to be considered in kernel balancing. See the manual of \code{KBAL} for more details.}
\item{method}{a string that specify the objective function of the balancing algorithm, either "ebal" (entropy balancing) or "el" (empirical likelihood).}
\item{whiten}{a logical flag controlling pre-whitening of the data prior to construction of K. See the manual of \code{KBAL} for more details.}
\item{test}{a logical flag indicating whether cross-validation will be performed to choose sigma, the tuning parameter.}
\item{nsigma}{a positive integer specifying the number of candidate sigmas to be considered in cross-validation.}
\item{print.baltable}{a logical flag that controls whether to print out the balance table after the algorithm is run.}
\item{bootstrap}{a logical flag indicating whether uncertainty estimates will be produced using bootstrap.}
\item{conf.lvl}{a positive number in the range of 0 and 1 specifying confidence levels for uncertainty estimates. The default value is 0.95.}
\item{nboots}{an integer specifying the number of bootstrap runs. Ignored if \code{bootstrap = FALSE}.}
\item{parallel}{a logical flag indicating whether parallel computing will be used in bootstrapping. Ignored if \code{bootstrap = FALSE}.}
\item{cores}{an integer indicating the number of cores to be used in parallel computing. If not specified, the algorithm will use the
  maximum number of logical cores of your computer (warning: this could prevent you from multi-tasking on your computer).}
\item{seed}{an integer that sets the seed in random number generation. Ignored if \code{bootstrap = FALSE}.}
}
\details{
  \code{tjbal} provides a general reweighting approach to causal inference with time-series cross-sectional (TSCS) data. It includes two estimators, mean balancing and kernel balancing. The former reweights control units such that the averages of the pre-treatment outcomes and covariates are approximately equal between the treatment and (reweighted) control groups. The latter relaxes the linearity assumption and seeks approximate balance on a kernel-based feature expansion of the pre-treatment outcomes and covariates. The resulting approach inherits the ability of synthetic control and latent factor models to tolerate time-varying confounders, but (1) improves feasibility and stability with reduced user discretion; (2) accommodates both short and long pre-treatment time periods with many or few treated units; and (3) balances on the high-order ``trajectory" of pre-treatment outcomes rather than their period-wise average. We illustrate this method with simulations and two empirical examples. 
}
\value{
  \item{DID}{TRUE if the timing of the treatment is the same.}
  \item{data.wide}{a matrix storing data in wide form -- each row represents one observation.}
  \item{id.tr}{a vector of row numbers for the treatment units.}
  \item{id.co}{a vector of row numbers for the control units.}
  \item{Y.tr}{data of the treated unit outcome.}
  \item{Y.co}{data of the control unit outcome.}
  \item{weights.co}{weights of the control unit; added up to 1.}
  \item{names.co}{IDs of the control units}
  \item{Ttot}{the set of all time periods.}
  \item{Tpre}{the set of pre-treatment periods.}
  \item{Tpre}{the set of post-treatment periods.}
  \item{T0}{a vector that stores the timing of the treatment.}
  \item{N}{the total number of units.}
  \item{Ntr}{the number of treated units.}
  \item{Nco}{the number of control units.}
  \item{matchvar}{the set of pre-treatment variables being matched on.}
  \item{Y.bar}{average values of treated outcomes, counterfactual outcomes, and control outcomes over time.}
  \item{att}{average treatment effect on the treated over time (it is averaged 
    based on the timing of the treatment if it is different for each unit).}
  \item{att.avg}{average treatment effect on the treated (averaged both across unit and over time.}
  \item{success}{a logical flag indicating whether significant improvement in covariate distance is made through balancing.}
  \item{L1.before}{L1 measure of distance prior to balancing.}
  \item{L1.after}{L1 measure of distance after balancing.}
  \item{L1.ratio}{L1.after divided by L1.before -- a smaller number means more improvement is made.}
  \item{min90}{the minimum number of dimensions that explains over 90\% of the variation in the kernel matrix.}
  \item{ndims}{the number of dimensions of the kernel matrix being balanced on.}
  \item{sigma.best}{the turning parameter being used in kernel balancing.}
  \item{bal.table}{the balance table}
  \item{est.att}{inference for \code{att}.}
  \item{est.att.avg}{inference for \code{att.avg}.}
  \item{att.boot}{bootstrap results for \code{att}.}
  \item{att.avg.boot}{bootstrap results for \code{att.avg}.}
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu (Maintainer) <yiqingxu@ucsd.edu>, UCSD
}
\references{
  Hazlett, Chad and Yiqing Xu, 2018. ``Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data.'' Working Paper, UCLA and UCSD.
}
\seealso{
  \code{\link{plot.tjbal}} and \code{\link{print.tjbal}}
}
\examples{
library(tjbal)
data(tjbal)
out <- tjbal(roa ~ treat + so_portion + rev2007, data = npc, 
       index = c("gvkey","fyear"), kernel = FALSE) 
print(out) 
}


