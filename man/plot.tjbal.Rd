\name{plot.tjbal}
\alias{plot.tjbal}
\title{Plotting}
\description{Visualizes estimation results of trajectory balancing.}
\usage{\method{plot}{tjbal}(x, type = "gap", xlim = NULL, ylim = NULL,
            xlab = NULL, ylab = NULL, legendOff = FALSE, raw = "none",
            main = NULL, nfactors = NULL, id = NULL, axis.adjust = FALSE, \dots)
}
\arguments{
  \item{x}{a \code{\link{tjbal}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the average treatment effect on
    the treated; "counterfactual" or "ct" for short (plotting predicted Y(0)'s); 
    "balance" (plotting balance of covariates and pre-treatment outcomes before and after reweighting); 
    "weights" (plotting the distribution of weights of the control units -- weights added up to the number of the treated units).}
    \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{legendOff}{a logical flag controlling whether to show the legend.}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{raw}{a string indicating whether or how raw data for the outcome
    variable will be shown in the "counterfactual" plot. Ignored if \code{type} is not
    "counterfactual". Must be one of the following: "none" (not showing
    the raw data); "band" (showing the middle 90 percentiles of the raw
    data); and "all" (showing the raw data as they are).}
  \item{wmin}{a numeric number specifying the minimum weight (in common logarithmic). Ignored if \code{type} is not"weights".}
  \item{stat}{a string that specifies the statistics in the \code{balance} plot, either "mean" (showing the difference in means of covariates between the treatment and control groups before and after reweighting, normalized by the standard deviation of the corresponding covariates among the treated units, default) or "sd" (showing the difference in standard deviations of covariates between the treatment and control groups before and after reweighting, normalized by the standard deviation of the corresponding covariates among the treated units). Ignored if \code{type} is not "balance".}  
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis. 
    Useful when class of time variable is string and data magnitude is large.}
  \item{\dots}{other argv.}  
}
\details{
  \code{plot.tjbal} visualizes the results of trajectory balancing.
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu (Maintainer) <yiqingxu@ucsd.edu>, UCSD
}
\references{
  Hazlett, Chad and Yiqing Xu, 2018. ``Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data.'' Working Paper, UCLA and UCSD.
}
\seealso{
  \code{\link{tjbal}} and \code{\link{print.tjbal}}
}



