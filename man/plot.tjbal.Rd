\name{plot.tjbal}
\alias{plot.tjbal}
\title{Plotting}
\description{Visualizes results from trajectory balancing.}
\usage{\method{plot}{tjbal}(x, type = "gap", xlim = NULL, ylim = NULL,
    xlab = NULL, ylab = NULL, count = TRUE, legendOff = FALSE, 
    legend.pos = "bottom", legend.ncol = NULL, legend.labs = NULL, 
    main = NULL, raw = "none", stat = "mean", trim = TRUE, trim.wtot = 0.9, 
    theme.bw = TRUE, cex.main = NULL, cex.axis = NULL, cex.lab = NULL, 
    cex.legend = NULL, cex.text = NULL, log.weights = FALSE, wmin = NULL, 
    axis.adjust = FALSE, \dots)
}
\arguments{
  \item{x}{a \code{\link{tjbal}} object.}
  \item{type}{a string that specifies the type of the plot. Must be one
    of the following: "gap" (plotting the average treatment effect on
    the treated; "counterfactual" or "ct" for short (plotting predicted average Y(0)'s); "balance" (plotting balance of covariates and pre-treatment outcomes before and after reweighting); "weights" (plotting the distribution of weights of the control units -- weights added up to the number of the treated units).}
  \item{subgroup}{a numeric value indicating the subgroup for which the results will be plotted.}  
  \item{xlim}{a two-element numeric vector specifying the range of x-axis. When 
    class of time varible is string, must specify not original value but a counting number e.g. \code{xlim=c(1,30)}.}
  \item{ylim}{a two-element numeric vector specifying the range of y-axis.}
  \item{xlab}{a string indicating the label of the x-axis.}
  \item{ylab}{a string indicating the label of the y-axis.}
  \item{count}{a logical flag indicating whether a bar plot of the number of treated units will be shown at the bottom of the plot. The default is TRUE. Ignored if \code{type = "balance"} or \code{type = "weights"}.}
  \item{legendOff}{a logical flag controlling whether to show the legends.}
  \item{legend.pos}{a character controlling the position of the legends. The default is "bottom".}
  \item{legend.ncol}{a numeric value controlling the number of columns of the legends.}
  \item{legend.labs}{a vector of characters specifying the labels of the legends.}
  \item{main}{a string that controls the title of the plot. If not
    supplied, no title will be shown.}
  \item{raw}{a string indicating whether or how raw data for the outcome
    variable will be shown in the "counterfactual" plot. Ignored if \code{type} is not "counterfactual" or treatment timing is different. Must be one of the following: "none" (not showing the raw data); "band" (showing the middle 90 percentiles of the raw data); and "all" (showing the raw data as they are).}
  \item{stat}{a string that specifies the statistics in the \code{balance} plot, either "mean" (showing the difference in means of covariates between the treatment and control groups before and after reweighting, normalized by the standard deviation of the corresponding covariates among the treated units, default) or "sd" (showing the difference in standard deviations of covariates between the treatment and control groups before and after reweighting, normalized by the standard deviation of the corresponding covariates among the treated units). Ignored if \code{type} is not "balance".}  
  \item{trim}{a logical flag controlling whether to trim the control group data in a counterfactual plot when \code{raw = "all"} or \code{raw = "band"}. The default is TRUE. }
  \item{trim.wtot}{a numeric value controlling the total sum of weights of the most heavily weighted controls to remain in the plot. The default is 0.9.}
  \item{theme.bw}{a logical flag controlling whether a black-white theme will be used. The default is TRUE.}  
  \item{cex.main}{a numeric value controlling the font size of the title.}
  \item{cex.axis}{a numeric value controlling the font size of the axis labels.}  
  \item{cex.lab}{a numeric value controlling the font size of the axis titles.}  
  \item{cex.legend}{a numeric value controlling the font size of the legends.}  
  \item{cex.text}{a numeric value controlling the font size of the text in the figure (such as Ntr = 47).}  
  \item{log.weights}{a logical flag controlling whether to take a common logarithmic transformation of the weights. The default is FALSE. Ignored if \code{type} is not"weights".}  
  \item{wmin}{a numeric number specifying the minimum weight (in common logarithmic). Ignored if \code{type} is not"weights".}
  \item{axis.adjust}{a logic flag indicating whether to adjust labels on x-axis. 
    Useful when class of time variable is string and data magnitude is large.}
  \item{\dots}{other argv.}  
}
\details{
  \code{plot.tjbal} visualizes the results of trajectory balancing.
}
\author{
  Chad Hazlett <chazlett@ucla.edu>, UCLA

  Yiqing Xu (Maintainer) <yiqingxu@stanford.edu>, Stanford
}
\references{
  Hazlett, Chad and Yiqing Xu, 2020. ``Trajectory Balancing: A General Reweighting Approach to Causal Inference with Time-Series Cross-Sectional Data.'' Working Paper, UCLA and Stanford.
}
\seealso{
  \code{\link{tjbal}} and \code{\link{print.tjbal}}
}



