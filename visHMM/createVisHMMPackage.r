rm(list = ls())
library('roxygen2')

#' ggHMM package for visualizations for non-homogeneous Hidden Markov Models. 
#'
#' \tabular{ll}{
#' Package: \tab ggHMM\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2013-10-28\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Static visualizations based on time series, components, and structure of non-homogenous Hidden Markov Models that depend on one covariate
#'
#' @name ggHMM
#' @aliases ggHMM
#' @docType package
#' @title Visualizations for non-homogeneous HMMs
#' @author Adrian Albert \email{adalbert -at- stanford.edu}
#' @keywords package
#' @include ggplot2
#' @include useful
#' @include grid
#' @include reshape2
#' @include timeDate
#' @include lubridate
#' @include RColorBrewer
roxygen()

files = c('acf_ggplot.r', 'plot_utils.r')
package.skeleton(name = 'ggHMM', code_files = files)
roxygenize("ggHMM",
           roxygen.dir="ggHMM",
           copy.package=FALSE,
           unlink.target=FALSE)