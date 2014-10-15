#' tsbugs.
#'
#' @name tsbugs
#' @docType package
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the tsbugs package")
  packageStartupMessage("Plesae cite as:")
  packageStartupMessage("Abel, G.J., Bijak, J., Forster, J.J., Raymer J., Smith P.W.F. and Wong J.S.T. (2013). Integrating uncertainty in time series population forecasts: An illustration using a simple projection model. Demographic Research. 43 (29) 187--1226")
}
