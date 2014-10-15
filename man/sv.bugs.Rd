\name{sv.bugs}
\alias{sv.bugs}
\title{
Create BUGS Script of a Stochastic Volatility (SV) Model
}
\description{
Create BUGS script of an SV time series model, similar those in Meyer and Yu (2002). Options allow for the inclusion of a different lag orders for the mean and variance terms, forecasts, posterior simulations from the model and alternative specification of prior distributions on each parameters.
}
\usage{
sv.bugs(y, ar.order = 0, k = NULL, sim = FALSE, mean.centre = FALSE, beg = ar.order + 1, 
		mean.prior = ar.prior, ar.prior = "dnorm(0,1)", 
		sv.order = 1, sv.beg = beg + sv.order,
		sv.mean.prior1 = "dnorm(0,0.001)", sv.mean.prior2 = NULL, 
		sv.ar.prior1 = "dunif(0,1)", sv.ar.prior2 = NULL, 
		sv.tol.prior = "dgamma(0.01,0.01)", 
space = FALSE)
}
\arguments{
  \item{y}{
Data to be used for the BUGS model.
}
  \item{ar.order}{
AR order of the mean process for BUGS model.
}
  \item{k}{
Length of forecast horizon to be included in the BUGS model.
}
  \item{sim}{
Enable posterior simulations to be included in the BUGS model. Default is \code{FALSE}.
}
  \item{mean.centre}{
Include a term to centre the data on its mean value. Default is \code{FALSE}.
}
  \item{beg}{
Starting value for which data are considered onwards (and including) in the likelihood of the BUGS model. By default this is the \code{ar.order+1} but if comparing models of different orders, users may wish to set all \code{beg} to the same value.
}
  \item{mean.prior}{
Prior for mean term (not used if mean.centre is not set to \code{TRUE}). The distribution should be stated in BUGS syntax. By default, the same prior as the autoregressive terms are used.
}
  \item{ar.prior}{
Prior for autoregressive terms. The distribution should be stated in BUGS syntax. By default set to a normal distribution with mean 0 and tolerance 1 (\code{dnorm(0,1)}). The same prior is used for all autoregressive terms.
}
  \item{sv.order}{
AR order of the volatility process for BUGS model.
}
  \item{sv.beg}{
Starting value for which the volatility process in the data are considered onwards (and including) in the likelihood of the BUGS model. Default is \code{beg + sv.order} but if comparing models of different orders of the volatility process, users may wish to set all \code{beg} to the same value.
}
  \item{sv.mean.prior1}{
Prior distribution for the mean volatility term. The distribution should be stated in BUGS syntax. By default set to a normal distribution with mean 0 and tolerance 0.001 (\code{dnorm(0,0.001)}), as in Meyer and Yu (2002).
}
  \item{sv.mean.prior2}{
Alternative prior for mean volatility term. The distribution set here will be transformed by taking the negative logarithm. The distribution should be stated in BUGS syntax.
}
  \item{sv.ar.prior1}{
Prior for autoregressive terms of the volatility process. The distribution set here will be transformed by doubling and then subtracting 1 as in Meyer and Yu (2002). The distribution syntax must be recognisable to BUGS. By default set to a uniform distribution with bounds 0 and 1 (\code{dunif(0,1)}). 
}
  \item{sv.ar.prior2}{
Alternative prior for autoregressive terms of the volatility process. The distribution should be stated in BUGS syntax. 
}
  \item{sv.tol.prior}{
Prior for the tolerance of the volatility. The distribution should be stated in BUGS syntax. By default set to a uninformative gamma distribution.
}
  \item{space}{
Include some additional empty lines to seperate the likelihood, priors, forcasts and simulation components of the BUGS model.
}
}
\details{
This function create BUGS scripts of an SV time series model, similar those in Meyer and Yu (2002). Prior distributions should be set up using BUGS syntax. For example, \code{dnorm} is a normal distribution with mean and tolerance (not variance) arguments. Prior distributions for the mean and autoregressive parameter in the volatility process can be specified directly in \code{sv.mean.prior1} and \code{sv.ar.prior1} respectively. In \code{sv.mean.prior2} the prior is given on the mean variance, and then transformed to the mean volatility. In \code{sv.ar.prior2} the prior is given on the positive scale as illustrated in Meyer and Yu (2002). Only one set of prior distributions are allowed (i.e. either \code{sv.ar.prior1} or \code{sv.ar.prior2} should specified by the user). 

The data \code{y}, can contain missing values. Note, if missing values are close the beginning of the series when a high order model for the mean process is specified (i.e. the second data point is missing and a AR(4) is specified) the user with have to set a high starting point for model to be fitted on (\code{beg}) for the BUGS model to function (i.e. 7).
}
\value{
\item{bug }{A BUGS model of type \code{tsbugs}.}
\item{data }{The data to be used with the model. This might extend the original data passed to the function with \code{k} unknown future values to be forecast.}
\item{info }{Additional information on the length of the data, variance type and line numbers of certain parts of the BUGS model.}
}
\references{
Meyer, R. and J. Yu (2002). BUGS for a Bayesian analysis of stochastic volatility models. \emph{Econometrics Journal} 3 (2), 198--215.
}
\author{
Guy J. Abel
}
\seealso{
\code{\link{ar.bugs}}, \code{\link{rv.bugs}}, \code{\link{garch.bugs}}
}
\examples{
# Create AR(0)-SV model for svpdx
y <- svpdx$pdx
plot(y, type = "l")
sv0 <- sv.bugs(y)
print(sv0)

# AR(1)-SV model with AR(2) structure in the volatility
sv1 <- sv.bugs(y, ar.order = 1, sv.order = 2)
print(sv1)

# AR(0)-SV model with alternative prior
sv0 <- sv.bugs(y, sv.ar.prior2 = "dunif(-1,1)")
print(sv0)

# AR(0)-SV model with forecast and posterior simulations
sv0 <- sv.bugs(y, k = 10, sim = TRUE)
print(sv0) 

\dontrun{
  # Run in OpenBUGS
  writeLines(sv0$bug, "sv0.txt")
  library("R2OpenBUGS")
  
  sv0.bug <- bugs(data = sv0$data,
  				  inits = list(inits(sv0)),
  				  param = c(nodes(sv0, "prior")$name, "y.sim", "y.new"), 
  				  model = "sv0.txt", 
  				  n.iter = 20000, n.burnin = 10000, n.chains = 1)
  
  # Plot the parameters posteriors and traces
  library("coda")
  param.mcmc <- as.mcmc(sv0.bug$sims.matrix[, nodes(sv0, "prior")$name])
  plot(param.mcmc)
  
  # Plot posterior simulations using fanplot
  library("fanplot")
  y.mcmc <- sv0.bug$sims.list$y.sim
  y.pn <- pn(y.mcmc, st = 1)
  plot(y, type = "n")
  fan(y.pn)
  lines(y)
  
  # Plot forecast using fanplot
  ynew.mcmc <- sv0.bug$sims.list$y.new
  ynew.pn <- pn(ynew.mcmc, st = sv0$info$n + 1)
  plot(y, type = "n", xlim = sv0$info$n + c(-100, 20))
  fan(ynew.pn)
  lines(y)
  
  # Plot volatility
  h.mcmc <- sv0.bug$sims.list$h
  h.pn <- pn(h.mcmc, st = 1)
  sigma.pn <- pn(sims = sqrt(exp(h.mcmc)), st = 1)
  par(mfrow = c(2, 1), mar = rep(2, 4))
  plot(NULL, type = "n", xlim = tsp(h.pn)[1:2], ylim = range(h.pn[, 5:95]), main = "h_t")
  fan(h.pn)
  abline(v = length(y))
  plot(NULL, type = "n", xlim = tsp(sigma.pn)[1:2], ylim = range(sigma.pn[, 1:95]), main = "sigma_t")
  fan(sigma.pn)
  abline(v = length(y)) 
}
}

