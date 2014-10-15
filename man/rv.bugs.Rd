\name{rv.bugs}
\alias{rv.bugs}
\title{
Create BUGS Script of a Random Variance Shift Model
}
\description{
Create BUGS script of an Random Variance (RV) shift model similar to that of McCulloch and Tsay (1991). Options allow for the inclusion of a different lag orders for the mean term, forecasts, posterior simulations from the model and alternative specification of prior distributions on some parameters.
}
\usage{
rv.bugs(y, ar.order = 0, k = NULL, sim = FALSE, mean.centre = FALSE, beg = ar.order + 1, 
		mean.prior = ar.prior, ar.prior = "dnorm(0,1)", 
		rv.tol0.prior = "dgamma(0.000001,0.000001)", rv.eps.prior = "dbeta(1, 100)", 
		rv.ilambda2.prior = "dgamma(0.01,0.01)",
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
Starting value for which data are considered onwards (and including) in the likelihood of the BUGS model. Default is \code{ar.order+1} but if comparing models of different orders, users may wish to set all \code{beg} to the same value.
}
  \item{mean.prior}{
Prior for mean term (not used if mean.centre is not set to \code{TRUE}). The distribution should be stated in BUGS syntax. By default, the same prior as the autoregressive terms are used.
}
  \item{ar.prior}{
Prior for autoregressive terms. The distribution should be stated in BUGS syntax. By default this is set to a normal distribution with mean 0 and tolerance 1 (\code{dnorm(0,1)}). The same prior is used for all autoregressive terms.
}
  \item{rv.tol0.prior}{
Prior distribution for the time-specific tolerance of the first data point, from which potential future shifts are based on. This must be a distribution of syntax recognisable to BUGS. By default this is set to a uninformative gamma distribution.
}
  \item{rv.eps.prior}{
Prior distribution for the epsilon term (the probability of the variance shift). This must be a distribution of syntax recognisable to BUGS and in theory restricted to generate values between 0 and 1, although no check is made to ensure users specify distributions as such. By default this argument is set to a beta distribution with a small probability for a variance shift (\code{dbeta(1,100)}).
}
  \item{rv.ilambda2.prior}{
Prior for the inverse of the squared lambda term of a variance shift model. The lambda represents the magnitude (on the log scale) of average variance shifts, thus \code{ilambda2} in a BUGS model is a similar to a tolerance on the shifts in the variance, and hence must be positive. By default this is set to a gamma distribution (\code{dgamma(0.01,0.01)}).
}
  \item{space}{
Include some additional empty lines to seperate the likelihood, priors, forcasts and simulation components of the BUGS model.
}
}
\details{
This function create BUGS scripts of an random variance shift model adapted from McCulloch and Tsay (1991). Prior distributions for the tolerance of the initial data point, the probability of the variance shift, and the magnitude of average variance shifts can be specified by users to differ from their default values. User specified prior distributions should be set up using BUGS syntax. For example, \code{dnorm} is a normal distribution with mean and tolerance (not variance) arguments.

The data \code{y}, can contain missing values. Note, if missing values are close the beginning of the series when a high order model for the mean process is specified (i.e. the second data point is missing and a AR(4) is specified) the user with have to set a high starting point for model to be fitted on (\code{beg}) for the BUGS model to function (i.e. 7).
}
\value{
\item{bug }{A BUGS model of type \code{tsbugs}.}
\item{data }{The data to be used with the model. This might extend the original data passed to the function with \code{k} unknown future values to be forecast.}
\item{info }{Additional information on the length of the data, variance type and line numbers of certain parts of the BUGS model.}
}
\references{
McCulloch, R. E. and R. S. Tsay (1993). Bayesian Inference and Prediction for Mean and Variance Shifts in Autoregressive Time Series. \emph{Journal of the American Statistical Association} 88 (423), 968--978.
}
\author{
Guy J. Abel
}
\seealso{
\code{\link{ar.bugs}}, \code{\link{sv.bugs}}, \code{\link{garch.bugs}}
}
\examples{
# Create AR(0)-SV model for population growth rate
r <- ts(ew[2:167]/ew[1:166]-1, start=1841)
y <- diff(r) 
plot(y, type = "l")
rv0 <- rv.bugs(y)
print(rv0)

# AR(2)-RV model with alternative priors
rv2 <- rv.bugs(y, rv.eps.prior = "dbeta(1,20)")
print(rv2)

# AR(0)-RV model with posterior simulations
rv0 <- rv.bugs(y, sim = TRUE)
print(rv0) 

\dontrun{
  # Run in OpenBUGS
  writeLines(rv0$bug, "rv0.txt")
  library("R2OpenBUGS")
  
  # Run model (can take some time, depending on data length)
  rv0.bug <- bugs(data = rv0$data,
  				 inits = list(inits(rv0)), 
  				 param = c(nodes(rv0, "prior")$name, "y.sim", "y.new"), 
  				 model = "rv0.txt", 
  				 n.iter = 11000, n.burnin = 1000, n.chains = 1)
  
  # Plot the parameters posteriors and traces
  library("coda")
  param.mcmc <- as.mcmc(rv0.bug$sims.matrix[, nodes(rv0, "prior")$name])
  plot(param.mcmc)
  
  # Plot posterior simulations using fanplot
  library("fanplot")
  y.mcmc <- rv0.bug$sims.list$y.sim
  y.pn <- pn(y.mcmc, st = start(y)[1])
  plot(y, type = "n")
  fan(y.pn)
  lines(y)
  
  # Plot volatility
  h.mcmc <- rv0.bug$sims.list$h
  h.pn <- pn(h.mcmc, st = start(y)[1])
  sigma.pn <- pn(sims = sqrt(exp(h.mcmc)), st = start(y)[1])
  par(mfrow = c(2, 1), mar = rep(2, 4))
  plot(NULL, type = "n", xlim = tsp(h.pn)[1:2], ylim = range(h.pn[, 5:95]), main = "h_t")
  fan(h.pn)
  plot(NULL, type = "n", xlim = tsp(sigma.pn)[1:2], ylim = range(sigma.pn[, 1:95]), main = "sigma_t")
  fan(sigma.pn) 
}
}
