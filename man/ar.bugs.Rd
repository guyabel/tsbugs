\name{ar.bugs}
\alias{ar.bugs}
\title{
Create BUGS Script of a Autoregressive (AR) Time Series Model
}
\description{
Create BUGS script of an AR time series model, where the data is assumed to be normally distributed. Options allow for the inclusion of a different lag orders for the mean term, forecasts, posterior simulations from the model and alternative specification of prior distributions on each parameter.
}
\usage{
ar.bugs(y, ar.order = 1, k = NULL, sim = FALSE, mean.centre = FALSE, beg = ar.order + 1, 
  mean.prior = ar.prior, ar.prior = "dnorm(0,1)", 
	tol.prior = "dgamma(0.000001,0.000001)", var.prior = NULL, sd.prior = NULL, 
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
Include a term to centre the data on its mean value.Default is \code{FALSE}.
}
  \item{beg}{
Starting value for which data are considered onwards (and including) in the likelihood of the BUGS model. Default is \code{ar.order+1} but if comparing models of different orders, users may wish to set all \code{beg} to the same value.
}
  \item{mean.prior}{
Prior for mean term (not used if mean.centre is not set to \code{TRUE}). The distribution should be stated in BUGS syntax. By default, the same prior as the autoregressive terms is used.
}
  \item{ar.prior}{
Prior for autoregressive terms. The distribution should be stated in BUGS syntax. By default this is set to a normal distribution with mean 0 and tolerance 1 (\code{dnorm(0,1)}). The same prior are used for all autoregressive terms.
}
  \item{tol.prior}{
Prior for the tolerance of the model. The distribution should be stated in BUGS syntax. By default this is set to a uninformative gamma distribution. 
}
  \item{var.prior}{
Prior for the variance of the model. This must be a distribution of syntax recognisable to BUGS. This is not set by default (where a gamma distribution for the tolerance is used).
}
  \item{sd.prior}{
Prior for the standard deviation of the model. This must be a distribution of syntax recognisable to BUGS. This is not set by default (where a gamma distribution for the tolerance is used).
}
  \item{space}{
Include some additional empty lines to seperate the likelihood, priors, forcasts and simulation components of the BUGS model.
}
}
\details{
This function creates BUGS scripts of an AR time series model. Prior distributions should be set up using BUGS syntax. For example, \code{dnorm(.,.)} is a normal distribution with mean and tolerance (not variance) arguments. Only one argument of the \code{tol.prior}, \code{var.prior} or \code{sd.prior} should be set. Also, any new prior for either of these parameters should be based on distributions that provide only positive values, otherwise the BUGS model is likely to fail.

The data, \code{y}, can contain missing values. Note, if missing values are close the beginning of the series when a high order model specified (i.e. the second data point is missing and a AR(4) is specified) the user with have to set a high starting point for model to be fitted on (\code{beg}) for the BUGS model to function (i.e. 7).

}
\value{
\item{bug }{A BUGS model of type \code{tsbugs}.}
\item{data }{The data to be used with the model. This might extend the original data passed to the function with \code{k} unknown future values to be forecast.}
\item{info }{Additional information on the length of the data, variance type and line numbers of certain parts of the BUGS model.}
}
\author{
Guy J. Abel
}
\seealso{
\code{\link{sv.bugs}}, \code{\link{rv.bugs}}
}
\examples{
# Create AR(1) model for Lake Huron data
plot(LakeHuron)
LH <- LakeHuron
ar1 <- ar.bugs(y = diff(LH), ar.order = 1)
print(ar1)

# AR(2) model with alternative prior
ar2 <- ar.bugs(y = diff(LH), ar.order = 2, ar.prior = "dunif(-1,1)", var.prior = "dgamma(0.001,0.001)")
print(ar2)

# AR(3) model with forecast and posterior simulations
ar3 <- ar.bugs(y = diff(LH), ar.order = 3, k = 10, sim = TRUE, mean.centre = TRUE)
print(ar3) 

\dontrun{
# Run in OpenBUGS
writeLines(ar3$bug, "ar3.txt")
library("R2OpenBUGS")

ar3.bug <- bugs(data = ar3$data, 
				inits = list(inits(ar3)), 
				param = c(nodes(ar3, "prior")$name, "y.sim", "y.new"), 
				model = "ar3.txt", 
				n.iter = 11000, n.burnin = 1000, n.chains = 1)

# Plot the parameters posteriors and traces
library("coda")
param.mcmc <- as.mcmc(ar3.bug$sims.matrix[, nodes(ar3, "prior")$name])
plot(param.mcmc)

# Plot posterior simulations using fanplot
library("fanplot")
y.mcmc <- ar3.bug$sims.list$y.sim
y.pn <- pn(y.mcmc, start = tsp(diff(LH))[1] + 1)
plot(diff(LH), type = "n")
fan(y.pn)
lines(diff(LH)) 
}
}
