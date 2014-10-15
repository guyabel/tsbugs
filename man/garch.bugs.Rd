\name{garch.bugs}
\alias{garch.bugs}
\title{
Create BUGS Script of a Gerneralised Autoregressive Conditional Hetroskedastic (GARCH) Time Series Model
}
\description{
Create BUGS script of an GARCH time series model, where the data is assumed to be normally distributed. Options allow for the inclusion of a different lag orders for the mean and variance term, forecasts, posterior simulations from the model and alternative specification of prior distributions on each parameter.
}
\usage{
garch.bugs(y, ar.order = 0, k = NULL, sim = FALSE, mean.centre = FALSE, beg = ar.order + 1, mean.prior = ar.prior, ar.prior = "dnorm(0,1)", arch.order = 1, garch.order = 1, vol.beg = max(arch.order, garch.order) + 1, h0.prior = "dgamma(0.000001,0.000001)", vol.mean.prior = "dgamma(0.000001,0.000001)", arch.prior = "dgamma(0.000001,0.000001)", garch.prior = "dgamma(0.000001,0.000001)", space = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{arch.order}{
MA order of the volatility process for BUGS model.
}
  \item{garch.order}{
AR order of the volatility process for BUGS model.
}
  \item{vol.beg}{
Starting value for which the volatility process in the data are considered onwards (and including) in the likelihood of the BUGS model. Default is \code{max(arch.order+garch.order + 1} but if comparing models of different orders, users may wish to set all \code{beg} to the same value.
}
  \item{h0.prior}{
Prior distribution for the inital latent volatility. This must be a distribution of syntax recognisable to BUGS. By default this is set to a uninformative gamma distribution.
}
  \item{vol.mean.prior}{
Prior distribution for the mean of the volatility. This must be a distribution of syntax recognisable to BUGS. By default this is set to a uninformative gamma distribution.
}
  \item{arch.prior}{
Prior distribution for the arch coefficients. This must be a distribution of syntax recognisable to BUGS. By default this is set to a uninformative gamma distribution.
}
  \item{garch.prior}{
Prior distribution for the garch coefficients. This must be a distribution of syntax recognisable to BUGS. By default this is set to a uninformative gamma distribution.
}
  \item{space}{
Include some additional empty lines to seperate the likelihood, priors, forcasts and simulation components of the BUGS model.
}
}
\details{
This function create BUGS scripts of a standard GARCH model. Prior distributions for the mean, autoregressive and moving average coefficients for the volatiltiy process can be specified by users to differ from their default values. User specified prior distributions should be character strings following BUGS syntax. For example, \code{dnorm} is a normal distribution with mean and tolerance (not variance) arguments.

The data \code{y}, can contain missing values. Note, if missing values are close the beginning of the series when a high order model for the mean process is specified (i.e. the second data point is missing and a AR(4) is specified) the user with have to set a high starting point for model to be fitted on (\code{beg}) for the BUGS model to function (i.e. 7).
}
\value{
\item{bug }{A BUGS model of type \code{tsbugs}.}
\item{data }{The data to be used with the model. This might extend the original data passed to the function with \code{k} unknown future values to be forecast.}
\item{info }{Additional information on the length of the data, variance type and line numbers of certain parts of the BUGS model.}
}
\references{
Congdon, P. (2006). \emph{Bayesian Statistical Modelling} 2nd edition, John Wiley.
}
\author{
Guy J. Abel
}
\seealso{
\code{\link{ar.bugs}}, \code{\link{sv.bugs}}, \code{\link{rv.bugs}}
}
\examples{
# Create AR(0)-GARCH(1,1) model for svpdx
y <- svpdx$pdx
plot(y, type = "l")
garch11 <- garch.bugs(y)
print(garch11)

# AR(1)-ARCH(2) structure in the volatility
arch20 <- garch.bugs(y, ar.order = 1, arch.order = 2, garch.order=0)
print(arch20)

# AR(0)-GARCH(1,1) model with alternative prior
garch11 <- garch.bugs(y, garch.prior = "dgamma(0.1,0.1)")
print(garch11)

# AR(0)-GARCH(1,1) model with forecast and posterior simulations
garch11 <- garch.bugs(y, k = 10, sim = TRUE)
print(garch11) 
}
