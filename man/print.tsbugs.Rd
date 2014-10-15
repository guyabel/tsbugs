\name{print.tsbugs}
\alias{print.tsbugs}
\title{
Prints tsbugs object
}
\description{
Prints the BUGS model in a \code{tsbugs} object.
}
\usage{
\method{print}{tsbugs}(x, ...)
}
\arguments{
  \item{x}{
A BUGS model of class \code{tsbugs}
}
  \item{\dots}{
Additional arguments passed to \code{print}
}
}
\details{
Prints \code{tsbugs} BUGS model. Note, in a \code{tsbugs} object there are additional details (the data used and other summary information) that are not printed.
}
\author{
Guy J. Abel
}
\examples{
# Create AR(1) model for Lake Huron data
LH <- LakeHuron
ar1 <- ar.bugs(y = diff(LH), ar.order = 1)

print(ar1) 
}
