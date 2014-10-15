\name{inits}
\alias{inits}
\title{
Produce a Set of Candidate Initial Values
}
\description{
Creates a list of initial values for a \code{tsbugs} object to help shorten your code when running models through R2WinBUGS or R2OpenBUGS.
}
\usage{
inits(bug, warn.mess = TRUE)
}
\arguments{
  \item{bug}{
A time series BUGS model created using the \code{tsbugs} package.
}
  \item{warn.mess}{
Print warning message notifying users that the initial values simple guesses and have nothing to do with data. Default is \code{TRUE}.
}
}
\details{
The \code{inits} function is intended to provide a set of candidate initial values for BUGS models created using the \code{tsbugs} package. The list provided from \code{inits} does not guarantee that BUGS model will run from an efficient set of starting points, or indeed that the BUGS model will run at all. Values are set to be in the parameter space and are not based on the data used to create the \code{tsbugs} model. If problems do occur running the BUGS model, users are advised to set their own initial values.
}
\value{
\code{list} with components dependent on model type.
}
\author{
Guy J. Abel
}

\seealso{
\code{\link{ar.bugs}}, \code{\link{sv.bugs}}, \code{\link{rv.bugs}}
}
\examples{
# Create AR(4) model for Lake Huron data
LH <- LakeHuron
ar4 <- ar.bugs(y = diff(LH), ar.order = 4)

# Candidate initial values
inits(ar4) 
}
