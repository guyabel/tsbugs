\name{nodes}
\alias{nodes}
\title{
Data Frame of Nodes within a BUGS model.
}
\description{
Provides a \code{data.frame} of both random and deterministic nodes in a \code{tsbugs} BUGS model. Can also be used to extract nodes from just a single part of the model, such as the prior distributions. This might be of particular use to users when setting up parameters to follow when running models through R2WinBUGS or R2OpenBUGS.
}
\usage{
nodes(bug, part = NULL)
}
\arguments{
  \item{bug}{
A time series BUGS model created using the \code{tsbugs} package.
}
  \item{part}{
A part of the \code{tsbugs} models. Must be one of \code{"likelihood"}, \code{"priors"}, \code{"forecast"} or \code{"simulation"}. Abbreviations of these parts are allowed (\code{pmatch} is used for the matching).
}
}
\details{
Provides a summary of nodes in all of, or part of, the  \code{tsbugs} BUGS model. Returns a data frame with details on each node. When nodes are within a loop additional details on the starting and ending time points for the node are also given.
}
\value{
A \code{data.frame} with columns:
 \item{name }{Name of node.}
 \item{type }{Stochastic node represented by \code{~}. Deterministic node represented by \code{<-}.}
 \item{dt }{Distribution or transformation in node.}
 \item{beg }{Beginning of node loop.}
 \item{end }{End of node loop.}
 \item{stoc }{Binary representation of type column, where 1 represents a stochastic node and 0 a deterministic node.}
 \item{id }{Line number of node in the BUGS model.}
 \item{dist }{Distribution in node.}
 \item{param1 }{First parameter of distribution.}
 \item{param2 }{Second parameter of distribution.}
}
\author{
Guy J. Abel
}
\examples{
# Create AR(1) model for Lake Huron data
LH <- LakeHuron
ar1 <- ar.bugs(y = diff(LH), ar.order = 1)

# All nodes in model
nodes(ar1)

# Priors in model
nodes(ar1, "priors") 
}

