\name{EC_T}
\alias{EC_T}
\title{
Compute the Euler characteristic for a given  field
}
\description{
It computes the Euler characteristic (EC) of a given  field above specified thresholds over a specified search area.}
\usage{
EC_T(ck, Ts, THETA)
}
\arguments{
  \item{ck}{
Vector of thresholds defining the excursions sets with respect to which the ECs are computed.
}
  \item{Ts}{
 Vector of values of the field for each grid point in \code{THETA}.
}
  \item{THETA}{
A vector or matrix of   grid values for the nuisance parameter with respect to which the search is performed.
}
}

\value{
Returns a vector of EC values with respect to the  thresholds specified in \code{ck}.
}
\references{
S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional case.   arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.}
\author{
Sara Algeri}



\seealso{
\code{\link{EC_LRT}}
}
\examples{
EC_T(ck=c(3,4),Ts=rnorm(10), THETA=cbind(1:10,21:30))
}
\keyword{ Euler characteristic }
\keyword{ Generic random field }
