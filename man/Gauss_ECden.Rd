\name{Gauss_ECden}
\alias{Gauss_ECden}

\title{
Compute Euler characteristic density for Gaussian random fields 
}
\description{
Computes the Euler characteristic (EC) density of a given order for Gaussian random fields.
}
\usage{
Gauss_ECden(c, j)
}

\arguments{
  \item{c}{
Value on which the EC  density is evaluated. 
}
  \item{j}{
Order of the EC density to be implemented.
}
}

\value{
Returns the value of the EC density of order \code{j} evaluated at \code{c} for a Gaussian random field.
}
\references{
R.J. Adler and J.E. Taylor. Random fields and geometry. Springer Science and Business
Media, 2009.}
\author{
Sara Algeri
}

\seealso{
\code{\link{chi2_ECden}}, \code{\link{ECden_vec}}
}
\examples{
c<-1
j<-2
Gauss_ECden(c,4)
}

\keyword{ Euler characteristic densities }
\keyword{ Gaussian random field }