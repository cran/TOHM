\name{chi2_ECden}
\alias{chi2_ECden}
\title{Compute Euler characteristic density for Chi-square random fields 
}
\description{
Computes the Euler characteristic (EC) density of a given order for Chi-squared random fields.
}
\usage{
chi2_ECden(c, k, j)
}
\arguments{
  \item{c}{Value on which the EC  density is evaluated. 
}
  \item{k}{Degrees of freedom of the Chi-square random field.
}
  \item{j}{
Order of the EC density to be implemented.
}
}

\value{
Returns the value of the EC density of order \code{j} evaluated at \code{c} for a Chi-square random field with \code{k}
degrees of freedom.
}
\references{
R.J. Adler and J.E. Taylor. Random fields and geometry. Springer Science and Business
Media, 2009.}
\author{
Sara Algeri
}

\seealso{
\code{\link{Gauss_ECden}}, \code{\link{ECden_vec}}
}
\examples{
c<-1
k<-1
j<-2
chi2_ECden(c,k,j)
}

\keyword{ Euler characteristic densities }
\keyword{ Chi-squared random field }
