\name{ECden_vec}
\alias{ECden_vec}
\title{Compute the Euler characteristic densities}
\description{Compute the Euler characteristic (EC) densities for Gaussian, Chi-square and Chi-bar-square random fields up to a given order.}
\usage{
ECden_vec(c, D, type = c("Gaussian", "Chi^2", "Chi-bar^2"),  
k = NULL, k_vec = NULL, weights = NULL)
}
\arguments{
  \item{c}{Value on which the EC  densities are evaluated. }
  \item{D}{Maximum order of the EC density to be computed.}
  \item{type}{Type of random field. The possible options are \code{"Gaussian"}, \code{"Chi^2"}, and \code{"Chi-bar^2"}. See details.}
  \item{k}{If \code{type="Chi^2"}, degrees of freedom of the Chi-square random field.}
  \item{k_vec}{If \code{type="Chi-bar^2"} the degrees of freedom  of the Chi-square random fields  in the mixture.}
  \item{weights}{ If \code{type="Chi-bar^2"} the   weights of the mixure of Chi-square random fields in the mixture. The ordering should be the same as in \code{k_vec}.}}
  
\details{If \code{type="Chi-bar^2"} the degrees of freedom  of the Chi-square random fields involved in the mixture, as well as the respective weights, have to be spefcified in the arguments \code{k_vec} and \code{weights}. }
\value{Returns the values of  the EC densities of order zero up to the dimension of the search area considered and evaluated at \code{c}.}
\references{
R.J. Adler and J.E. Taylor. Random fields and geometry. Springer Science and Business
Media, 2009.

J.E. Taylor and K.J.Worsley. Detecting sparse cone alternatives for gaussian random
fields, with an application to fmri. Statistica Sinica, 2013.}
\author{Sara Algeri}
\seealso{\code{\link{chi2_ECden}}, \code{\link{Gauss_ECden}}}
\examples{ECden_vec(12,2,"Chi-bar^2",k_vec=c(0,1),weights=c(0.5,0.5))}
\keyword{Euler characteristic densities}