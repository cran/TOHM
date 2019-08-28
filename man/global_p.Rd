\name{global_p}
\alias{global_p}

\title{Compute global p-values}
\description{It computes the global p-value for a given value of the test statistic.}
\usage{global_p(c, ck, type = c("Gaussian", "Chi^2", "Chi-bar^2"),  
k = NULL, k_vec = NULL, weights = NULL, ECdensities = NULL, ECs)}
\arguments{
  \item{c}{Observed value of the test statistic.}
  \item{ck}{Vector of thresholds defining the excursions sets with respect to which the ECs are computed.
}
  \item{type}{Type of random field. The possible options are \code{"Gaussian"}, \code{"Chi^2"}, and \code{"Chi-bar^2"}. See details.}
   \item{k}{If \code{type="Chi^2"}, degrees of freedom of the Chi-square random field.
}
  \item{k_vec}{If \code{type="Chi-bar^2"} the degrees of freedom  of the Chi-square random fields  in the mixture.
}
  \item{weights}{ If \code{type="Chi-bar^2"} the   weights of the mixure of Chi-square random fields in the mixture. The ordering should be the same as in \code{k_vec}.
}
  \item{ECdensities}{See datails.}
  \item{ECs}{A vector or matrix containing the Euler characteristics (ECs) computed over a Monte Carlo simulation of the random field under the null model. Each colum correspond to the ECs obtained for each of the thresholds in \code{ck}.}
}
\details{If \code{type="Chi-bar^2"} the degrees of freedom  of the Chi-square random fields involved in the mixture, as well as the respective weights, have to be spefcified in the arguments \code{k_vec} and \code{weights}.  
If the distribution of the random field is not available in \code{type}, the user can specify in \code{ECdensities} a function taking \code{c} as argument and returning
the vector of the desired EC densities  to be evaluated at \code{c}. Notice that the length of the vector returned by the function specified in \code{ECdensities} must corresponds to one plus the dimension of the search area (since the first value should correspond to the EC density of order zero (see \code{\link{ECden_vec}})).
}
\value{
\item{global_p }{Global p-value.}
\item{MCerror }{Monte Carlo error associated to the global p-vaue.}
}
\references{
S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional
case.   arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.}
\author{
Sara Algeri}



\seealso{
\code{\link{find_max}},\code{\link{TOHM_LRT}},\code{\link{ECden_vec}}}
\examples{
ck<-c(1,2)
ECs<-cbind(rpois(100,1.5),rpois(100,1))
global_p(c=12,ck=ck,type="Gaussian",ECs=ECs)
}

\keyword{Global p-value }
\keyword{Generic random field }
