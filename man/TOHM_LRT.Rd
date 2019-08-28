\name{TOHM_LRT}
\alias{TOHM_LRT}

\title{Compute the Likelihood Ratio Test under non-standard conditions.}

\description{It implements the procedure described in Algeri and van Dyk (2018) to perform tests of hypothesis under non-regular conditions, and which can be formulated as test of hypothesis where a nuisance parameter is present only under the alternative.}

\usage{TOHM_LRT(x, mll, null0, init, lowlim, uplim, THETA, ck,
type = c("Chi^2", "Chi-bar^2"), k = NULL, k_vec = NULL,
weights = NULL, ECdensities = NULL, ECs = NULL)}

\arguments{
\item{x}{A vector or matrix collecting the data. }
    \item{mll}{
A function specifying the negative (profile) log-likelihood. See details. }
  \item{null0}{
A   vector or scalar of the free parameters under the null hypothesis. See details.}
  \item{init}{
A   vector or scalar  of initial values for the MLE.}
  \item{lowlim}{
A   vector or scalar of lower bounds for the MLE.}
  \item{uplim}{
A   vector or scalar of upper bounds for the MLE.
}
  \item{THETA}{
A vector or matrix of   grid values of the nuisance parameter with respect to which the search is performed.
}
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

\details{\code{mll}  takes as first argument the vector of the parameters for which   the  MLE is generated. Other arguments of \code{mll}  are the data vector or matrix (\code{x}) and a scalar or vector corresponding to the fixed value for the nuisance parameter with respect to which the profilying is computed (\code{theta}, see \code{\link{gLRT}}). If the latter is a vector it must be of same length of the rows in \code{THETA}.
If the null model has nuisance parameters, \code{null0} takes as arguments the values of the parameters being tested under the null hypothesis, followed by the estimates of the nuisance parameters obtained assuming that the null hypothesis is true.
If \code{type="Chi-bar^2"} the degrees of freedom  of the Chi-square random fields involved in the mixture, as well as the respective weights, have to be spefcified in the arguments \code{k_vec} and \code{weights}.
If the distribution of the random field is not available in \code{type}, the user can specify in \code{ECdensities} a function taking \code{c} as argument and returning
the vector of the desired EC densities  to be evaluated at \code{c}. Notice that the length of the vector returned by the function specified in \code{ECdensities} must corresponds to one plus the dimension of the search area (since the first value should correspond to the EC density of order zero (see \code{\link{ECden_vec}})). }

\value{
\item{max_gLRT }{Maximum observed of the LRT field.}
\item{theta_max }{Value of \code{THETA} at which the maximum is observed.}
\item{global_p }{Global p-value.}
\item{MCerror }{Monte Carlo error associated to the global p-vaue.}
}
\references{S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional
case.   arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.}
\author{Sara Algeri}


\seealso{
\code{\link{find_max}}, \code{\link{global_p}},\code{\link{EC_T}}.
}
\examples{
#generating data of interest
N<-100
x<-as.matrix(cbind(runif(N*2,172.5,217.5),runif(N*2,-2,58)))
x2<-x[(x[,1]<=217.5)&(x[,1]>=172.5),]
x_sel<-x2[(x2[,2]<=(28+sqrt(30^2-(x2[,1]-195)^2)))&(x2[,2]>=(28-
sqrt(30^2-(x2[,1]-195)^2))),]
data<-x_sel[sample(seq(1:(dim(x_sel)[1])),N),]

#Specifying minus-log-likelihood
kg<-function(theta){integrate(Vectorize(function(x) {
exp(-0.5*((x-theta[1])/0.5)^2)*integrate(function(y) {
exp(-0.5*((y-theta[2])/0.5)^2) }, 28-sqrt(30^2-(x-195)^2),
28+sqrt(30^2-(x-195)^2))$value}) , 172.5, 217.5)$value}
mll<-function(eta,x,theta){
  -sum(log((1-eta)/(pi*(30)^2)+eta*exp(-0.5*((x[,1]-
  theta[1])/0.5)^2-
  0.5*((x[,2]-theta[2])/0.5)^2)/kg(theta)))}

#Specifying search region
theta1<-seq(172.5,217.5,by=15)
theta2<-seq(-2,58,by=10)
THETA<-as.matrix(expand.grid(theta1,theta2))
originalR<-dim(THETA)[1]
rownames(THETA)<-1:(dim(THETA)[1])
THETA2<-THETA[(THETA[,1]<=217.5)&(THETA[,1]>=172.5),]
THETA_sel<-THETA2[(THETA2[,2]<=(28+sqrt(30^2-(THETA2[,1]-
195)^2)))&(THETA2[,2]>=(28-sqrt(30^2-(THETA2[,1]-195)^2))),]

#Generating toy EC
ECs<-cbind(rpois(100,1.5),rpois(100,1))

TOHM_LRT(data,mll,null0=0,init=c(0.1),lowlim=c(0),uplim=c(1),
THETA=THETA_sel,ck=c(1,8),type=c("Chi-bar^2"),
k=NULL,k_vec=c(0,1),weights=c(0.5,0.5),
ECdensities=NULL,ECs=ECs)
  }

\keyword{ Testing One HYpothesis Multiple times }
\keyword{ Non regularity }
