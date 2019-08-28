\name{gLRT}
\alias{gLRT}
\title{
Compute the generalized Likelihood Ratio Test}
\description{
Compute the generalized Likelihood Ratio Test (LRT) for a specified value of the nuisance parameter.}
\usage{
gLRT(theta, mll, x, init, lowlim, uplim, null0)
}
\arguments{
  \item{theta}{
A vector or scalar of the value of the nuisance parameter with respect to which the LRT is computed.}
  \item{mll}{
A function specifying the negative (profile) log-likelihood. See details. }
  \item{x}{
A vector or matrix collecting the data.
}
  \item{init}{
A   vector or scalar  of initial values for the MLE.}
  \item{lowlim}{
A   vector or scalar of lower bounds for the MLE.}
  \item{uplim}{
A   vector or scalar of upper bounds for the MLE.
}
  \item{null0}{
A   vector or scalar of the free parameters under the null hypothesis. See details.}}
\details{
\code{mll}  takes as first argument the vector of the parameters for which   the  MLE is generated. Other arguments of \code{mll}  are the data vector or matrix (\code{x}) and a scalar or vector corresponding to the fixed value for the nuisance parameter with respect to which the profilying is computed (\code{theta}, see \code{\link{gLRT}}). If the latter is a vector it must be of same length of the rows in \code{THETA}.
If the null model has nuisance parameters, \code{null0} takes as arguments the values of the parameters being tested under the null hypothesis, followed by the estimates of the nuisance parameters obtained assuming that the null hypothesis is true.}
\value{
The value of the generalized LRT for a specified value of \code{theta}.
}
\references{
S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional
case.   arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.

A.C. Davison. Statistical models, volume 11. Cambridge University Press, 2003.
}
\author{
Sara Algeri
}

\seealso{
\code{\link{find_max}}, \code{\link{TOHM_LRT}}
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

gLRT(theta=c(200,30),mll=mll,init=0.1,lowlim=0,uplim=1,null0=0,x=data)
}

\keyword{Generalized Likelihood Ratio Test }
\keyword{Profile Likelihood }
