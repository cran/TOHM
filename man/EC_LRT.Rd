\name{EC_LRT}
\alias{EC_LRT}
\title{
Compute the Euler characteristic for the generalized Likelihood Ratio Test field
}
\description{
It computes the Euler characteristic (EC) of the generalized  Likelihood Ratio test (LRT)  field above  specified thresholds over a given search area.}
\usage{
EC_LRT(ck, x, mll, null0, init, lowlim, uplim, THETA)
}
\arguments{
  \item{ck}{
Vector of thresholds defining the excursions sets with respect to which the ECs are computed.
}
  \item{x}{
A vector or matrix collecting the data on which the LRT is computed.
}
  \item{mll}{
A function specifying the negative (profile) log-likelihood. See details. }
  \item{null0}{
A scalar or vector collecting the values of the free parameters under the null hypothesis. See details.}
  \item{init}{
Vector of initial values for the MLE.}
  \item{lowlim}{
Vector of lower bounds for the MLE.}
  \item{uplim}{
Vector of upper bounds for the MLE.
}
  \item{THETA}{
A vector or matrix of   grid values for the nuisance parameter with respect to which the search is performed.
}}
\details{
\code{mll}  takes as first argument the vector of the parameters for which   the  MLE is generated. Other arguments of \code{mll}  are the data vector or matrix (\code{x}) and a scalar or vector corresponding to the fixed value for the nuisance parameter with respect to which the profilying is computed (\code{theta}, see \code{\link{gLRT}}). If the latter is a vector it must be of same length of the rows in \code{THETA}.
If the null model has nuisance parameters, \code{null0} takes as arguments the values of the parameters being tested under the null hypothesis, followed by the estimates of the nuisance parameters obtained assuming that the null hypothesis is true.}
\value{
Returns a vector of EC values with respect to the  thresholds specified in \code{ck}.
}
\references{
S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional case.   arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.}
\author{
Sara Algeri}

\seealso{
\code{\link{EC_T}}
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

EC_LRT(ck=c(1,8),x=data,mll=mll,null0=0,init=c(0.1),
lowlim=c(0),uplim=c(1), THETA_sel)
}
\keyword{ Euler characteristic }
\keyword{ generalized Likelihood Ratio random field }
