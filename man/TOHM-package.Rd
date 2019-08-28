\name{TOHM-package}
\alias{TOHM-package}
\alias{TOHM}
\docType{package}
\title{\packageTitle{TOHM}}
\description{\packageDescription{TOHM}}
\details{
The functions collected in \code{TOHM} mainly focus on the implementation of the Likelihood Ratio Tests (see \code{\link{TOHM_LRT}}). However, several functions (e.g.,\code{\link{EC_T}}, \code{\link{global_p}} ) can be used to obtain global p-values for other test statistics and to compute the Euler characteristic using the graph algorithm described in Algeri and van Dyk (2018).}
\author{
\packageAuthor{TOHM}
Maintainer: \packageMaintainer{TOHM}
}
\references{S. Algeri and D.A. van Dyk. Testing one hypothesis multiple times: The multidimensional case. arXiv:1803.03858, submitted to the Journal of Computational and Graphical Statistics, 2018.}
\keyword{ Testing One Hypothesis Multiple TImes }
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
