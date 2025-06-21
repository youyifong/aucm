\name{ramp.f}
\alias{ramp.f}
\alias{ramp.f}
\title{ramp.f}
\description{
computes ramp function value from paired difference of linear combinations
}
\usage{
ramp.f (eta,s,loss=TRUE)
}
\arguments{
  \item{eta}{a vector of paired difference of linear combinations}
  \item{s}{absolute value of the slope parameter}
  \item{loss}{a boolean. If TRUE, return loss function i.e. 1 - RAUC. If FALSE, rturn RAUC. Default to TRUE, because we minimize loss.}
}

\value{
  A vector of sames size as eta
}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}

\examples{

dat = sim.dat.1(n=100,seed=1)
X1 = as.matrix(subset(dat, y==0, select=c(x1,x2)))
X2 = as.matrix(subset(dat, y==1, select=c(x1,x2)))
X.diff = get.X.diff (X1, X2)
dim(X1)
dim(X2)
dim(X.diff)
aux = ramp.f(X.diff \%*\% c(1,1), s=1)
length(aux)
mean(aux)
aux = ramp.f(X.diff \%*\% c(1,1), s=1, loss=FALSE)
length(aux)
mean(aux)

}
