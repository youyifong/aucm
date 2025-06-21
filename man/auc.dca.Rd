\name{dcsauc}
\alias{dcsauc}
\alias{srauc}
\alias{auc.dca}
\alias{smooth.rauc}
\alias{dcsauc.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AUC optimization with DCA
}
\description{
optimizes SAUC using a smoothed DCA algorithm.
}
\usage{
dcsauc (formula, data, ...)
srauc (formula, data, ...)
auc.dca (formula, data, 
    type="srauc",
    kernel="linear", para=NULL,  
    lambda=.1, zeta=.1, b=10, s=1, epsilon=1e-3, 
    method="tron", decomposition=TRUE,
    dca.control =  list(maxit=1e3, abstol=1e-5, coef.init=NULL, lincomb.init=NULL),
    tron.control = list(q=50, maxfev=1e3, gtol=1e-2, frtol=1e-12, K.thresh=1, verbose=0),
    return.K=FALSE, verbose = FALSE
)
}
\arguments{
  \item{formula}{formula, e.g. y~x1+x2}

  \item{data}{a data frame}
  
  \item{type}{string. Either srauc or dcsauc}
  \item{kernel}{See \link{getK} for more details}
  
  \item{para}{See \link{getK} for more details}

  \item{lambda}{scale parameter of the penalty function, defaults to 1}

  \item{zeta}{parameter (->0+) in writing sigmoid function as differene of two convex functions.}
  \item{b}{'decay rate' parameter in sigmoid function \eqn{1/(exp(bx))}}
  \item{s}{the parameter in rauc}
  \item{epsilon}{the parameter in the approximation of a hinge function}
  \item{method}{the optimizer to use, "tron", or an \code{optim} method}
  \item{decomposition}{Boolean. If TRUE, decomposition strategy is used if tron is the method}
  \item{dca.control}{list of control parameters for the DCA algorithm}
  \item{tron.control}{list of control parameters to 'tron' optimizer}
  \item{return.K}{logical, whether to return the Kernel matrix}
  \item{verbose}{logical, whether to print info as alg. progresses}
  \item{...}{parameters passed to auc.dca}
  
}
\details{
    dcsauc and srauc pass directly to auc.dca with the name-sake type. 
}
\examples{
#
#
#dat = sim.dat.1(n=100,seed=1)
#dat.test = sim.dat.1(n=1e3,seed=1000)
#
#t.1 = system.time({
#    fit.1=sauc.dca(y~x1+x2, dat, zeta=.1)
#})
#
#t.2 = system.time({
#    fit.2=sauc.dca(y~x1+x2, dat, zeta=1)
#})
#
## compare time
#rbind(t.1, t.2)[,3]
#
## compare performance
#RUnit::checkEqualsNumeric(
#    c(fit.1$train.auc, fit.2$train.auc)
#, c(0.7291917, 0.7282913), tolerance=1e-6)
#



}
