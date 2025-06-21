\name{sauc.phi}
\alias{sauc.phi}
\alias{phi.f}
\alias{vcov.sauc}
\title{sauc.phi}
\description{
\code{sauc.phi} optimizes Normal CDF approximation of AUC using Newton Raphson
}
\usage{
sauc.phi (formula,dat,constrain.method="L2",h.method="Lin",start.method="rlogit",
opt.method = "Lin", upper = NULL, verbose = FALSE, ret.vcov = FALSE, truth = NULL, 
    beta.init=NULL)
}
\arguments{
    \item{formula}{a formula}
    \item{dat}{a data frame}
    \item{constrain.method}{a string. If "L2", L2 norm is constrained to 1. If "beta1", beta1 is fixed to 1. Default "L2".}
    \item{h.method}{a string. If "Lin", Lin et al, data dependent. If "Vexler", (n1*n2)^(-0.1) Vexler et al (2006). If "MH", Ma and Huang. Default "Lin".}
    \item{start.method}{a string. If "rlogit", robust logistic fit is used as beta.init If "1", a vector of 1 is used as beta.init. Default "rlogit".}
    \item{opt.method}{character string, possible values are "truth","YH","Lin", please see code for more details}
    \item{upper}{required for opt.method = 'YH'}
    \item{verbose}{logical}
    \item{ret.vcov}{logical, whether to return an estimate of the covariance matrix of 'beta' for normal or logistic sigmoid functions.}
    \item{truth}{numeric, it will be returned as the result of the fit, please see code for more details}
    \item{beta.init}{vector. Initial values for coefficients.}
}

\details{
If an error happens during optimization (typically due to solve()), the errors are catched and NAs are returned.
}

\author{
  Shuxin Yin \email{}\cr
  Ying Huang \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
}

\examples{

seed=26
seed=16
seed=3
dat.train = sim.dat.1(n=200, seed=seed, add.outliers=TRUE) 
fits=list()
fits[[1]]=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="Lin")
fits[[2]]=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="MH")
fits[[3]]=sauc.phi(y~x1+x2, dat.train,constrain.method="beta1",h.method="Lin")
fits[[4]]=sauc.phi(y~x1+x2, dat.train,constrain.method="beta1",h.method="MH") 
# not a good combination of constrain.method and h.method
sapply(fits, function(x) ratio(x)[2])


# explosion
seed=954
dat.train = sim.dat.1(n=200, seed=seed, add.outliers=TRUE) 
fit.1 = sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="Lin")
ratio(fit.1)


}
