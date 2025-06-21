\name{tune.it}
\alias{tune.it}
\alias{get.gamma0}
\title{tune.it}
\description{
Tuning methods.
}
\usage{
tune.it (formula, dat.train, dat.tune, method, kernel, verbose=TRUE, step.size = 2)
get.gamma0 (formula, dat.train)

}
\arguments{
  \item{formula}{a formula object.}
  \item{dat.train}{a data frame. Training data}
  \item{dat.tune}{a data frame. Tuning data. If NULL, gacv is done.}
  \item{method}{a string. "svm" or "rauc"}
  \item{kernel}{a string. "rbf" or "linear"}
  \item{verbose}{logical}
  \item{step.size}{step size for lambda as in step.size^seq(-3,3,1), please see code for more details}
}

\value{
  \code{tune.it} returns a vector of one or two elements: lamcost and gamma, depending on the kernel. \cr
  \code{get.gamma0} returns .5/quantile(dist., c(.9,.75,.5)), where dist. is the Eucledian distance between objects from two classes.
}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{ksebesty@gmail.com}\cr
}
