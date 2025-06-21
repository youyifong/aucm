\name{rlogit}
\alias{rlogit}
\alias{coef.rlogit}
\alias{predict.rlogit}
\alias{ratio.rlogit}
\alias{trainauc.rlogit}
\alias{logistic.f}
\title{rlogit}
\description{
Robust logistic regression estimator of Bianco and Yohai  
}
\usage{
rlogit (formula, dat, const=0.5, kmax=1e3, maxhalf=10, verbose=FALSE)
\S3method{coef}{rlogit}(object,\dots)
\S3method{trainauc}{rlogit}(fit, training.data=NULL, \dots)
\S3method{predict}{rlogit}(object, newdata, \dots)
\S3method{ratio}{rlogit}(fit)
logistic.f(eta,h,loss=TRUE)
}
\arguments{

  \item{formula}{a formula specifying the model to be fit.}
  \item{dat}{a data frame containing the outcome and covariates in the model}
  \item{const}{tuning constant used in the computation of the estimator, defaults to 0.5}
  \item{kmax}{maximum number of iterations before convergence, defaults to 1000}
  \item{maxhalf}{max number of step-halving ,defaults to 10}
  \item{verbose}{logical}
  \item{object}{an object of class 'rlogit'}
  \item{fit}{an object that inherits from class 'auc' such as 'rauc' or 'sauc'}
  \item{newdata}{data at which to predict}
  \item{training.data}{data frame used to compute auc based on a fit obtained by a call to  \code{rauc}, \code{sauc} or \code{sauc.dca}}
  \item{eta,h}{logistic.f computes for loss = FALSE expit(eta/h) or expit(-eta/h) for loss = TRUE}
  \item{loss}{a boolean. if TRUE (default) logistic loss is assumed.}
  \item{\dots}{arguments passed to or from methods}
  
}

\value{
 A list with the follwoing components:
 \item{convergence}{logical, was convergence achieved}
 \item{objective}{value of the objective function at the minimum}
 \item{coef}{estimates for the parameters}
 \item{sterror}{standard errors of the parameters (if convergence is TRUE)}
 }
 
\details{
    This program computes the estimator of Bianco and Yohai (1996) in
    logistic regression. By default, an intercept term is included
    and p parameters are estimated. The outcome is coded as a 0/1 binomial variable. 

    If initwml == TRUE, a weighted ML estimator is
    computed with weights derived from the MCD estimator
    computed on the explanatory variables. If initwml == FALSE,
    a classical ML fit is perfomed.
    When the explanatory variables contain binary observations,
    it is recommended to set initwml to FALSE or to modify the
    code of the algorithm to compute the weights only on the
    continuous variables. 
}
\references{
\emph{Implementing the Bianco and Yohai estimator for Logistic Regression}\cr
Croux, C., and Haesbroeck, G. (2003)\cr
Computational Statistics and Data Analysis, 44, 273-295 \cr
}
\author{
Christophe Croux, Gentiane Haesbroeck. 
Thanks to Kristel Joossens and Valentin Todorov for improving the code.
}

\examples{

set.seed(1)
x0 <- matrix(rnorm(100,1))
y  <- as.numeric(runif(100)>0.5)        # numeric(runif(100)>0.5)
dat=data.frame(y=y, x=x0)
rlogit(y~x, dat)




}
