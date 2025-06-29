\name{grid.auc}
\alias{eauc}
\alias{grid.auc}
\alias{plot.gridsearch}
\alias{vcov.gridsearch}
\title{grid.auc}
\description{
grid search for beta that maximize (penalized, partial) auc/sauc/rauc
eauc maximizes empirical AUC, but only works with two covariates
}
\usage{
grid.auc (formula, dat, beta, approx.type=NULL, approx.param=1, lambda=0, loss=TRUE, 
    t0=NULL, t1=NULL,ret.vcov = FALSE)
eauc (formula, dat,t0 = NULL, t1 = NULL)
}
\arguments{
  \item{formula}{a formula}
  \item{dat}{a data frame}
  \item{beta}{a matrix of coefficients}
  \item{approx.type}{a string. If NULL, AUC is computed. If "phi", normal CDF approximation SAUC is computed. If "logistic", logistic approximation SAUC is computed. If "rauc", ramp AUC approximation is computed. Defaults to NULL}
  \item{approx.param}{'s' for \code{rauc}, 'h' for \code{sauc}}
  \item{loss}{a boolean. TRUE is default and means 1-(p)RAUC is computed. If lambda is not 0, loss is forced to be TRUE internally.}
  \item{t0}{a number between 0 and 1 that is the lower boundary of pAUC}
  \item{t1}{a number between 0 and 1 that is the upper boundary of pAUC}
  \item{lambda}{a number that scales the L2 penalty, default to 0, meaning no penalty. If lambda is not 0, loss is forced to be TRUE.}
  \item{ret.vcov}{logical, whether to return an estimate of the covariance matrix of 'beta' for normal or logistic sigmoid functions.}
}

\value{
  A n x n matrix
}
\details{
eauc is a shortcut for grid.auc when empirical AUC is the objective function.
When loss is FALSE, the criterion function is mean_i_j(loss)
When loss is TRUE, including when lambda is not 0, the criterion function is sum_i_j(loss) + 0.5 * lambda * pen, i.e. the penalty is added to the sum of penalty and not mean of the penalty
}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}

\examples{

library(aucm)

dat = sim.dat.1(n=200,seed=1)
beta=cbind(4, 4*seq(-1,0,length=100))
dim(beta)

fit = eauc(y~x1+x2, dat)

\dontrun{
# not run due to r cmd check requirements

par(mfrow=c(3,2))

out1 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=FALSE, t0=NULL, t1=NULL, lambda=0, 
    approx.type="rauc")
plot(out1$pars[2,]/out1$pars[1,],out1$vals,type="l",xlab=expression(beta[2]/beta[1]),main="RAUC")

# penalized RAUC
out2 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=TRUE, t0=NULL, t1=NULL, lambda=0, 
    approx.type="rauc")
out3 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=TRUE, t0=NULL, t1=NULL, lambda=30, 
    approx.type="rauc")
plot(out2$pars[2,]/out2$pars[1,],out2$vals,type="l",xlab=expression(beta[2]/beta[1]),
    main="penalized RAUC loss")
lines(out3$pars[2,]/out3$pars[1,],out3$vals,type="l",col=2)
out2$par
out3$par

# pRAUC
out4 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=FALSE, t0=0, t1=0.5, lambda=0, 
    approx.type="rauc")
out5 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=FALSE, t0=0.5, t1=1, lambda=0, 
    approx.type="rauc")
plot(out4$pars[2,]/out4$pars[1,],out4$vals,type="l",xlab=expression(beta[2]/beta[1]),main="pRAUC")
plot(out5$pars[2,]/out5$pars[1,],out5$vals,type="l",xlab=expression(beta[2]/beta[1]),main="pRAUC")
out4$par
out5$par


# penalized pRAUC
out6 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=TRUE, t0=0, t1=0.5, lambda=0, 
    approx.type="rauc")
out7 = grid.auc (y~x1+x2, dat, beta, approx.param=1, loss=TRUE, t0=0, t1=0.5, lambda=10, 
    approx.type="rauc")
plot(out6$pars[2,]/out6$pars[1,],out6$vals,type="l",xlab=expression(beta[2]/beta[1]),
    main="penalized pRAUC loss")
lines(out7$pars[2,]/out7$pars[1,],out7$vals,type="l",col=2)
out3$par
out7$par

}


}
