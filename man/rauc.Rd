\name{rauc}
\alias{rauc}
\alias{nsv}
\alias{nsv.rauc}
\alias{plot.rauc}
\alias{vcov.rauc}
\alias{print.minQuad}
\title{rauc}
\description{
minimizes 1 - (p)AUC plus a penalty
}
\usage{
rauc (formula, dat, s = 1,lambda=1, kernel="linear", para=NULL, start.method="rlogit", 
eta0.init=NULL,beta.init = NULL, eta.diff.init=NULL, 
maxit=50, tol=1e-5,minQuad.control = control.minQuad(),
init.alpha.from.previous = TRUE,mem.efficient = TRUE,
ret.vcov = FALSE, garbage.collection = TRUE, verbose = FALSE, ...
)
}
\arguments{
  \item{formula}{formula, e.g. y~x1+x2}
  \item{dat}{Data frame}
  \item{s}{absolute value of the slope, default to 1 - REMOVE THIS, the pair (s,lambda) is redundant}
  \item{lambda}{scale parameter in front of the penalty function, default to 1}
  \item{kernel}{See \link{getK} for more details}
  \item{para}{See \link{getK} for more details}
  \item{start.method}{a string. When kernel is linear: If "rlogit", robust logistic fit is used as beta.init. If "1", a vector of 1 is used as beta.init. If "0", a vector of 0 is used as beta.init.}
  \item{eta0.init}{a vector of the same length as the number of rows in dat}
  \item{beta.init}{a vector of length equal to no. of covariates (without intercept) of initial values for linear kernel.}
  \item{eta.diff.init}{a vector of the same length as the number of rows in dat}
  \item{maxit}{maximum number of iterations in the DCA algorithm}
  \item{tol}{absolute tolerance in RAUC if kernel is not linear, relative tolerance in coefficients if kernel is linear.}
  \item{minQuad.control}{control parameters passed to method minQuad, please see \code{\link{minQuad}}.}  
  \item{init.alpha.from.previous}{defaults to TRUE, if TRUE then after the first iteration \code{\link{minQuad}} 
               receives as the initial "alpha" the estimate of "alpha" from the previous iteration in dca algorithm.}
   
  \item{mem.efficient}{if TRUE, the small matrix 'K' instead of 'Q' is used in computations, defaults to TRUE.}
  \item{ret.vcov}{logical, whether to return an estimate of the covariance matrix of 'beta' for normal or logistic sigmoid functions.}
  \item{garbage.collection}{logical, whether to call \code{\link{gc}} at end of each DCA iteration}
  \item{verbose}{prints information at each iteration, defaults to FALSE}
  \item{\dots}{for debugging purposes only}
}

\value{
  A list with the following elements:
   \item{convergence}{0 if converged, 1 if maximum iteration is reached.}
   \item{value}{value of the objective function.}
   \item{iterations}{number of iterations until convergence or 'maxit' reached.}
 }

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{ksebestyen@gmail.com}\cr
}

\examples{

\dontrun{

# options(path.svml = 'D:/downloaded_scientific_programs/svmlight') 
# options(path.svml ='~/bin/svmlight')

###########################################################
# a linear example

dat = sim.dat.1(n=200,seed=1)

# convergence takes long, to pass CRAN check, set maxit=1 

fit1 = rauc (y~x1+x2, dat, lambda=2, kernel="linear", maxit=2)
#fit2 = rauc.linear (y~x1+x2, dat, lambda=2, verbose=TRUE)
#aux2=fit2$X \%*\% fit2$coefficients
#all(fit1$linear.combination-aux2<1e-2)
fit1$train.auc # 0.7206015


fit3 = rauc (y~x1+x2, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
fit3$train.auc # 0.7773434



fit4 = svml (y~x1+x2, dat, kernel="r", fitted=FALSE, cost=1e4) 
fast.auc(predict(fit4, dat)$posterior[,1], dat$y) # 0.7921805
tune.svml(y~x1+x2, dat, kernel="r")
#        1        10       100      1000     10000     1e+05
#0.7027569 0.7254135 0.7517794 0.7653133 0.7921805 0.6674687

# glm derived score for comparision
fit.glm=glm(y~x1+x2, dat, family="binomial")
fast.auc(fit1$X \%*\% fit.glm$coef[-1], fit1$y) # 

# add outliers
dat = sim.dat.1(n=200,seed=1, add.outliers=TRUE)

fit3 = rauc (y~x1+x2, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
fit3$train.auc # 0.7066667

fit4 = svml (y~x1+x2, dat, kernel="r", fitted=FALSE, cost=1e4) 
fast.auc(predict(fit4, dat)$posterior[,1], dat$y) # 0.6910101
tune.svml(y~x1+x2, dat, kernel="r")
#        1        10       100      1000     10000     1e+05 
#0.6485859 0.6705051 0.6722222 0.6767677 0.6910101 0.5007071



###########################################################
# a nonlinear example

dat=skin.orange (n=100,seed=1,noise=FALSE)
dim(dat)

# nonlinear kernel fit
fit1 = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
# glm fit
fit.glm=glm(y~x1+x2+x3+x4, dat, family="binomial")
# linear kernel fit
fit2 = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="linear", start.method = "rlogit", verbose=TRUE)

# training data prediction
fast.auc(fit1$linear.combination, fit1$y)
fast.auc(fit1$X \%*\% fit.glm$coef[-1], fit1$y)
fast.auc(fit2$linear.combination, fit2$y)

# test data prediction
newdata=skin.orange (n=1000,seed=2,noise=FALSE)
fast.auc(predict(fit1, newdata), newdata$y)
fast.auc(as.matrix(subset(newdata, select=c(x1,x2,x3,x4))) \%*\% fit.glm$coef[-1], newdata$y)
fast.auc(predict(fit2, newdata), newdata$y)



###### IMPROVEMENTS ####################################################

 
## rank = 2 problem 
dat = sim.dat.1(n=300,seed=1,add.outliers = TRUE,std.dev = 1.0);fm = y~x1+x2

## linear kernel and random working set selection - low rank (2) problem
## setting initial alpha (to be passed to minQuad at each iteration in dca-loop) 
# to estimate from previous dca() iteration 
## size of working set is automatically set
set.seed(100) 
fit.lin = rauc (fm, dat,lambda=.1,kernel="linear",
verbose=TRUE,maxit = 100,tol = 1e-5,
init.alpha.from.previous = TRUE,mem.efficient = TRUE,
minQuad.control = control.minQuad(
                            verbose = 1,maxit = 1e6,tol = 1e-4,
                            method = "tron",                            
                            working.set= "rv2wg")
)

## 'rbf' kernel and random working set selection
## low rank mapped to possibly infinite rank problem try larger working set 'q' set.seed(100) 
## size of working set is set to q = 100
fit.rbf = rauc (fm, dat,lambda=.1,kernel="rbf",para = 1, verbose=TRUE,maxit = 100,tol = 1e-5,
init.alpha.from.previous = TRUE,mem.efficient = TRUE,
minQuad.control = control.minQuad(
                            verbose = 1,maxit = 1e6,tol = 1e-4,
                            q = 100,
                            method = "tron",                            
                            working.set= "rv2wg")
)

}
 
}

\keyword{auc}
\keyword{roc}
\keyword{area}
\keyword{receiver}
