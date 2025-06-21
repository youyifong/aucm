\name{get.X.diff}
\alias{get.X.diff}
\alias{get.X.diff.default}
\alias{get.X.diff.formula}
\title{get.X.diff}
\description{
computes X.diff matrix
}
\usage{
get.X.diff (x1,...)
\S3method{get.X.diff}{default}(x1,x2,...)
\S3method{get.X.diff}{formula}(formula, data,...)
}
\arguments{
  \item{x1}{data matrix from the case group, dimension n1 x d}
  \item{x2}{data matrix from the non-case group, dimension n2 x d}
  \item{formula}{a formula}
  \item{data}{a data frame}  
  \item{...}{arguments passed 'to' or 'from' methods}
}

\value{
  A (n1*n2) x d matrix
}

\details{
In get.X.diff.formula, x is the case predictors and x2 control.
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


}
