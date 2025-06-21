\name{getQ}
\alias{getQ}
\alias{getQ}
\title{getQ}
\description{
\code{getQ} calculates Q or Q.pred matrix depending on the value of do.pred.
}
\usage{
getQ (K,n1,n2,call.C=TRUE,do.pred=FALSE)
}
\arguments{
  \item{K}{kernel matrix of dimension (n1+n2) by (n1+n2) or n.pred by (n1+n2). The (n1+n2) observations must be ordered case followed by non-case}
  \item{n1}{number of cases}
  \item{n2}{number of non-cases}
  \item{call.C}{boolean. If TRUE, make .C call, otherwise compute Q in R.}
  \item{do.pred}{boolean. If TRUE, K is a n.pred by (n1+n2) matrix; otherwise, it is a (n1+n2) by (n1+n2) matrix. }
}

\value{
  A n1*n2 by n1*n2 matrix if do.pred is FALSE, or n.pred by n1*n2 matrix if do.pred is TRUE
}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}

\examples{

dat = sim.dat.1(n=100,seed=1)
dat=rbind(subset(dat, y==1), subset(dat, y==0))
X = as.matrix(subset(dat, select=c(x1,x2)))
n1=sum(dat$y)
n2=sum(1-dat$y)

K = kyotil::getK(X,"linear", 1)
Q1 = getQ(K,n1=n1,n2=n2,call.C=FALSE)
Q2 = getQ(K,n1=n1,n2=n2,call.C=TRUE)
all(Q2-Q1<1e-6)

# compare to a direct computation
X.diff=get.X.diff(X[1:n1,], X[1:n2+n1,])
Q3 = tcrossprod(X.diff, X.diff)
all(Q3-Q1<1e-6)

# two printouts of Q2 should not be different
Q2[1:3,1:3]
K = kyotil::getK(X,"rbf", 1)
Q4 = getQ(K,n1=n1,n2=n2,call.C=TRUE)
Q2[1:3,1:3]
Q4[1:3,1:3]

K = kyotil::getK(X[1:10,],"linear", 1, X2=X)
Q5 = getQ(K,n1=n1,n2=n2,call.C=FALSE,do.pred=TRUE)
Q6 = getQ(K,n1=n1,n2=n2,call.C=TRUE,do.pred=TRUE)
dim(Q5)
dim(Q6)
all(Q5-Q6<1e-6)



}
