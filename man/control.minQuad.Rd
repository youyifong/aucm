\name{control.minQuad}
\alias{control.minQuad}
\title{control.minQuad}
\description{
Control function to \code{\link{minQuad}} 
}
\usage{       
control.minQuad(
    maxit = 1e4, tol = 1e-04, 
    q = 0,
    ws = c("v","v2","greedy","rv2wg","rvwg","rv","rv2"),
    method = c("default","tron","loqo","exhaustive","x"),
	optim.control = list(),
    rank = 0,
	DUP = FALSE,
	NAOK = FALSE,
	verbose = FALSE, 
	ret.ws = FALSE,
	ret.data = FALSE
)

}
\arguments{
  \item{rank}{a nonnegative integer indicating the 'rank' of the design matrix, only used by  by method 'exhaustive' or 'x'. if zero it is
              estimated by the singular value decomposition of each 'sub-matrix' associated with the sub-optimization problem
  in the decomposition method.}
  \item{method}{a character string (first letter is sufficient) indicating which quadratic optmizer to use, defaults to 'default'. See details.} 
  \item{optim.control}{a list of control parameters to methods 'tron' or 'loqo';\cr
   'tron' : list(maxfev = 1000,fatol = tol,frtol = tol,cgtol=tol,gtol=tol,fmin = -.Machine$double.xmax),\cr
   'loqo' : list(bound = 10,margin=0.05,maxiter=40,sigfig = 7,inf = 1e6)\cr   
   } 
  \item{q}{size of the working set, will be set to 2 for all methods except for method = 'tron' when it defaults to NULL. In that case
   workings set size is automatically chosen to be sqrt(#violators) at each iteration.}
  \item{ws}{a character string indicating the strategy of how to select the working set, defaults to "rv2wg", see details.}
  \item{maxit}{maximum number of iterations whose \code{\link{typeof}} returns "integer".}
  \item{tol}{tolerance for termination criterion whose \code{\link{typeof}} returns "double".}
  \item{DUP}{should arguments be passed by reference ? defaults to FALSE.}
  \item{NAOK}{should NA's,NaN's be allowed to be passed to C code (no checking) ? defaults to FALSE.}
  \item{verbose}{some output at each iteration, possible values are FALSE/TRUE or and integer if more details are wanted, defaults to FALSE.}
  \item{ret.ws}{defaults to FALSE, indicates whether to return the working set selected at each iteration.}
  \item{ret.data}{defaults to FALSE, indicates whether to return the data passed to minQuad.}
}



\value{
  A list with the following elements:
   \item{convergence}{0 if converged, 1 if maximum iteration is reached.}
   \item{alpha}{estimated vector of coefficients.}
   \item{value}{value of the objective function.}
   \item{iterations}{number of iterations until convergence or "maxit" reached.}
   \item{epsilon}{stopping rule value to be compared to "tol".}
   \item{n}{}
   \item{nSV}{#\{0 < a\}, no. of support vectors.}
   \item{nBSV}{#\{a==C\}, no. of bounded support vectors.}
   \item{nFSV}{#\{0 < alpha < C\}, no. of unbounded support vectors.} 
   \item{control}{the control argument.}
   \item{ws}{if requested, the working set selected at current iteration.}
 }
\details{
Four quadratic optimizers are available within \code{\link{minQuad}}, "default", "tron", "loqo" and "exhaustive" (optimizer 'x' is 
a slightly faster implementation of the exhaustive method). For working set size q = 2, the 'default' option is a fast implementation that
loosely minimizes the quadratic objective function, which is often sufficient to achieve convergence in the DCA-loop in \code{\link{rauc}}.
For working set size q = 2, the 'default' option minimizes the quadratic objective function by "solving" an associated equation at each data point.   
The "exhaustive" method is a brute-force method that gives an exact solution to each quadratic sub-problem in \code{\link{minQuad}} 
and should probably not be used beyond working set size q = 8,10 on most computers. Method 'tron' is a positive semidefinite 
quadratic optimizer and thus well suited for low-rank problems - for this method 'q' can be larger, ~100 or perhaps even ~1000. 
Method 'loqo' is a positive definite quadratic optimizer that accepts 'm' constraints specified by (m x n) matrix A in the form v <= A*x <= v+r with both
v and r finite.  
The default value of the size of the working set 'q' is 0. This means that if 'method' is 'tron' then 
'q' is automatically set to the sqrt(no. of violators) at the current iteration in \code{\link{minQuad}} (rounded up). Otherwise 'q'
defaults to 2 but may be set to any nonzero integer that is greater than 1. 
The "ws" argument sets the type of strategy to select the working set. 
Denote the two sets of violators as V0 = \{1 if (a[p] > 0.0),(df[p] > 0.0), 0 ow.\},  VC = \{1 if (a[p] < C),(df[p] < 0.0), 0 ow.\}
where "df[a]" stands for the gradient of the objective function at 'a'.
"greedy" selects the extremes pairs (max,min) from the two sets (M,m) where M = \{-df[i] , a[i] < C\} and m = \{-df[j] | a[j] > 0\}.
"v" selects from violators V = V0 U VC ranked by |df|. 
"v2" selects separately from V0 and VC separately, ranked by |df|.
"rv" selects without replacement (WOR) from all violators. 
"rvwg" selects WOR from all violators V with probability ~ |df[V]|. 
"rv2wg" selects WOR from the two sets of violators V0 and VC with probability ~ |df[V]|. 
}

\references{
 \emph{Combining Biomarkers Nonlinearly for Classification Using the Area Under the ROC Curve}
  Y. FONG, S. YIN, Y. HUANG
  Biometrika (2012), pp 1-28

 \emph{Newton's method for large bound-constrained optimization problems}
 Lin, C.J. and More, J.J.
 SIAM Journal on Optimization (1999), volume 9, pp 1100-1127.

  \emph{kernlab - An S4 Package for Kernel Methods in R.} 
  Alexandros Karatzoglou, Alex Smola, Kurt Hornik, Achim Zeileis 
  Journal of Statistical Software (2004) 11(9), 1-20. 
  URL http://www.jstatsoft.org/v11/i09/ 
}

\author{
  Krisztian Sebestyen \email{ksebestyen@gmail.com}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Shuxin Yin \email{}\cr
}
\seealso{\code{\link{rauc}}, \code{\link{minQuad}}}

\keyword{auc}
\keyword{roc}
\keyword{area}
\keyword{receiver}



