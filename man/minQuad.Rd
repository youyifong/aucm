\name{minQuad}
\alias{minQuad}
\alias{minQuad}
\title{minQuad}
\description{
minimizes the objective function 0.5a'Qa + b'a with respect to "a" subject to 0 <= a <= C, 
or to the constraints lower <= a <= upper, v <= Ax <= v + r.
}
\usage{       
minQuad(H,b,C = 1.0,n1=0,n2=0, 
    mem.efficient = FALSE,alpha = NULL,
    lower = NULL,upper = NULL,mat.constr = NULL, lhs.constr = NULL,rhs.constr = NULL,
    control = list(DUP = TRUE,maxit = 1e4, tol = 1e-04, 
        verbose = FALSE, ret.ws = FALSE,ret.data = FALSE,
        rank = 0,
        method = c("default","tron","loqo","exhaustive","x"),
        optim.control = list(),
        q = 2, 
        ws =  c("v","v2","greedy","rv2wg","rvwg","rv","rv2")
    )
) 
}
\arguments{
  \item{H}{A symmetric matrix whose \code{\link{typeof}} returns "double" of dimension (n x n). If mem.efficient = FALSE 
   n = n1*n2 matches the length of the vector 'b' else n = n1 + n2, see details, defaults to NULL.}  
  \item{b}{a numeric vector of length 'n' whose \code{\link{typeof}} returns "double".}
  \item{C}{a numeric variable whose \code{\link{typeof}} returns "double", defaults to 1.0 . It is the upper bound on alpha's, where the lower bound is 0.0 .}
  \item{n1,n2}{integer variables giving the specific values for n1 = #\{diseased\}, n2 = #\{non-diseased\}} subjects if mem.efficient = TRUE.
  \item{mem.efficient}{logical, if FALSE then 'H' is represented by the (n1n2 x n1n2) matrix 'Q' else by the (n1+n2 x n1+n2) matrix 'K', defaults to FALSE.}
  \item{alpha}{a length-n1n2 vector vector of initial values for "alpha", whose \code{\link{typeof}} returns "double", 
            defaults to NULL, in which case it is set to 0.5C.}
  \item{control}{a list with control parameters. See \code{\link{control.minQuad}}}
  \item{mat.constr}{m x n constraint matrix for loqo optimizer}
  \item{lhs.constr}{numeric of length 'm', the left hand side constraints for loqo optimizer}
  \item{rhs.constr}{numeric of length 'm', theleft hand side for constraints for loqo optimizer}  
  \item{lower}{numeric of length 'n', the lower bounds on primal valriables for loqo optimizer}
  \item{upper}{numeric of length 'n', the upper bounds on primal valriables for loqo optimizer}
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

The function \code{minQuad} passes its arguments by "reference" via the 
\code{\link{.C}} function call to C if DUP = FALSE to avoid copying the large matrix "Q". 
When 'H' = 'Q', 'Q' is a symmetric matrix and should have numeric type "double", 
be of type "matrix" not of "data.frame": is.matrix(.) should return "TRUE". 
We do not make an extra copy by tranposing Q but access the 'flattened' vector in C directly since 'Q' is symmetric.
When 'mem.efficient' = TRUE 'H' = K_\{n1+n2 x n1+n2\} and may be obtained by the function \code{\link{getK}}. 'K' is relevant to AUC
estimation, see \code{\link{rauc}} for more details.
The "ws" argument sets the type of strategy to select the working set. 
Denote the two sets of violators as V0 = \{1 if (a[p] > 0.0),(df[p] > 0.0), 0 ow.\},  VC = \{1 if (a[p] < C),(df[p] < 0.0), 0 ow.\}
where "df[a]" stands for the gradient of the objective function at 'a'.
"greedy" selects the extremes pairs (max,min) from the two sets (M,m) where M = \{-df[i] , a[i] < C\} and m = \{-df[j] | a[j] > 0\}.
"v" selects from violators V = V0 U VC ranked by |df|. 
"v2" selects separately from V0 and VC separately, ranked by |df|.
"rv" selects without replacement (WOR) from all violators. 
"rvwg" selects WOR from all violators V with probability ~ |df[V]|. 
"rv2wg" selects WOR from the two sets of violators V0 and VC with probability ~ |df[V]|. 
Three methods are available, "tron","hideo" and "exhaustive". Optimizer 'x' is a slightly faster implementation of the exhaustive method, whereas 
'default' is a fast implementation for q = 2 only. The "exhaustive" method should probably not be used beyond working set size 
q = 8,10 on most computers.
The 'loqo' optimizer accepts constraints of the form v <= A*x <= v + r, lower <= x <= upper. v = lhs.constr and
r = rhs.constr - lhs.constr . The entries in 'v','r' and 'A' must be finite. 
When verbose is TRUE, each DCA iteration prints one line. Delta means, for the linear kernel, max(abs((beta.new-beta.init)/beta.init)), and for the nonlinear kernel the difference in penalized RAUC loss. epsilon is the KKT criterion of minQuad.

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
\seealso{\code{\link{control.minQuad}},\code{\link{rauc}}}

\keyword{auc}
\keyword{roc}
\keyword{area}
\keyword{receiver}
