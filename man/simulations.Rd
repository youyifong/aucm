\name{simulations}
\alias{simulations}
\alias{sim.dat.1}
\alias{sim.MH1}
\alias{sim.MH11}
\alias{sim.MH2}
\alias{sim.NL1}
\alias{sim.YH1}
\alias{sim.YH2}
\alias{sim.d20c}
\alias{sim.disc}
\alias{sim.easy}
\alias{sim.p1}
\alias{sim.ring}
\alias{skin.orange}
\alias{sim.NL}
\alias{sim.mix}
\alias{sim.pepe.2}
\title{Simulate datasets}
\description{
\code{sim.dat.1} simulates a dataset with two covariates to reproduce Pepe Figure 1.
\code{skin.orange} simulates a skin of orange dataset as in Hastie et al. 
}
\usage{
sim.dat.1(n, seed, add.outliers=FALSE, std.dev = 0.2)
}
\arguments{
  \item{n}{sample size}
  \item{seed}{seed for random number generator}
  \item{add.outliers}{boolean. If TRUE, 10\% of data are replaced by a contaminating distribution}
  \item{std.dev}{standard deviation in data generating process}
}

\value{
  A data frame with n rows, and 4 columns: y, x1, x2, and eta, where eta is the linear combination X*beta. 
}

\author{
  Shuxin Yin \email{}\cr
  Youyi Fong \email{youyifong@gmail.com}\cr
  Krisztian Sebestyen \email{}\cr
}

\examples{

dat = sim.dat.1(n=100,seed=1)
nrow(dat)

dat = sim.dat.1(n=100,seed=1,add.outliers=TRUE)
nrow(dat)

}
