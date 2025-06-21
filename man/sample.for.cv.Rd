\name{sample.for.cv}
\alias{sample.for.cv}
\alias{sample.for.cv}
\title{sample.for.cv}
\description{
\code{sample.for.cv} does sampling for cross validation
}
\usage{
sample.for.cv (dat, v, seed)
}
\arguments{
  \item{dat}{a data frame. One of the columns must be named y and y should be 0/1 with 1 for case and 0 for control}
  \item{v}{v-fold cross validation}
  \item{seed}{seed for random number generators}
}

\value{
  A list of two vector of integers: train and test, which refer to the rows of dat
}
\details{
    case and controls are sampled separately
}

\author{
  Youyi Fong \email{youyifong@gmail.com}\cr
}

