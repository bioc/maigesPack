\name{bootstrapT}
\alias{bootstrapT}

\title{
  Calculate bootstrap p-values for t statistics
}

\description{
  This function takes a numerical matrix and column indexes for two
  groups to calculate bootstrapped (by re-sampling) p-values comparing the
  equality of means from the two groups.
}

\usage{
bootstrapT(x, k=20000, obs1, obs2, \dots)
}

\arguments{
  \item{x}{numerical matrix to be bootstrapped. The t statistics is
    calculated by row using the column indexes given by \code{obs1} and
    \code{obs2} for the two groups tested.}
  \item{k}{number of bootstrap re-samplings to be done. Defaults to 20000.}
  \item{obs1}{logical or numerical column indexes of the first group.}
  \item{obs2}{logical or numerical column indexes of the second group.}
  \item{\dots}{additional parameters for \code{\link[stats]{t.test}}
    function from package \emph{stats}.}
}

\value{
  The result of this function is a numerical matrix with number of rows
  given by the rows of the argument \code{x} and 3 columns. The first
  column contain the difference of means between the two groups, the
  second one contain the original t statistic and the last one gives the
  bootstrapped p-values, for all rows of the matrix \code{x}.
}

\seealso{
  \code{\link[stats]{t.test}} from package \emph{stats}.
}

\examples{
z <- matrix(rnorm(100, 0, 1), 4, 25)
bootstrapT(z, k=100, obs1=1:14, obs2=15:25)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
