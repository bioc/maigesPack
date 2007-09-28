\name{robustCorr}
\alias{robustCorr}

\title{
  Calculate a robust correlation value
}

\description{
  This function is intended to calculate robust correlation values
  between pairs of rows of numerical matrix or between two numerical
  vectors.
}

\usage{
robustCorr(x, y=NULL)
}

\arguments{
  \item{x}{numerical matrix or vector. If a matrix the method calculates
    the robust correlations between all pairs of rows. If \code{x} is a
    vector, \code{y} must be specified as another vector of same length
    as \code{x} and the robust correlation between them is calculate.}
  \item{y}{optional numeric vector, must be specified if \code{x} is a vector.}
}

\value{
  If \code{x} is a matrix, the method return a list with two square
  matrices, the first one containing the robust correlation values
  between all pairs of rows from \code{x} and the second containing the
  index of the point removed from calculation. If \code{x} is a vector,
  \code{y} must be specified and the function return a list with the
  robust correlation value between them and the index of the point removed.
}

\details{
  This function calculates a robust correlation value in a procedure
  similar to the leave-one-out used for cross-validation of
  classification results. The algorithm removes one point at a time and
  calculates a usual Pearson correlation value. Them, with a vector
  \code{\bold{r}} of correlation values that has the same length as the
  columns of \code{x} (or vectors \code{x} and \code{y}), the algorithm
  decides by the \code{min}(\bold{r}) or \code{max}(\bold{r}), according
  with that one that is more distant from the median value.
}

\seealso{
  \code{\link[stats]{cor}} from package \emph{stats}.
}

\examples{
x <- runif(50, 0, 1)
y <- rbeta(50, 1, 2)
robustCorr(x, y)

z <- matrix(rnorm(100, 0, 1), 4, 25)
robustCorr(z)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
