\name{MI}
\alias{MI}

\title{
  Calculate Mutual Information
}

\description{
  Function to calculate the mutual information of 2 random variables, or
  between all pairs of rows of a numerical matrix.
}

\usage{
MI(x, y=NULL, k=1)
}

\arguments{
  \item{x}{numerical matrix to calculate the MI between all pairs of
    rows from \code{x}. Also, \code{x} must be a numerical vector and
    \code{y} must be specified as another numerical vector of same
    lenght as \code{x} and the MI value between both them are calculated.}
  \item{y}{optional numerical vector that must be specified if \code{x}
    is a vector. Defaults to NULL.}
  \item{k}{integer specifying the number of the neighbours to be used in the
    calculation of the MI value.}
}

\value{
  If \code{x} is a matrix, the function return a square matrix with
  lenght equal to the number of rows of \code{x} with MI values between
  all pairs of rows from \code{x}. If \code{x} is a numerical vector,
  \code{y} must be specified and the function returns a positive real
  number with the MI value between the two vectors.
}

\details{
  This function implements an algorithm proposed by Kraskov et
  al. (2004) that don't use estimator of the entropy.
}

\references{
  Kraskov, A.; Stogbauer, H. and Grassberger, P. Estimating mutual
  information, \bold{Physical Review E}, 69, 066138, 2004 (\url{http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PLEEE8000069000006066138000001&idtype=cvips&gifs=yes}).
}

\examples{
x <- runif(50, 0, 1)
y <- rbeta(50, 1, 2)
MI(x, y)

z <- matrix(rnorm(100, 0, 1), 4, 25)
MI(z)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
