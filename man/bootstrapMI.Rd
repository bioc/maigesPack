\name{bootstrapMI}
\alias{bootstrapMI}

\title{
  Calculate bootstrap p-values for mutual information (MI) measures
}

\description{
  This function takes a numerical matrix (or two vectors) and calculates
  bootstrapped (by permutation) p-values to test if the mutual
  information value is equal to zero. If the first argument is a matrix,
  the p-values are calculated between all pairs of rows of the matrix.
}

\usage{
bootstrapMI(x, y=NULL, bRep, ret="p-value")
}

\arguments{
  \item{x}{numerical matrix or vector to be analysed. If a vector, the
    argument \code{y} must be informed.}
  \item{y}{numerical vector. Must be informed if \code{x} is a
    vector. If \code{x} is a matrix, this argument is ignored. Defaults
    to NULL.}
  \item{bRep}{number of permutation to be done in the test.}
  \item{ret}{character string with the value to return. Must be
    'p-value' (default) for the usual p-value or 'max', to return the
    maximum absolute correlation value obtained by the permutation.}
}

\value{
  The result of this function is a square matrix (length equal to the
  number of rows of \code{x}) if \code{x} is a matrix or a numerical
  value if \code{x} and \code{y} are vectors. The result is the p-values
  or maximum MI values calculated by permutation tests.
}

\details{
  The method implemented in this function is proposed by Butte and
  Kohane (2000). The MI value is calculated using the function \code{\link{MI}}.
}

\references{
  Butte, A.J. and Kohane, I.S. Mutual information relevance networks:
  functional genomic clustering using pairwise entropy measurements. In
  Pacific Symposium on Biocomputing, 5, 415-426, 2000
  (\url{http://psb.stanford.edu/psb-online/proceedings/psb00/})
}  

\seealso{
  \code{\link{MI}}
}

\examples{
x <- runif(50, 0, 1)
y <- rbeta(50, 1, 2)
bootstrapMI(x, y, bRep=100)

z <- matrix(rnorm(100, 0, 1), 4, 25)
bootstrapMI(z, bRep=100)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
