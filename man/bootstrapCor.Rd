\name{bootstrapCor}
\alias{bootstrapCor}

\title{
  Calculate bootstrap p-values for correlation measures
}

\description{
  This function takes a numerical matrix (or two vectors) and calculates
  bootstrapped (by permutation) p-values to test if the correlation value
  is equal to zero. If the first argument is a matrix, the p-values are
  calculated between all pairs of rows of the matrix.
}

\usage{
bootstrapCor(x, y=NULL, bRep, type="Rpearson", ret="p-value",
             alternative="two.sided")
}

\arguments{
  \item{x}{numerical matrix or vector to be analysed. If a vector, the
    argument \code{y} must be informed.}
  \item{y}{numerical vector. Must be informed if \code{x} is a
    vector. If \code{x} is a matrix, this argument is ignored. Defaults
    to NULL.}
  \item{bRep}{number of permutation to be done in the test.}
  \item{type}{character string specifying the type of correlation
    statistic to be used. Possible values are 'Rpearson', 'pearson',
    'spearman' or 'kendall'.}
  \item{ret}{character string with the value to return. Must be
    'p-value' (default) for the usual p-value or 'max', to return the
    maximum absolute correlation value obtained by the permutation.}
  \item{alternative}{character specifying the type of test to do, must be
    'two.sided' (default), 'less' or 'greater'.}
}

\value{
  The result of this function is a square matrix (length equal to the
  number of rows of \code{x}) if \code{x} is a matrix or a numerical
  value if \code{x} and \code{y} are vectors. The result is the p-values
  or maximum correlation values calculated by permutation tests.
}

\details{
  Pearson, spearman and kendall types of correlation values are
  calculated by \code{\link[stats]{cor}} function from package
  stats. The method Rpearson was developed in this package and is a
  generalisation of the \emph{jackniffe} correlation proposed by Heyer
  et al. (1999), it
  is calculated using the function \code{\link{robustCorr}}.
}

\references{
  Heyer, L.J.; Kruglyak, S. and Yooseph, S. Exploring expression data:
  identification and analysis of coexpressed genes, \bold{Genome
    Research}, 9, 1106-1115, 1999 (\url{http://www.genome.org/cgi/content/full/9/11/1106})
}

\seealso{
  \code{\link[stats]{cor}}, \code{\link{robustCorr}}
}

\examples{
x <- runif(50, 0, 1)
y <- rbeta(50, 1, 2)
bootstrapCor(x, y, bRep=100)

z <- matrix(rnorm(100, 0, 1), 4, 25)
bootstrapCor(z, bRep=100)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
