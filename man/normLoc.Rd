\name{normLoc}
\alias{normLoc}

\title{
  Normalise a cDNA Microarray Object
}

\description{
  This function loads a \code{\link{maigesRaw}} object and corrects for
  location bias.
}

\usage{
normLoc(obj=NULL, \dots)
}

\arguments{
  \item{obj}{object of type \code{\link{maigesRaw}} to be normalised.}
  \item{\dots}{additional parameters for function
    \code{\link[limma]{normalizeWithinArrays}} from \emph{limma} package.}
}

\details{
  This function for normalisation is entirely based on the function
  \code{\link[limma]{normalizeWithinArrays}} from \emph{limma}
  package. See their help page to known how to setup the parameters
  correctly. The parameters \code{layout} and \code{weights}, are
  automatically specified by the object \code{obj} and must not be specified.
}

\value{
  This function returns a \code{\link{maiges}} object.
}

\seealso{
  \code{\link[limma]{normalizeWithinArrays}} from \emph{limma} package.
}

\examples{
## Loading the dataset
data(gastro)

## Do the normalization by loess method and span 0.4
gastro.norm = normLoc(gastro.raw2, span=0.4, method="loess")

## Do the same normalization without background subtraction
gastro.norm = normLoc(gastro.raw2, span=0.4, method="loess", bc.method="none")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
