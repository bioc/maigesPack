\name{calcW}
\alias{calcW.default}
\alias{calcW.maigesRaw}
\alias{calcW}

\title{
  Method calcW to calculate W values
}

\description{
  Generic function \code{\link{calcW}} to calculate W values from
  classes of  microarray data objects defined in this package.
}

\usage{
calcW(object, \dots)

\method{calcW}{default}(object, \dots)

\method{calcW}{maigesRaw}(object, bkgSub="subtract", \dots)

}

\arguments{
  \item{object}{object of any class. But only some methods are defined
    in this moment.}
  \item{bkgSub}{character string indicating the type of background
    subtraction. May be 'none', 'subtract', 'half', 'minimum',
    'movingmin', 'edwards', 'normexp' or 'rma'. Uses limma and defaults
    to 'subtract'.}
  \item{\dots}{additional parameters for \code{\link{calcW}} method.}
}

\details{
  This method receive an object (at moment of class
  \code{\link{maiges}}, \code{\link{maigesRaw}} or
  \code{\link{maigesANOVA}}) and returns the matrix of W values. For
  objects of class \code{\link{maigesRaw}} it uses the function
  \code{\link[limma]{backgroundcorrect}} from \emph{limma} package to do
  background correction before the calculation of W values.
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\seealso{
  \code{\link[limma]{backgroundcorrect}} in the limma package and \code{\link{calcA}}.
}

\examples{
## Loading the dataset
data(gastro)

calcW(gastro.raw)
calcW(gastro.norm)
}

\keyword{array}
