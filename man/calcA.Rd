\name{calcA}
\alias{calcA.default}
\alias{calcA.maigesRaw}
\alias{calcA}

\title{
  Method calcA to calculate A values
}

\description{
  Generic function \code{\link{calcA}} to calculate A values from
  classes of  microarray data objects defined in this package.
}

\usage{
calcA(object, \dots)

\method{calcA}{default}(object, \dots)

\method{calcA}{maigesRaw}(object, bkgSub="subtract", \dots)

}

\arguments{
  \item{object}{object of any class. But only some methods are defined
    in this moment.}
  \item{bkgSub}{character string indicating the type of background
    subtraction. May be 'none', 'subtract', 'half', 'minimum',
    'movingmin', 'edwards', 'normexp' or 'rma'. Uses limma and defaults
    to 'subtract'.}
  \item{\dots}{additional parameters for \code{\link{calcA}} method.}
}

\details{
  This method receive an object (at moment of class
  \code{\link{maiges}}, \code{\link{maigesRaw}} or
  \code{\link{maigesANOVA}}) and returns the matrix of A values. For
  objects of class \code{\link{maigesRaw}} it uses the function
  \code{\link[limma]{backgroundcorrect}} from \emph{limma} package to do
  background correction before the calculation of A values.
}

\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\examples{
## Loading the dataset
data(gastro)

calcA(gastro.raw)
calcA(gastro.norm)
}

\seealso{
  \code{\link[limma]{backgroundcorrect}} in the limma package and \code{\link{calcW}}.
}

\keyword{array}
