\name{[-method}
\docType{methods}

\alias{[,maiges-method}
\alias{[,maigesANOVA-method}
\alias{[,maigesRaw-method}
\alias{[,maigesPreRaw-method}
\alias{[}

\title{Sub-setting methods for maiges objects}

\description{
  Sub-setting methods were defined for the classes presented
  in this package, \code{\link{maiges}}, \code{\link{maigesANOVA}}, 
  \code{\link{maigesRaw}} and \code{\link{maigesPreRaw}}. These methods
  create instances of the given class, for a subset of spots and/or
  arrays in a batch.
}

\section{Methods}{
  \describe{
    \item{x = ANY}{generic method.}
    \item{x = maiges}{\code{x[i, j]} extract object of class
      \code{\link{maiges}} for spots with indexes i and samples with
      indexes j.}
    \item{x = maigesANOVA}{\code{x[i, j]} extract object of class
      \code{\link{maigesANOVA}} for spots with indexes i and samples
      with indexes j.}
    \item{x = maigesRaw}{\code{x[i, j]} extract object of class
      \code{\link{maigesRaw}} for spots with indexes i and arrays with
      indexes j.}
    \item{x = maigesPreRaw}{\code{x[i, j]} extract object of class
      \code{\link{maigesPreRaw}} for spots with indexes i and arrays with
      indexes j.}
  }
}

\seealso{
  \code{\link{maiges}}, \code{\link{maigesANOVA}},
  \code{\link{maigesRaw}} and \code{\link{maigesPreRaw}}.
}

\examples{
## Loading the dataset
data(gastro)

gastro[1:10,]
gastro[1,1]

gastro.raw[rep(TRUE, 15),]

gastro.norm[c(1,4,6), c(10, 18)]
}

\keyword{methods}
