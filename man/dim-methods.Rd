\name{dim}
\alias{dim.maiges}
\alias{dim.maigesRaw}
\alias{dim.maigesPreRaw}
\alias{dim.maigesANOVA}

\title{
  Retrieve the dimension of microarray objects
}

\description{
  Generic function \code{\link[base]{dim}} to retrieve the number of
  rows (genes) and columns (arrays) for classes \code{\link{maigesPreRaw}},
  \code{\link{maigesRaw}}, \code{\link{maiges}} and \code{\link{maigesANOVA}}.
}

\usage{
\method{dim}{maigesPreRaw}(x)

\method{dim}{maigesRaw}(x)

\method{dim}{maiges}(x)

\method{dim}{maigesANOVA}(x)
}

\arguments{
  \item{x}{an object of class \code{\link{maigesPreRaw}},
    \code{\link{maigesRaw}}, \code{\link{maiges}} or \code{\link{maigesANOVA}}}
}

\details{
  This is a function to retrieve the dimensions of the dataset - number of
  genes (rows) and samples (columns) studied. Note that the commands
  \code{\link[base]{nrow}(x)}, \code{\link[base:nrow]{ncol}(x)} and
  commands related with matrix dimensions on also work.
}

\value{
  Vector with 2 numbers, the first one is the number of rows (genes) and the
  second one the number of columns (samples).
}

\author{
  modified from the \code{marray} package
}

\examples{
## Loading the dataset
data(gastro)

dim(gastro)
dim(gastro.raw)
dim(gastro.norm)
dim(gastro.summ)
}

\seealso{
  \code{\link[base]{dim}} in the base package.
}

\keyword{array}
