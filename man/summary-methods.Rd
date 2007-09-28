\name{summary}
\alias{summary.maigesPreRaw}
\alias{summary.maigesRaw}
\alias{summary.maiges}
\alias{summary.maigesANOVA}
\alias{summary.maigesDE}
\alias{summary.maigesDEcluster}
\alias{summary.maigesClass}
\alias{summary.maigesRelNetB}
\alias{summary.maigesRelNetM}
\alias{summary.maigesActMod}
\alias{summary.maigesActNet}
\alias{summary}

\title{
  Method summary for the object from this package
}

\description{
  Generic function \code{\link[base]{summary}} to display a nice summary
  of the contents from classes of objects define in this package.
}

\usage{
\method{summary}{maigesPreRaw}(object, \dots)

\method{summary}{maigesRaw}(object, \dots)

\method{summary}{maiges}(object, \dots)

\method{summary}{maigesDE}(object, \dots)

\method{summary}{maigesDEcluster}(object, \dots)

\method{summary}{maigesClass}(object, \dots)

\method{summary}{maigesRelNetB}(object, \dots)

\method{summary}{maigesRelNetM}(object, \dots)

\method{summary}{maigesActMod}(object, \dots)

\method{summary}{maigesActNet}(object, \dots)
}

\arguments{
  \item{object}{an object of any class defined in this package}
  \item{\dots}{further arguments passed to or from other methods}
}


\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\examples{
## Loading the dataset
data(gastro)

summary(gastro)
summary(gastro.raw)
summary(gastro.norm)
summary(gastro.summ)
}

\seealso{
  \code{\link[base]{summary}} in the base package.
}

\keyword{array}
