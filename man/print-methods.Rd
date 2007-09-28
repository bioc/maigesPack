\name{print}
\alias{print.maigesPreRaw}
\alias{print.maigesRaw}
\alias{print.maiges}
\alias{print.maigesANOVA}
\alias{print.maigesDE}
\alias{print.maigesDEcluster}
\alias{print.maigesClass}
\alias{print.maigesRelNetB}
\alias{print.maigesRelNetM}
\alias{print.maigesActMod}
\alias{print.maigesActNet}
\alias{print}

\title{
  Method to print a nice visualisation of the objects defined in this package
}

\description{
  Generic function \code{\link[base]{print}} to display a nice (and
  simple) visualisation of the objects define in this package.
}

\usage{
\method{print}{maigesPreRaw}(x, \dots)

\method{print}{maigesRaw}(x, \dots)

\method{print}{maiges}(x, \dots)

\method{print}{maigesDE}(x, \dots)

\method{print}{maigesDEcluster}(x, \dots)

\method{print}{maigesClass}(x, \dots)

\method{print}{maigesRelNetB}(x, \dots)

\method{print}{maigesRelNetM}(x, \dots)

\method{print}{maigesActMod}(x, \dots)

\method{print}{maigesActNet}(x, \dots)
}

\arguments{
  \item{x}{an object of any class defined in this package}
  \item{\dots}{further arguments passed to or from other methods}
}


\author{
  Gustavo H. Esteves <gesteves@vision.ime.usp.br>
}

\seealso{
  \code{\link[base]{print}} in the base package.
}

\examples{
## Loading the dataset
data(gastro)

print(gastro)
print(gastro.raw)
print(gastro.norm)
print(gastro.summ)
}

\keyword{array}
