\name{show-method}
\docType{methods}
\alias{show,maigesPreRaw-method}
\alias{show,maigesRaw-method}
\alias{show,maiges-method}
\alias{show,maigesANOVA-method}
\alias{show,maigesDE-method}
\alias{show,maigesDEcluster-method}
\alias{show,maigesClass-method}
\alias{show,maigesRelNetB-method}
\alias{show,maigesRelNetM-method}
\alias{show,maigesActMod-method}
\alias{show,maigesActNet-method}
\alias{show}

\title{
  Show a nice visualisation of the objects defined in this package
}

\usage{
\S4method{show}{maigesPreRaw}(object)

\S4method{show}{maigesRaw}(object)

\S4method{show}{maiges}(object)

\S4method{show}{maigesANOVA}(object)

\S4method{show}{maigesDE}(object)

\S4method{show}{maigesDEcluster}(object)

\S4method{show}{maigesClass}(object)

\S4method{show}{maigesRelNetM}(object)

\S4method{show}{maigesRelNetB}(object)

\S4method{show}{maigesActMod}(object)

\S4method{show}{maigesActNet}(object)
}

\arguments{
  \item{object}{object of any class defined in this package, to be
    printed in a nice way}
}

\description{
  S4 method to show the objects defined in this package without list all
  the slots of the object. 
}

\note{It is possible to use Package \emph{convert} to convert objects between
  classes defined in packages \emph{limma}, \emph{marray} and \emph{Biobase}.
}

\examples{
## Loading the dataset
data(gastro)

gastro
gastro.raw
gastro.norm
gastro.summ
}

\seealso{
\code{\link[methods]{show}} in the \emph{methods} package.
}

\keyword{methods}


