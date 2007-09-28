\name{coerce-method}
\docType{methods}
\alias{coerce,maiges,marrayNorm-method}
\alias{coerce,marrayNorm,maiges-method}
\alias{coerce,maiges,MAList-method}
\alias{coerce,MAList,maiges-method}
\alias{coerce,maigesRaw,marrayRaw-method}
\alias{coerce,marrayRaw,maigesRaw-method}
\alias{coerce,maigesRaw,RGList-method}
\alias{coerce,RGList,maigesRaw-method}
\alias{as,maiges,marrayNorm-method}
\alias{as,marrayNorm,maiges-method}
\alias{as,maiges,MAList-method}
\alias{as,MAList,maiges-method}
\alias{as,maigesRaw,marrayRaw-method}
\alias{as,marrayRaw,maigesRaw-method}
\alias{as,maigesRaw,RGList-method}
\alias{as,RGList,maigesRaw-method}
\alias{as}

\title{
  Coerce a maiges object to classes defined by packages limma and marray
}

\usage{
\S4method{as}{maiges,marrayNorm}(from, to)

\S4method{as}{marrayNorm,maiges}(from, to)

\S4method{as}{maiges,MAList}(from, to)

\S4method{as}{MAList,maiges}(from, to)

\S4method{as}{maigesRaw,marrayRaw}(from, to)

\S4method{as}{marrayRaw,maigesRaw}(from, to)

\S4method{as}{maigesRaw,RGList}(from, to)

\S4method{as}{RGList,maigesRaw}(from, to)
}

\arguments{
  \item{from}{object of class \code{\link{maiges}},
    \code{\link{maigesRaw}}, \code{\link[limma:malist]{MAList}},
    \code{\link[limma:rglist]{RGList}}, \code{\link[marray:marrayNorm-class]{marrayNorm}} or
    \code{\link[marray:marrayRaw-class]{marrayRaw}}.}
  \item{to}{character string specifying the class of object to which the
    object \code{from} will be coerced}
}

\description{
  Coercing methods were defined to convert maiges objects of classes
  \code{\link{maiges}} and \code{\link{maigesRaw}} into objects of classes
  \code{\link[marray:marrayNorm-class]{marrayNorm}} and \code{\link[marray:marrayRaw-class]{marrayRaw}}
  from package \emph{marray} or classes \code{\link[limma:malist]{MAList}} and
  \code{\link[limma:rglist]{RGList}} from package \emph{limma} and vice-versa. 
}

\details{
  When converting from objects of classes \code{\link[limma:rglist]{RGList}} or
  \code{\link[marray:marrayRaw-class]{marrayRaw}} to class \code{\link{maigesRaw}}, the
  slot \code{Sf} and \code{Sb} will always receive the channel 2 (red)
  values, and \code{Rf} and \code{Rb} will always receive channel 1
  (green) values. For the normalized ojects, the slot \code{W} will be
  equivalent to the M values.

  When converting from objects of classes \code{\link{maigesRaw}} or
  \code{\link{maiges}} to classes \code{\link[limma:rglist]{RGList}},
  \code{\link[marray:marrayRaw-class]{marrayRaw}},
  \code{\link[limma:malist]{MAList}} or
  \code{\link[marray:marrayNorm-class]{marrayNorm}} the correct values
  of the two channels or of M values are calculated.
}

\section{Methods}{
  \describe{
    \item{from = maiges, to = marrayNorm}{convert an object of class
      \code{\link{maiges}} into an object of class \code{\link[marray:marrayNorm-class]{marrayNorm}}.} 
    \item{from = marrayNorm, to = maiges}{convert an object of class
      \code{\link[marray:marrayNorm-class]{marrayNorm}} into an object of class \code{\link{maiges}}.} 
    \item{from = maiges, to = MAList}{convert an object of class
      \code{\link{maiges}} into an object of class \code{\link[limma:malist]{MAList}}.} 
    \item{from = MAList, to = maiges}{convert an object of class
      \code{\link[limma:malist]{MAList}} into an object of class \code{\link{maiges}}.} 
    \item{from = maigesRaw, to = marrayRaw}{convert an object of class
      \code{\link{maigesRaw}} into an object of class \code{\link[marray:marrayRaw-class]{marrayRaw}}.} 
    \item{from = marrayRaw, to = maigesRaw}{convert an object of class
      \code{\link[marray:marrayRaw-class]{marrayRaw}} into an object of class \code{\link{maigesRaw}}.} 
    \item{from = maigesRaw, to = RGList}{convert an object of class
      \code{\link{maigesRaw}} into an object of class \code{\link[limma:rglist]{RGList}}.} 
    \item{from = RGList, to = maigesRaw}{convert an object of class
      \code{\link[limma:rglist]{RGList}} into an object of class \code{\link{maigesRaw}}.} 
  }
}

\note{It is possible to use Package \emph{convert} to convert objects between
  classes defined in packages \emph{limma}, \emph{marray} and \emph{Biobase}.
}

\examples{
## Loading the dataset
data(gastro)

## Converting a maigesRaw class object into marrayRaw object
as(gastro.raw, "marrayRaw")

## Converting a maigesRaw class object into RGList
as(gastro.raw, "RGList")

## Converting a maiges class object into marrayNorm object
as(gastro.norm, "marrayNorm")

## Converting a maiges class object into MAList object
as(gastro.summ, "MAList")
}

\seealso{
\code{\link[methods]{as}} in the \emph{methods} package.
}

\keyword{methods}


