\name{normScaleMarray}
\alias{normScaleMarray}

\title{
  Scale adjust a cDNA Microarray Object
}

\description{
  This function loads a \code{\link{maigesRaw}} or \code{\link{maiges}}
  object and scale adjust (normalise between arrays) the data using
  functions from \emph{marray} package.
}

\usage{
normScaleMarray(obj=NULL, \dots)
}

\arguments{
  \item{obj}{object of type \code{\link{maigesRaw}} or
    \code{\link{maiges}} to be  normalised.}
  \item{\dots}{additional parameters for function
    \code{\link[marray]{maNormScale}}.}
}

\details{
  This function for scale adjustment is entirely based on function
  \code{\link[marray]{maNormScale}} from \emph{marray} package. See the
  help page for this function to see how to set the parameter. Pay
  attention to the \code{subset} argument that is fixed directly from
  the \code{UseSpots} and \code{BadSpots} from \code{obj} object, and
  must not be specified in the additional arguments.

  The functionality of the scale adjustment function from \emph{marray}
  package was added because it uses an estimator of MAD different from
  that one used in \emph{limma} package. Also, using 
  \code{\link[marray]{maNormScale}} function it is possible to do print
  tip scale adjustment.
}

\value{
  This function returns a \code{\link{maiges}} object.
}

\seealso{
  \code{\link[marray]{maNormScale}} from \emph{marray} package.
}

\examples{
## Loading the dataset
data(gastro)

## Doing global MAD scale adjustment
gastro.norm = normScaleMarray(gastro.norm, norm="globalMAD")
boxplot(gastro.norm) ## To see the effect of MAD adjustment

## For print tip MAD use the following command
\dontrun{
gastro.norm = normScaleMarray(gastro.norm, norm="printTipMAD")
}
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
