\name{normScaleLimma}
\alias{normScaleLimma}

\title{
  Scale adjust a cDNA Microarray Object
}

\description{
  This function loads a \code{\link{maigesRaw}} or \code{\link{maiges}}
  object and scale adjust (normalise between arrays) the data using
  functions from \emph{limma} package.
}

\usage{
normScaleLimma(obj=NULL, \dots)
}

\arguments{
  \item{obj}{object of type \code{\link{maigesRaw}} or
    \code{\link{maiges}} to be  normalised.}
  \item{\dots}{additional parameters for function
    \code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}}.}
}

\details{
  This function for scale adjustment is entirely based on function
  \code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}} from
  \emph{limma} package. See the help page for this function to see how
  to set the parameters. Pay attention to the 'vsn' method of scale
  adjustment, that must be used alone.
}

\value{
  This function returns a \code{\link{maiges}} object.
}

\seealso{
  \code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}}
  from \emph{limma} package.
}

\examples{
## Loading the dataset
data(gastro)

## Doing the scale adjustment from median-absolute-value method (from
## limma)
gastro.norm = normScaleLimma(gastro.norm, method="scale")
boxplot(gastro.norm) ## To see the efect of adjustment

## To do VSN scale adjustment (from vsn package) use the command. Be
## carefull that this method adjust the variance along A values and not
## between chips!!
gastro.norm = normScaleLimma(gastro.raw2, method="vsn")
boxplot(gastro.norm) ## See the efect
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
