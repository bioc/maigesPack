\name{normOLIN}
\alias{normOLIN}

\title{
  Normalise a cDNA Microarray Object
}

\description{
  This function loads a \code{\link{maigesRaw}} object and normalise it
  using methods OLIN or OSLIM.
}

\usage{
normOLIN(obj=NULL, \dots)
}

\arguments{
  \item{obj}{object of type \code{\link{maigesRaw}} to be normalised.}
  \item{\dots}{additional parameters for function \code{\link[OLIN]{olin}}.}
}

\details{
  This function for normalisation is entirely based on function
  \code{\link[OLIN]{olin}} from \emph{OLIN} package. This function
  implements the methods OLIN and OSLIM, proposed by Futschich and
  Crompton (2004). See help page for this function to discover how to
  set the parameters.
}

\references{
  Futschik, M. and Crompton, T. Model selection and efficiency testing
  for normalization of cDNA microarray data, \bold{Genome Biology}, 5,
  R60, 2004 (\url{http://genomebiology.com/2004/5/8/R60}).
}

\value{
  This function returns a \code{\link{maiges}} object.
}

\examples{
## Loading the dataset
data(gastro)

## Doing the OLIN normalization. Pay attention because, this methods are
## very time consuming!
\dontrun{
gastro.norm = normOLIN(gastro.raw2) ## Without scale adjustment

gastro.norm = normOLIN(gastro.raw2, OSLIN=TRUE) ## With scale adj
}

}

\seealso{
 \code{\link[OLIN]{olin}} from \emph{OLIN} package.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
