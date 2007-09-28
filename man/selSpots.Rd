\name{selSpots}
\alias{selSpots}

\title{
  Select spots to use in normalisation
}

\description{
  Function to select spots to be used in microarray normalisation.
}

\usage{
selSpots(obj=NULL, sigNoise=1, rmFlag=NULL, gLabelsID=c("Name"),
         remove=list(c("BLANK","DAP","LYS","PHE","Q_GENE","THR","TRP")),
         badSpots=NULL, badLabel=NULL)
}

\arguments{
  \item{obj}{object of class \code{\link{maigesRaw}}.}
  \item{sigNoise}{positive real number indicating the cutoff to remove
    spots with signal-to-noise ratio below it.}
  \item{rmFlag}{vector of flag symbols to be removed (for
    normalisation). These flags are stored in the slot \code{Flags} from
    \code{\link{maigesRaw}} class.}
  \item{gLabelsID}{character vector indicating the gene labels to be
    searched to exclude that ones specified in \code{remove} argument.}
  \item{remove}{list of same length as \code{GlabelsID} containing character
    vector indicating the symbols of spots to be removed, according to
    the \code{GlabelsID} argument.}
  \item{badSpots}{index of bad spots (numeric or logical) identifying
    bad spots. May be the gene labels, with label ID specified by the
    argument \code{badLabel}.}
  \item{badLabel}{character string specifying the gene label ID for
    remove \code{badSpots}.}
}
 
\details{
  This function takes the object of class \code{\link{maigesRaw}} and
  actualise the slot \code{UseSpots} according with the arguments passed
  to the function. This slot is read by the normalisation functions
  \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}} and
  \code{\link{normScaleMarray}} to use only the spots that passed the
  criteria specified here.
}

\value{
  This function returns another object of class \code{\link{maigesRaw}}
  with the \code{UseSpots} slot actualised.
}

\examples{
## Loading the dataset
data(gastro)

## Filtering all spots with signal2noise ratio (Sf/Sb or Rf/Rb) greater
## or equal to 1 and that have 'Name' label as 'BLANK', 'DAP', ..., 'TRP'.
gastro.raw2 = selSpots(gastro.raw, sigNoise=1, rmFlag=NULL, gLabelsID="Name",
  remove=list(c("BLANK","DAP","LYS","PHE","Q_GENE","THR","TRP")))

## To see the number of spots that suvived the filtering above do
apply(gastro.raw2@UseSpots, 2, sum)

## To do the same filtering as above, also filtering flags marcke as 1
## and 4 do
gastro.raw2 = selSpots(gastro.raw, sigNoise=1, rmFlag=c(1,4), gLabelsID="Name",
  remove=list(c("BLANK","DAP","LYS","PHE","Q_GENE","THR","TRP")))

apply(gastro.raw2@UseSpots, 2, sum)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
