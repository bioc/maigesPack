\name{summarizeReplicates}
\alias{summarizeReplicates}

\title{
  Summarise microarray objects
}

\description{
  Function to summarise measures of a \code{\link{maiges}} class object,
  both by samples and genes.
}

\usage{
summarizeReplicates(object=NULL, gLabelID="GeneName", sLabelID="Sample",
                    funcS="mean", funcG="mean", rmBad=TRUE, keepEmpty=TRUE)
}

\arguments{
  \item{object}{object of class \code{\link{maiges}}.}
  \item{gLabelID}{character string giving the gene label ID to be used
    to summarise the data by rows.}
  \item{sLabelID}{character string giving the sample label ID to be used
    to summarise the data by columns.}
  \item{funcS}{character string specifying the function to be applied
    for sample replicates. Defaults to 'mean'. If NULL, no resume is
    done for samples.}
  \item{funcG}{character string specifying the function to be applied
    for genes (spots) replicates. Defaults to 'mean'. If NULL, no resume is
    done for genes..}
  \item{rmBad}{logical indicating if you want to remove or not bad
    spots, given by the slot \code{BadSpots} in argument \code{object}.}
  \item{keepEmpty}{logical indicating if you want to maintain spots
    with empty information.}
}
 
\details{
  This function takes the object of class \code{\link{maiges}} and
  resume the data by spots (rows) and samples (columns).
}

\value{
  This function returns another object of class \code{\link{maiges}}
  with replicates summarised to only one observation.
}

\examples{
## Loading the dataset
data(gastro)

## Summarising the data (maigesNorm class), replicated samples will be
## summarised by mean and genes by median
gastro.summ = summarizeReplicates(gastro.norm, gLabelID="GeneName",
  sLabelID="Sample", funcS="mean", funcG="median",
  keepEmpty=FALSE, rmBad=FALSE)

## To summarise genes by mean and keep the blank spots use
gastro.summ = summarizeReplicates(gastro.norm, gLabelID="GeneName",
  sLabelID="Sample", funcS="mean", funcG="mean",
  keepEmpty=TRUE, rmBad=FALSE)
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
