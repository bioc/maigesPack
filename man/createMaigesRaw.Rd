\name{createMaigesRaw}
\alias{createMaigesRaw}

\title{
  Function to create objects of class maigesRaw
}

\description{
  Function to create objects of class \code{\link{maigesRaw}} from
  objects of class \code{\link{maigesPreRaw}}.
}

\usage{
createMaigesRaw(PreRaw, greenDataField, greenBackDataField, redDataField,
                redBackDataField, flagDataField, gLabelGrp, gLabelPath)
}

\arguments{
  \item{PreRaw}{object of class \code{\link{maigesPreRaw}} to be used to
    generate another object of class \code{\link{maigesRaw}}.}
  \item{greenDataField}{character string specifying the name of the
    \code{Data} slot from \code{PreRaw} that will be used to read the
    spot intensity values for green channel.}
  \item{greenBackDataField}{character string specifying the name of the
    \code{Data} slot from \code{PreRaw} that will be used to read the
    background intensity values for green channel.}
  \item{redDataField}{character string specifying the name of the
    \code{Data} slot from \code{PreRaw} that will be used to read the
    spot intensity values for red channel.}
  \item{redBackDataField}{character string specifying the name of the
    \code{Data} slot from \code{PreRaw} that will be used to read the
    background intensity values for red channel.}
  \item{flagDataField}{character string specifying the name of the
    \code{Data} slot from \code{PreRaw} that will be used to read the
    flag values.}
  \item{gLabelGrp}{character string with the gene label to match gene groups.}
  \item{gLabelPath}{character string with the gene label to match gene networks.}
}

\value{
  The result of this function is an object of class \code{\link{maigesRaw}}.
}

\seealso{
  \code{\link{loadData}}, \code{\link{maigesPreRaw}},
}

\examples{
## Load a little dataset
data(gastro)

## See the object
gastro

## Transform gastro into a maigesRaw object
gastro.raw = createMaigesRaw(gastro, greenDataField="Ch1.Mean",
  greenBackDataField="Ch1.B.Mean", redDataField="Ch2.Mean",
  redBackDataField="Ch2.B.Mean", flagDataField="Flags",
  gLabelGrp="GeneName", gLabelPath="GeneName")
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
