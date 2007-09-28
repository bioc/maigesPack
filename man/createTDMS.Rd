\name{createTDMS}
\alias{createTDMS}

\title{
  Create a tab delimited file for TIGR MeV
}

\description{
  This function gets an object of class \code{\link{maiges}} an save a
  tab delimited file containing the W values to be load by TIGR MeV
  software, in the format TDMS file.
}

\usage{
createTDMS(data=NULL, sLabelID=names(data@Slabels)[1], file="data.txt")
}

\arguments{
  \item{data}{object of class \code{\link{maiges}} to be saved as a TDMS
    file.}
  \item{sLabelID}{character string giving the sample label ID to be used
    to label the samples in the TDMS file. Defaults to the first
    sample label ID.}
  \item{file}{character string specifying the file name to where TDMS
    file must be saved. Defaults to 'data.txt'.}
}

\value{
  This function save an ASCII file and do not return any object or value.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
