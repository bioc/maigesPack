\name{maigesRelNetM-class}
\docType{class}
\alias{maigesRelNetM-class}
\alias{maigesRelNetM}

\title{
  maigesRelNetM class, store results of relevance network analysis
  (comparing two sample types)
}

\description{
  This class defines a structure to store the results of relevance network analysis,
  using a method that compare the correlation values between pairs of
  genes in two different sample types by a
  \emph{z}-transformation. Objects of this type are generated calling
  the function \code{\link{relNetworkM}}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{numeric matrix giving the W values of the genes
      tested. This information is useful for doing plots representing the
      correlation coefficient values.}
    \item{\code{Corr1}:}{numerical square matrix that store the
      correlation coefficient values between all pair of genes tested
      for the first sample type.}
    \item{\code{Corr2}:}{numerical matrix similar to the previous one,
      storing the correlation values for the second sample type.}
    \item{\code{DifP}:}{numerical matrix that store the p-value for the
      comparison between the correlation values of both sample types tested.}
    \item{\code{types}:}{character vector giving the names of sample
      types tested.}
    \item{\code{Slabel}:}{string with the sample label used to define
      the two sample types tested.}
    \item{\code{Date}:}{character string giving the date and time that
      the object was generated.}
    \item{\code{V.info}:}{list containg three characters. The first one is
      a string containing the R version used when the object was
      created. The second is a char vector with base packages and the
      last one is another char vector with additional packages and
      version numbers.}
  }
}

\details{
  Objects of this class are produced by calling the functions
  \code{\link{relNetworkM}} over an object of class \code{\link{maiges}}
  to search for pairs of genes with altered correlation values between
  two sample types tested.
}

\section{Methods}{
  \describe{
    \item{image}{\code{signature(x = 'maigesRelNetM')}: image method for
      \code{\link{maigesRelNetM}} class. Display colour representation of
      the correlation and p-values matrices.}
    \item{plot}{\code{signature(x = 'maigesRelNetM')}: plot method for
      \code{\link{maigesRelNetM}} class. Display circular graphs
      representing the relevance networks.}
    \item{print}{\code{signature(x = 'maigesRelNetM')}: print method for
      \code{\link{maigesRelNetM}} class.}
    \item{show}{\code{signature(x = 'maigesRelNetM')}: show method for
      \code{\link{maigesRelNetM}} class.}
    \item{summary}{\code{signature(x = 'maigesRelNetM')}: summary method for
      \code{\link{maigesRelNetM}} class.}
  }
}

\seealso{
  \code{\link{relNetworkM}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
