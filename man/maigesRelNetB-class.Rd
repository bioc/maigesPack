\name{maigesRelNetB-class}
\docType{class}
\alias{maigesRelNetB-class}
\alias{maigesRelNetB}

\title{
  maigesRelNetB class, store results of relevance network analysis (Butte's method)
}

\description{
  This class defines a structure to store the results of relevance network analysis,
  using a method that was proposed by Butte et al. (2000), see the
  reference below. Objects of this type are generated calling
  the function \code{\link{relNetworkB}}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{numeric matrix giving the W values of the genes
      tested. This information is useful for doing plots representing the
      correlation coefficient values.}
    \item{\code{Corr}:}{numerical square matrix that store the
      correlation coefficient values between all pair of genes tested.}
    \item{\code{Pval}:}{numerical square matrix that store the p-value
      for each correlation value between the pairs of genes.}
    \item{\code{maxB}:}{numerical square matrix that store the maximum
      (in absolute value) of the correlation values calculated by
      permutation values, as described by Butte et al. (2000).}
    \item{\code{types}:}{character vector giving the name of sample
      type tested.}
    \item{\code{Slabel}:}{string with the sample label used to define
      the sample types tested.}
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
  Objects of this class are produced by calling the function
  \code{\link{relNetworkB}} over an object of class \code{\link{maiges}}
  to search for pairs of genes with significant correlation values in a
  specific sample type, according to the method proposed by Butte et
  al. (2000).
}

\section{Methods}{
  \describe{
    \item{image}{\code{signature(x = 'maigesRelNetB')}: image method for
      \code{\link{maigesRelNetB}} class. Display colour representation of
      the correlation matrix.}
    \item{plot}{\code{signature(x = 'maigesRelNetB')}: plot method for
      \code{\link{maigesRelNetB}} class. Display a circular graph
      representing the relevance network.}
    \item{print}{\code{signature(x = 'maigesRelNetB')}: print method for
      \code{\link{maigesRelNetB}} class.}
    \item{show}{\code{signature(x = 'maigesRelNetB')}: show method for
      \code{\link{maigesRelNetB}} class.}
    \item{summary}{\code{signature(x = 'maigesRelNetB')}: summary method for
      \code{\link{maigesRelNetB}} class.}
  }
}

\references{
  Butte, A.J.; Tamayo, P.; Slonim, D.; Golub, T.R. and Kohane,
  I.S. Discovering functional relationships between RNA expression and
  chemotherapeutic susceptibility using relevance networks. \bold{PNAS}, 97,
  12182-12186, 2000. (\url{http://www.pnas.org/cgi/content/full/97/22/12182})
}

\seealso{
  \code{\link{relNetworkB}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
