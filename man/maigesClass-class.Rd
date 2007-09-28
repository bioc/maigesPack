\name{maigesClass-class}
\docType{class}
\alias{maigesClass-class}
\alias{maigesClass}

\title{
  maigesClass class, store results of discrimination (or classification) analysis
}

\description{
  This class defines a structure to store the results from
  discrimination analysis. This type of analysis can be done using the
  functions \code{\link{classifyLDA}}, \code{\link{classifySVM}},
  \code{\link{classifyKNN}}, \code{\link{classifyLDAsc}},
  \code{\link{classifySVMsc}} or \code{\link{classifyKNNsc}}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{numeric matrix giving the W values of the genes
      tested. This information is useful for doing plots of the cliques.}
    \item{\code{CV}:}{numeric vector that store the number of correct
      classifications in the leave-one-out cross validation procedure.}
    \item{\code{SVD}:}{numeric vector that store the singular value
      decomposition from Fisher linear discriminant analysis.}
    \item{\code{cliques}:}{character matrix that gives the genes that
      constitute the cliques returned. The rows of the matrix represent
      the cliques while the columns represent the genes that form the clique.}
    \item{\code{cliques.idx}:}{numeric matrix similar to the above,
      storing the indexes (onto \code{W} slot) of the genes.}
    \item{\code{method}:}{character string giving the method of
      discrimination analysis used.}
    \item{\code{Date}:}{character string giving the date and time that
      the object was generated.}
    \item{\code{V.info}:}{list containg three characters. The first one is
      a string containing the R version used when the object was
      created. The second is a char vector with base packages and the
      last one is another char vector with additional packages and
      version numbers.}
  }
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = 'maigesClass')}: plot method for
      \code{\link{maigesClass}} class. Display dispersion plots.}
    \item{print}{\code{signature(x = 'maigesClass')}: print method for
      \code{\link{maigesClass}} class.}
    \item{show}{\code{signature(x = 'maigesClass')}: show method for
      \code{\link{maigesClass}} class.}
    \item{summary}{\code{signature(x = 'maigesClass')}: summary method for
      \code{\link{maigesClass}} class.}
  }
}

\details{
  Objects of this class are produced by calling the functions
  \code{\link{classifyLDA}}, \code{\link{classifySVM}},
  \code{\link{classifyKNN}}, \code{\link{classifyLDAsc}},
  \code{\link{classifySVMsc}} or \code{\link{classifyKNNsc}} over an
  object of class \code{\link{maiges}} to search for cliques satisfying the
  criteria specified for classification.
}

\seealso{
  \code{\link{classifyLDA}}, \code{\link{classifySVM}},
  \code{\link{classifyKNN}}, \code{\link{classifyLDAsc}},
  \code{\link{classifySVMsc}} and \code{\link{classifyKNNsc}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
