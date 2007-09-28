\name{maiges-class}
\docType{class}
\alias{maiges-class}
\alias{maiges}

\title{
  maiges class, store normalised microarray datasets
}

\description{
  This class describes objects to handle ratio of intensities (\emph{W})
  and average of intensities (\emph{A}) values values and information
  about genes and samples used in the data. Objects of this class are
  created from class \code{\link{maigesRaw}} using the functions
  \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}} and/or
  \code{\link{normScaleMarray}} to do the normalisation.

  Here, the \emph{M=log(R/G)} value of intensity ratio was redefined as
  \emph{W=log(Test/Ref)}, where \emph{Test} and \emph{Ref} are the test
  and reference samples.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{numeric matrix containing the ratio values (in log2
      scale) between the test and reference sample intensities (W
      values). Spots are indexed by rows and samples by columns.}
    \item{\code{A}:}{numeric matrix containing the mean intensity
      values between test and reference samples (also in log2 scale). 
      Spots corresponding to rows and samples (or chips) corresponding
      to columns, too.}
    \item{\code{SD}:}{numeric matrix containing the standard deviation
      of W values when the lowess step is repeated several times.}
    \item{\code{IC1}:}{numeric matrix containing the left margin of
      confidence interval defined by repeated lowess during the
      normalisation step.}
    \item{\code{IC2}:}{numeric matrix containing the right margin of
      confidence interval defined above.}
    \item{\code{BadSpots}:}{logical vector specifying spots that was
      judged as bad ones. By default this slot is created as a vector of
      FALSEs with same length as number of spots.}
    \item{\code{UseSpots}:}{logical matrix indexing the spots to be used
      for normalisation.}
    \item{\code{GeneGrps}:}{a logical matrix with rows representing the
      spots and columns representing different gene groups. Each column
      give the index of spots in that gene group.}
    \item{\code{Paths}:}{list containing \code{\link[graph:graphNEL-class]{graphNEL}} objects
      specifying gene regulatory networks (or pathways). The first
      object in this list is a char string giving the gene label used to
      match the genes.}
    \item{\code{Layout}:}{a list containing the number of rows (\code{gridR}) and
      columns (\code{gridC}) of grids, the number of rows (\code{spotR})
      and columns (\code{spotC}) of spots inside each grid and the total
      number of spots.}
    \item{\code{Glabels}:}{data frame giving the gene labels. These
      labels are generally used during the data analysis.}
    \item{\code{Slabels}:}{data frame giving the sample labels. These
      labels are generally used during the data analysis.}
    \item{\code{Notes}:}{char string that receives any comment about the
      dataset. The dataset description is stored in this slot.}
    \item{\code{Date}:}{char string giving the date and hour that the
      object was created.}
    \item{\code{V.info}:}{list containg three characters. The first one is
      a string containing the R version used when the object was
      created. The second is a char vector with base packages and the
      last one is another char vector with additional packages and
      version numbers.}
  }
}

\details{
  This defines the main class of objects defined in
  this package. It is created from \code{\link{maigesRaw}} class
  using the normalisation functions \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}} and
  \code{\link{normScaleMarray}}. From this class
  of objects it is possible to do any type of analysis defined by
  several functions in this package. Also, it is possible to summarise
  spots (or samples) information using the function
  \code{\link{summarizeReplicates}}.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = 'maiges')}: subsetting operator for
      spots on the array or arrays in the batch, ensures that all slots
      are subset properly.}
    \item{boxplot}{\code{signature(x = 'maiges')}: boxplot method for
      \code{\link{maiges}} class. Display boxplots of the slides and
      print tip groups using package \emph{marray} or boxplots of
      one gene previously defined.}
    \item{dim}{\code{signature(x = 'maiges', value = 'numeric')}: get
      the dimensions of the object, numeric vector of length two.}
    \item{image}{\code{signature(x = 'maiges')}: image method for
      \code{\link{maiges}} class. Display colour representation of
      the slides using package \emph{marray}.}
    \item{plot}{\code{signature(x = 'maiges')}: plot method for
      \code{\link{maiges}} class. Display 'MA' plots.}
    \item{print}{\code{signature(x = 'maiges')}: print method for
      \code{\link{maiges}} class.}
    \item{show}{\code{signature(x = 'maiges')}: show method for
      \code{\link{maiges}} class.}
    \item{summary}{\code{signature(x = 'maiges')}: summary method for
      \code{\link{maiges}} class.}
  }
}

\seealso{
  \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}},
  \code{\link{normScaleMarray}} and \code{\link{summarizeReplicates}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
