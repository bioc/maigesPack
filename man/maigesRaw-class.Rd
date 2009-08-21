\name{maigesRaw-class}
\docType{class}
\alias{maigesRaw-class}
\alias{maigesRaw}

\title{
  maigesRaw class, store raw microarray datasets
}

\description{
  This class describes objects to handle intensity values and information
  about genes and samples used in the data. Objects of this class are
  obtained from class \code{\link{maigesPreRaw}} using function
  \code{\link{createMaigesRaw}}.
}

\section{Slots}{
  \describe{
    \item{\code{Sf}:}{numeric matrix containing the test samples
      spots intensity values, spots corresponding to rows and data
      tables (or chips) corresponding to columns.}
    \item{\code{Sb}:}{numeric matrix containing the test samples
      local background values, spots corresponding to rows and data
      tables (or chips) corresponding to columns.}
    \item{\code{Sdye}:}{character vector with length equal the length of
      data files (chips), specifying the channel ('ch1' or 'ch2') used
      to label each test sample in each chip.}
    \item{\code{Rb}:}{numeric matrix containing the reference samples
      spots intensity values, spots corresponding to rows and data
      tables (or chips) corresponding to columns.}
    \item{\code{Rf}:}{numeric matrix containing the reference samples
      spots intensity values, spots corresponding to rows and data
      tables (or chips) corresponding to columns.}
    \item{\code{Rdye}:}{character vector with length equal the length of
      data files (chips), specifying the channel ('ch1' or 'ch2') used
      to label each reference sample in each chip.}
    \item{\code{Flag}:}{matrix specifying the flags for the spots, as in
      the matrices above, rows and columns correspond to genes and
      samples, respectively. Type o value may be numeric or character,
      depending of the flags used.}
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
  This class of objects defines a real raw object that is used to do the
  normalisation. Objects of this class are generated from objects of
  class \code{\link{maigesPreRaw}} using the function
  \code{\link{createMaigesRaw}}. Here it is possible to do several
  plots for exploratory analysis using functions from
  \emph{marray} package. Using the function
  \code{\link{selSpots}}, you select spots to use in the normalisation
  method, that is done by the functions \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}} and
  \code{\link{normScaleMarray}}.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = 'maigesRaw')}: subsetting operator for
      spots on the array or arrays in the batch, ensures that all slots
      are subset properly.}
    \item{boxplot}{\code{signature(x = 'maigesRaw')}: boxplot method for
      \code{\link{maigesRaw}} class. Display boxplots of the slides and
      print tip groups using package \emph{marray}.}
    \item{dim}{\code{signature(x = 'maigesRaw', value = 'numeric')}: get
      the dimensions of the object, numeric vector of length two.}
    \item{image}{\code{signature(x = 'maigesRaw')}: image method for
      \code{\link{maigesRaw}} class. Display colour representation of
      the slides using package \emph{marray}.}
    \item{plot}{\code{signature(x = 'maigesRaw')}: plot method for
      \code{\link{maigesRaw}} class. Display \emph{MA} plots.}
    \item{print}{\code{signature(x = 'maigesRaw')}: print method for
      \code{\link{maigesRaw}} class.}
    \item{show}{\code{signature(x = 'maigesRaw')}: show method for
      \code{\link{maigesRaw}} class.}
    \item{summary}{\code{signature(x = 'maigesRaw')}: summary method for
      \code{\link{maigesRaw}} class.}
  }
}

\seealso{
  \code{\link{createMaigesRaw}}, \code{\link{selSpots}},
  \code{\link{normLoc}}, \code{\link{normOLIN}},
  \code{\link{normRepLoess}}, \code{\link{normScaleLimma}} and
  \code{\link{normScaleMarray}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
