\name{maigesPreRaw-class}
\docType{class}
\alias{maigesPreRaw-class}
\alias{maigesPreRaw}

\title{
  maigesPreRaw class, store pre raw microarray datasets
}

\description{
  This class describes objects to handle data values and information
  about genes and samples used in the data. Here, it is possible to put
  any field from data tables that you want.
}

\section{Slots}{
  \describe{
    \item{\code{Data}:}{contains a list with matrices. Each matrix has
      \code{nrow} = number of spots and \code{ncol} = number of data
      tables (or chips) with the numerical values of each data field
      specified by the user in the configuration file used as parameter
      for \code{\link{loadData}} function.}
    \item{\code{GeneGrps}:}{a list containing character vectors. These
      vectors specify genes from the gene groups that must be
      studied. Each vector is used as a gene group. The names of the
      groups are catch from the names of the vectors from the list,
      that is catch from file names.}
    \item{\code{Paths}:}{list containing \code{\link[graph:graphNEL-class]{graphNEL}} objects
      specifying gene regulatory networks (or pathways). As in the
      previous slot, the names are catch from file names and stored as
      names of the elements from the list.}
    \item{\code{Layout}:}{a list containing the number of rows (\code{gridR}) and
      columns (\code{gridC}) of grids, the number of rows (\code{spotR})
      and columns (\code{spotC}) of spots inside each grid and the total
      number of spots.}
    \item{\code{Glabels}:}{data frame giving the gene labels. These
      labels are generally used during the data analysis.}
    \item{\code{Slabels}:}{data frame giving the sample labels. These
      labels are generally used during the data analysis.}
    \item{\code{BadSpots}:}{logical vector specifying spots that was
      judged as bad ones. By default this slot is created as a vector of
      FALSEs with same length as number of spots.}
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
  This class of objects is intended to define an intermediate step
  between the external world (data tables, gene groups, gene networks,
  etc) and real raw objects given by the class
  \code{\link{maigesRaw}}. In this intermediate step the user can do any
  type of exploratory analysis and set bad spots (at slot
  \code{BadSpots}), or to do any type of calculation they judge
  important in the fields of \code{Data} slot.

  This class of objects is created by a call from
  \code{\link{loadData}}. The functions \code{\link{addGeneGrps}} and
  \code{\link{addPaths}} may also be used to add information for gene
  groups and gene networks. After the exploratory analysis it must be
  converted in an object of class \code{\link{maigesRaw}} using the
  function \code{\link{createMaigesRaw}}.

  If you have the package \emph{Rgraphviz} installed and working it is
  possible to display the graphs stored in the slot \code{Paths} using
  the command \code{\link[Rgraphviz:plot-methods]{plot}}.
}

\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = 'maigesPreRaw')}: print method for
      \code{\link{maigesPreRaw}} class.}
    \item{show}{\code{signature(x = 'maigesPreRaw')}: show method for
      \code{\link{maigesPreRaw}} class.}
    \item{summary}{\code{signature(x = 'maigesPreRaw')}: summary method for
      \code{\link{maigesPreRaw}} class.}
  }
}

\seealso{
  \code{\link{loadData}}, \code{\link{addGeneGrps}},
  \code{\link{addPaths}}, \code{\link{createMaigesRaw}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
