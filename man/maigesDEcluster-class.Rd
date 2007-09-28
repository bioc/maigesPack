\name{maigesDEcluster-class}
\docType{class}
\alias{maigesDEcluster-class}
\alias{maigesDEcluster}

\title{
  maigesDEcluster class, store results of differential gene expression
  analysis
}

\description{
  This class extends the class \code{\link{maigesDE}} to store results
  from differential expression analysis. It is produced in the same way
  as the latter, but has one more slot containing the W values. This is
  useful to do cluster analysis. Together the slots of the
  \code{\link{maigesDE}} class, it has one more given in \emph{Slots} section.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{numeric matrix containing the W values for the
      genes and samples used in the analysis.}
  }
}

\details{
  Objects of this class are produced  in the same way as
  \code{\link{maigesDE}}. The addition of the slot \code{W} turn
  possible to do cluster analysis in this class of objects using the
  functions \code{\link{hierMde}}, \code{\link{somMde}} and
  \code{\link{kmeansMde}}, selecting genes according to the results of gene
  expression analysis.
}

\section{Methods}{
  \describe{
    \item{boxplot}{\code{signature(x = 'maigesDEcluster')}: boxplot method for
      \code{\link{maigesDEcluster}} class. Display boxplots of one gene
      previously defined.}
    \item{plot}{\code{signature(x = 'maigesDEcluster')}: plot method for
      \code{\link{maigesDEcluster}} class. Display \emph{MA} plots.}
    \item{print}{\code{signature(x = 'maigesDEcluster')}: print method for
      \code{\link{maigesDEcluster}} class.}
    \item{show}{\code{signature(x = 'maigesDEcluster')}: show method for
      \code{\link{maigesDEcluster}} class.}
    \item{summary}{\code{signature(x = 'maigesDEcluster')}: summary method for
      \code{\link{maigesDEcluster}} class.}
  }
}

\seealso{
  \code{\link{deGenes2by2Ttest}}, \code{\link{deGenes2by2Wilcox}},
  \code{\link{deGenes2by2BootT}}, \code{\link{deGenesANOVA}},
  \code{\link{hierMde}}, \code{\link{somMde}}, \code{\link{kmeansMde}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
