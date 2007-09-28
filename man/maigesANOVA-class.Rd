\name{maigesANOVA-class}
\docType{class}
\alias{maigesANOVA-class}
\alias{maigesANOVA}

\title{
  maigesANOVA class, extend maiges class to fit ANOVA models
}

\description{
  This class extends the class \code{\link{maiges}} adding two numerical
  matrices to fit ANOVA models and estimate parameters of
  interest. Additional to the existing slots of class
  \code{\link{maiges}} this class defines other two, described in \emph{Slots} section.
}

\section{Slots}{
  \describe{
    \item{\code{Dmatrix}:}{numeric matrix describing the design matrix
      to fit an ANOVA model.}
    \item{\code{Cmatrix}:}{numeric matrix describing the contrasts
      matrix to estimate parameters and contrasts of interest, after the
      model fitting.}
  }
}

\details{
  This class of objects is exactly the same as \code{\link{maiges}} with
  the two slots described in \emph{Slots} section. Objects of this class are generated from
  \code{\link{maiges}} class using the function \code{\link{designANOVA}}. So,
  it is possible to fit the model and estimate parameters or contrasts
  using the function \code{\link{deGenesANOVA}}.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = 'maigesANOVA')}: sub-setting operator for
      spots on the array or arrays in the batch, ensures that all slots
      are subset properly.}
    \item{boxplot}{\code{signature(x = 'maigesANOVA')}: boxplot method for
      \code{\link{maigesANOVA}} class. Display boxplots of the slides and
      print tip groups using package \emph{marray} or boxplots of
      one gene previously defined.}
    \item{dim}{\code{signature(x = 'maigesANOVA', value = 'numeric')}: get
      the dimensions of the object, numeric vector of length two.}
    \item{image}{\code{signature(x = 'maigesANOVA')}: image method for
      \code{\link{maigesANOVA}} class. Display colour representation of
      the slides using package \emph{marray}.}
    \item{plot}{\code{signature(x = 'maigesANOVA')}: plot method for
      \code{\link{maigesANOVA}} class. Display 'MA' plots.}
    \item{print}{\code{signature(x = 'maigesANOVA')}: print method for
      \code{\link{maigesANOVA}} class.}
    \item{show}{\code{signature(x = 'maigesANOVA')}: show method for
      \code{\link{maigesANOVA}} class.}
    \item{summary}{\code{signature(x = 'maigesANOVA')}: summary method for
      \code{\link{maigesANOVA}} class.}
  }
}

\seealso{
  \code{\link{designANOVA}}, \code{\link{deGenesANOVA}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
