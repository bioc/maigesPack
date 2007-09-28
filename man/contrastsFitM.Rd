\name{contrastsFitM}
\alias{contrastsFitM}

\title{
  Compute Contrasts from Linear Model Fit
}

\description{
  Given a linear model fit to microarray data, compute estimated
  coefficients and standard errors for a given set of contrasts. This
  function was adapted from \code{\link[limma]{contrasts.fit}} of the
  \emph{limma} package.
}

\usage{
contrastsFitM(fit, contrasts)
}

\arguments{
  \item{fit}{an \code{\link[limma:marraylm]{MArrayLM}} object or a list
    object produced by the function \code{\link[limma]{lm.series}} or
    equivalent.}
  \item{contrasts}{numeric matrix with row corresponding to coefficients
    in fit and columns containing contrasts. May be a vector if there is
    only one contrast.}
}

\value{
  The result of this function is an object of class
  \code{\link[limma:marraylm]{MArrayLM}}.
}

\details{
  This function was adapted from the equivalent
  \code{\link[limma]{contrasts.fit}} limma's function to do the linear
  model fit without use the empirical Bayes method given by the function
  \code{\link[limma:ebayes]{eBayes}}.
}

\seealso{
  \code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}},
  \code{\link[limma:ebayes]{eBayes}}, \code{\link[limma:marraylm]{MArrayLM}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>, adapted from
  the limma's function.
}

\keyword{methods}
