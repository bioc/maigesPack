\name{maigesDE-class}
\docType{class}
\alias{maigesDE-class}
\alias{maigesDE}

\title{
  maigesDE class, store results of differential gene expression
  analysis
}

\description{
  This class defines a structure to store the results from differential
  expression analysis. This type of analysis can be done using the
  functions \code{\link{deGenes2by2Ttest}},
  \code{\link{deGenes2by2Wilcox}} and \code{\link{deGenes2by2BootT}}
  for analysis between two biological sample types, or function
  \code{\link{deGenesANOVA}}, for analysis between more than two sample
  types using ANOVA models.
}

\section{Slots}{
  \describe{
    \item{\code{GeneInfo}:}{data frame containing the information
      regarding the genes from the dataset.}
    \item{\code{SampleInfo}:}{data frame similar to the \code{GeneInfo}
      above, but containing information about the biological sample
      types used in the analysis.}
    \item{\code{fold}:}{numerical matrix containing the fold values
      (mean difference between the sample types) when two sample types
      were compared. Each matrix column gives one specific test.}
    \item{\code{stat}:}{numerical matrix giving the statistic of the
      tests that were done. Each column represents on test.}
    \item{\code{p.value}:}{numerical matrix giving the p-values of the
      statistical tests. Again, each column represents different tests.}
    \item{\code{factors}:}{character string giving the biological factors
      used in ANOVA model (when ANOVA models were used).}
    \item{\code{test}:}{character string describing the test done.}
    \item{\code{Date}:}{character string giving date and time that the
      object were created.}
    \item{\code{V.info}:}{list containg three characters. The first one is
      a string containing the R version used when the object was
      created. The second is a char vector with base packages and the
      last one is another char vector with additional packages and
      version numbers.}
  }
}

\details{
  Objects of this class are produced by calling the functions \code{\link{deGenes2by2Ttest}},
  \code{\link{deGenes2by2Wilcox}} and \code{\link{deGenes2by2BootT}}
  over an object of class \code{\link{maiges}} or by calling the
  function \code{\link{deGenesANOVA}} over an object of class
  \code{\link{maigesANOVA}} fitting an ANOVA model to the microarray dataset.
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = 'maigesDE')}: plot method for
      \code{\link{maigesDE}} class. Display volcano plots.}
    \item{print}{\code{signature(x = 'maigesDE')}: print method for
      \code{\link{maigesDE}} class.}
    \item{show}{\code{signature(x = 'maigesDE')}: show method for
      \code{\link{maigesDE}} class.}
    \item{summary}{\code{signature(x = 'maigesDE')}: summary method for
      \code{\link{maigesDE}} class.}
  }
}

\seealso{
  \code{\link{deGenes2by2Ttest}}, \code{\link{deGenes2by2Wilcox}},
  \code{\link{deGenes2by2BootT}}, \code{\link{deGenesANOVA}}.
}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{classes}
