\name{addGeneGrps}
\alias{addGeneGrps}

\title{
  Function to load gene groups into maigesPreRaw object
}

\description{
  This function read a directory and read files containing the genes for
  specific gene groups. This files must have one gene per line. This
  function stores the gene groups read in the slot \code{GeneGrps} into
  objects of class \code{\link{maigesPreRaw}}.
}

\usage{
addGeneGrps(data, folder="./", ext=".txt")
}

\arguments{
  \item{data}{object of \code{\link{maigesPreRaw}} class.}
  \item{folder}{char string specifying the directory of gene groups. The
    function tests the presence or not of the final bar.}
  \item{ext}{string giving the extension of the files, defaults to
    '.txt'. The function also tests the presence of the initial dot.}
}

\value{
  This function returns another object of class
  \code{\link{maigesPreRaw}}, with the slot \code{GeneGroups} actualised.
}

\details{
  If the \code{data} object already has gene groups with names equal
  to some someones that are been read, the groups with repeated names
  are not added. Warning messages are printed for every repeated group
  name.

  The \code{folder} directory must contain only one file for each gene
  group of interest. These files must discriminate one gene per line. The
  identification of the genes must be done by one of the gene labels
  given by \code{genemap} (see \code{\link{loadData}}).
}

\seealso{
  \code{\link{maigesPreRaw}}, \code{\link{addPaths}}
}

\examples{
## Don't run because you don't have the gene sets in a readable folder.
\dontrun{
gastro = addGeneGrps(gastro, folder="geneGrps", ext="txt")
}}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
