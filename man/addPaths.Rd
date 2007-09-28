\name{addPaths}
\alias{addPaths}

\title{
  Function to load gene pathways into maigesPreRaw object
}

\description{
  This function read a directory and read files containing the gene
  pathways in TGF format. This format must have the genes of the pathway
  sequentially in lines numbered from 1, followed by a '\#' character
  that separate the nodes (given by genes) from edges, that must be
  specified as number of the origin gene followed by a space, the number
  of the final gene, another space and the weight of the iteration. This
  function stores the gene networks read in the slot \code{Paths} into
  objects of class \code{\link{maigesPreRaw}}.
}

\usage{
addPaths(data, folder="./", ext=".tgf")
}

\arguments{
  \item{data}{object of \code{\link{maigesPreRaw}} class.}
  \item{folder}{char string specifying the directory of gene groups. The
    function tests the presence or not of the final bar.}
  \item{ext}{string giving the extension of the files, defaults to
    '.tgf'. The function also tests the presence of the initial dot.}
}

\value{
  This function return another object of class
  \code{\link{maigesPreRaw}}, with the slot \code{Paths} actualised.
}

\details{
  If the \code{data} object already has gene networks with names equal
  to some someones that are been read, the nets with repeated names
  are not added. Warning messages are printed for every repeated group
  name.

  The \code{folder} directory must contain only one file for each pathway
  of interest. These files must be done in TGF format, as described into
  description above. The gene identification are matched with some
  column from \code{genemap} (see \code{\link{loadData}}).
}

\seealso{
  \code{\link{maigesPreRaw}}, \code{\link{addGeneGrps}}
}

\examples{
## Don't run because you don't have the pathways in a readable folder.
\dontrun{
gastro = addPaths(gastro, folder="geneNets", ext="tgf")
}}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}

\keyword{methods}
