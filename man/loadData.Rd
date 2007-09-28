\name{loadData}
\alias{loadData}
\title{
  Load cDNA microarray data tables
}
\description{
  This function loads a cDNA microarray dataset into a temporary
  \code{\link{maigesPreRaw}} object.
}
\usage{
loadData(fileConf=paste(R.home(), "library/maiges/doc/gastro/load_gastro.conf", sep="/"))
}
\arguments{
  \item{fileConf}{string specifying a file name containing the
    parameters to load data. This file must contain all the information
    necessary to load the data, which are the following:\cr
    \describe{
      \item{\emph{dataDir}:}{specify a folder name containing the data files to
	be loaded. The function tests the presence or not of the final bar.}
      \item{\emph{ext}:}{string specifying the extension of the tables (if the
	\emph{sampleFile} below don't contain this information). You don't need to
	put the \emph{dot} onto string beginning, the function tests this automatically.}
      \item{\emph{sampleFile}:}{string containing the file name with
	the descriptions of the biological samples hybridised, including 
	the respective intensity data files. This file must be 
	spreadsheet-like separated by tabs in a plain text
	format. The column fields \bold{'File'} and \bold{'Ref'} are
	mandatory (with exactly these names). The first one describes the files containing the
	numerical data and the second one describes the channel
	used to label the reference sample, must be 'green' or 'red' and
	they are not case sensitive.}
      \item{\emph{datasetId}:}{string with a dataset identification.}
      \item{\emph{geneMap}:}{as in \emph{sampleFile}, this item is a character
	string giving a file name. This file must describe the genes on the
	slides. Also it must be a plain text spreadsheet-like separated by
	tabs. There are no mandatory field, but it is strongly recommended that you
	specify some fields containing gene names, genbank ID, cluster
	ID and gene annotations for a nice gene identification.}
      \item{\emph{headers}:}{character string (in the R format) specifying the
	column fields from data files you want to load.}
      \item{\emph{skip}:}{number of lines to be skipped in the numeric tables.}
      \item{\emph{sep}:}{character that separates the fields on the numeric tables.}
      \item{\emph{gridR}:}{number of print tip rows inside the slide.}
      \item{\emph{gridC}:}{number of print tip columns inside the slide.}
      \item{\emph{printTipR}:}{number of rows inside each print tip.}
      \item{\emph{printTipC}:}{number of columns inside each print
	tip.}
    }
    
    You can see an example of this configuration file in
    RHOME/library/maiges/doc/gastro/load\_gastro.conf.}
}

\details{
  This function takes the file name with initial arguments and load the
  dataset specified by this config file. It generate a
  \code{\link{maigesPreRaw}} object. During the process the function
  writes a file named load.out on your working folder, that is a log of
  the process, that you can check and verify if all was done
  correctly. Obviously, the parameters \emph{dataDir},
  \emph{sampleFile}, \emph{geneMap}, \emph{sep}, \emph{gridR},
  \emph{gridC}, \emph{printTipR}, \emph{printTipC} and \emph{headers}
  must be specified. All other parameters may be specified as NULL and,
  if so, they are ignored. It is possible to specify any fields that you
  want in the headers parameter, but it is  strongly recommended that
  you specify the fields of spot intensity and background for both
  channels and the filed giving quality weights for all spots.
}
\value{
  This function returns a \code{\link{maigesPreRaw}} object containing the dataset loaded.

  Once an object of class \code{\link{maigesPreRaw}} was generated, you
  may use the functions \code{\link{addGeneGrps}} and
  \code{\link{addPaths}} to load informaation about gene groups and gene
  networks, respectively.
}

\examples{
## Don't run because you don't have data tables.
\dontrun{
gastro = loadData(fileConf="load_gastro.conf")
}}

\author{
  Gustavo H. Esteves <\email{gesteves@vision.ime.usp.br}>
}
\keyword{file}
