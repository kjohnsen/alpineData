\name{alpineData}
\docType{data}
\alias{alpineData}
\alias{alpineData-package}
\alias{ERR188297}
\alias{ERR188088}
\alias{ERR188204}
\alias{ERR188317}
\title{GAlignmentPairs for four samples from the GEUVADIS project}
\description{
  FASTQ files for four RNA-seq samples from the GEUVADIS project
  were downloaded, aligned and constructed as \code{GAlignmentPairs},
  according to the steps presented in \code{inst/scripts/make-data.R}
}
\usage{
ERR188297(metadata=FALSE)
ERR188088(metadata=FALSE)
ERR188204(metadata=FALSE)
ERR188317(metadata=FALSE)
}
\arguments{
  \item{metadata}{
    \code{logical} value indicating whether metadata only should be returned
    or if the resource should be loaded. Default behavior(metadata=FALSE) 
    loads the data.
  }
}
\examples{
  ERR188297()
}
\value{ These accessor functions return GAlignmentPairs objects }
\format{GAlignmentPairs}
\source{FASTQ files for four RNA-seq samples from the GEUVADIS project.}
\references{
  Lappalainen, T., et al. Transcriptome and genome sequencing uncovers
  functional variation in humans. Nature. 501(7468): 506-511 (2013)
  doi: 10.1038/nature12531
}
\keyword{datasets}
