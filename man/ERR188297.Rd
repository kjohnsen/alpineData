\name{ERR188297}
\docType{data}
\alias{ERR188088}
\alias{ERR188204}
\alias{ERR188317}
\title{GAlignmentPairs for four samples from the GEUVADIS project}
\description{ }
\usage{
gap <- readRDS("ERR188297.rds")
gap <- readRDS("ERR188088.rds")
gap <- readRDS("ERR188204.rds")
gap <- readRDS("ERR188317.rds")
}
\format{GAlignmentPairs}
\source{FASTQ files for four RNA-seq samples from the GEUVADIS project
  were downloaded, aligned and constructed as \code{GAlignmentPairs},
  according to the steps presented in \code{inst/scripts/make-data.R}
}
\references{
  Lappalainen, T., et al. Transcriptome and genome sequencing uncovers
  functional variation in humans. Nature. 501(7468): 506–511 (2013)
  doi: 10.1038/nature12531
}
\keyword{datasets}
