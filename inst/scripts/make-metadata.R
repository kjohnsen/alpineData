# the GEUVADIS sample metadata is downloaded from ArrayExpress:
# browse: http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/
# download: http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
samps <- read.delim("E-GEUV-1.sdrf.txt",header=TRUE)
samps <- samps[!duplicated(samps$Derived.Array.Data.File),]
samps$date <- factor(sub(".*_(.*)_.*", "\\1", samps$Assay.Name))
table(samps$date, samps$Performer, samps$Characteristics.population.)
samps <- samps[samps$Characteristics.population. == "TSI",]
table(samps$date, samps$Performer)
perfs <- c("UNIGE","CNAG_CRG")
n <- 2
idx <- c(head(which(samps$Performer == perfs[1] & samps$date == "111124"), n),
         head(which(samps$Performer == perfs[2] & samps$date == "111215"), n))
samps <- samps[idx,]
samps$date <- droplevels(samps$date)
samps$Performer <- droplevels(samps$Performer)
table(samps$date, samps$Performer)
run <- samps$Comment.ENA_RUN.
N <- 2*n

metadata <- data.frame(
  Title=run,
  Description=paste("Subset of aligned reads from sample",run),
  BiocVersion=rep("3.4",N),
  Genome=rep("GRCh38",N),
  SourceType=rep("FASTQ",N),
  SourceUrl=samps$Comment.FASTQ_URI.,
  SourceVersion=rep("1",N),
  Species=rep("Homo sapiens",N),
  TaxonomyId=rep("9606",N),
  Coordinate_1_based=rep(TRUE,N),
  DataProvider=rep("GEUVADIS",N),
  Maintainer=rep("Michael Love <michaelisaiahlove@gmail.com>",N),
  RDataClass=rep("GAlignmentPairs",N),
  DispatchClass=rep("Rda",N),
  ResourceName=paste0(run, ".rda"),
  Performer=samps$Performer,
  Date=samps$date,
  Population=samps$Characteristics.population.)

write.csv(metadata, file="../extdata/metadata.csv", row.names=FALSE)
