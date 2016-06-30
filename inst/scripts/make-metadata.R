# the GEUVADIS sample metadata is downloaded from ArrayExpress:
# browse: http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/
# download: http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt
samps <- read.delim("E-GEUV-1.sdrf.txt",header=TRUE)
samps <- samps[!duplicated(x$Derived.Array.Data.File),]
samps$date <- factor(sub(".*_(.*)_.*", "\\1", samps$Assay.Name))
table(samps$date, samps$Performer, samps$Characteristics.population.)
samps <- samps[samps$Characteristics.population. == "TSI",]
table(samps$date, samps$Performer)
perfs <- c("UNIGE","CNAG_CRG")
n <- 15
idx <- c(head(which(samps$Performer == perfs[1] & samps$date == "111124"), n),
         head(which(samps$Performer == perfs[2] & samps$date == "111215"), n))
samps <- samps[idx,]
samps$date <- droplevels(samps$date)
samps$Performer <- droplevels(samps$Performer)
table(samps$date, samps$Performer)
