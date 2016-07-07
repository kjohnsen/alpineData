# download the following FASTQ files to local dir
metadata <- read.csv("../extdata/metadata.csv", stringsAsFactors=FALSE)
files <- metadata$SourceUrl
all.files <- c(files, sub("_1.fastq.gz","_2.fastq.gz",files))

# run HISAT2
# https://ccb.jhu.edu/software/hisat2/index.shtml
# with the genome: H. sapiens, Ensembl GRCh38 genome_tran
# ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz

