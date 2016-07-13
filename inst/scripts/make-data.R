# download the following FASTQ files to local dir
metadata <- read.csv("../extdata/metadata.csv", stringsAsFactors=FALSE)
files <- metadata$SourceUrl
all.files <- c(files, sub("_1.fastq.gz","_2.fastq.gz",files))

# run HISAT2 to align paired-end reads to genome
# https://ccb.jhu.edu/software/hisat2/index.shtml
# with the genome: H. sapiens, Ensembl GRCh38 genome_tran (this is version 84)
# ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz

# hisat2 -x $genome -1 $dir/$f\_1.fastq.gz -2 $dir/$f\_2.fastq.gz -p 10 > out/$f.sam
# samtools view -@ 4 -b out/$f.sam | samtools sort -@ 4 -O bam -T tmp - > out/$f.bam
# samtools index out/$f.bam

# resulting in BAM files
bam.files <- paste0("out/",metadata$Title,".bam")

# download the corresponding Ensembl GTF file 
# ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gtf.file <- "Homo_sapiens.GRCh38.84.gtf"

# count paired-end reads to genes using featureCounts
library(Rsubread)
if (!file.exists("featurecounts.rda")) {
  fc <- featureCounts(files=bam.files,
                      annot.ext=gtf.file,
                      isGTFAnnotationFile=TRUE,
                      isPairedEnd=TRUE,
                      autosort=FALSE)
  save(fc, file="featurecounts.rda")
} else {
  load("featurecounts.rda")
}

# select a set of genes with moderately high counts, but not too high
idx <- rowMeans(fc$counts) > 200 & rowMeans(fc$counts) < 10000
sum(idx)
mid.count.genes <- rownames(fc$counts)[idx]

# load an Ensembl TxDb
library(ensembldb)
if (!file.exists(basename(gtf.file))) {
  ensDbFromGtf(gtf.file, outfile=basename(gtf.file))
}
txdb <- EnsDb(basename(gtf.file))
txdf <- transcripts(txdb, return.type="DataFrame")

# names of single isoform genes
tab <- table(txdf$gene_id)
single.iso.genes <- names(tab)[tab == 1]

# subset to 200 genes with moderate counts and
# possessing a single isoform
g <- genes(txdb)
g <- keepSeqlevels(g, sort(c(as.character(1:22),"X","Y","MT")))
intersect.genes <- intersect(intersect(mid.count.genes, single.iso.genes), names(g))
length(intersect.genes)

set.seed(1)
g.sub <- sort(g[sample(intersect.genes,100)])
write(names(g.sub), file="selected.genes.txt")

# extract paired-end reads covering these genes for each BAM
library(GenomicAlignments)
for (i in seq_len(nrow(metadata))) {
  pt <- proc.time()
  gap <- readGAlignmentPairs(bam.files[i],
                             param=ScanBamParam(which=g.sub))
  et <- unname((proc.time() - pt)[3])
  print(paste(i,round(et),length(gap)))
  saveRDS(gap, file=paste0("out/",metadata$Title[i],".rds"))
}
