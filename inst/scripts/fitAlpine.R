library(ensembldb)
txdb <- EnsDb(basename(gtf.file))

metadata <- read.csv("../extdata/metadata.csv", stringsAsFactors=FALSE)
bam.files <- paste0("out/",metadata$Title,".bam")
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bamfiles)))

txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]
one.iso.txs <- txdf$tx_id[txdf$gene_id %in% one.iso.genes]

# pre-selected genes
selected.genes <- scan("selected.genes.txt",what="char")

ebt <- exonsBy(txdb, by="tx")
ebt <- ebt[intersect(one.iso.txs, selected.genes)]

# more than 1 exon
ebt <- ebt[elementLengths(ebt) > 1]

# filter small genes and long genes
min.bp <- 800
max.bp <- 5000
gene.lengths <- sum(width(ebt))
summary(gene.lengths)
ebt <- ebt[gene.lengths > min.bp & gene.lengths < max.bp]

models <- list(
  "GC" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen")),
  "all" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen","vlmm"))
  )

# need to load alpine
# need genome

# ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

library(Rsamtools)
FaFile()
seqinfo()
getSeq()

# library(BSgenome.Hsapiens.NCBI.GRCh38)

getReadlength(bamfiles)

options(mc.cores=4)

minsize <- 100
maxsize <- 300
readlength <- 75 

gene.names <- names(ebt)
names(gene.names) <- gene.names
fragtypes <- mclapply(gene.names, function(gene.name) {
                        buildFragtypes(exons=ebt[[gene.name]],
                                       genome=Hsapiens,
                                       readlength=readlength,
                                       minsize=minsize,
                                       maxsize=maxsize)
                      })

fitpar <- mclapply(bamfiles, function(bf) {
                     fitBiasModel(genes=ebt,
                                  bamfile=bf,
                                  fragtypes=fragtypes,
                                  genome=Hsapiens,
                                  models=models,
                                  readlength=readlength,
                                  minsize=minsize,
                                  maxsize=maxsize)
                   })

save(fitpar, file="fitpar.rda")

