metadata <- read.csv("../extdata/metadata.csv", stringsAsFactors=FALSE)

library(ensembldb)
gtf.file <- "Homo_sapiens.GRCh38.84.gtf"
txdb <- EnsDb(basename(gtf.file))

bam.files <- paste0("out/",metadata$Title,"_galignpairs.bam")
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bam.files)))

txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]

# pre-selected genes
selected.genes <- scan("selected.genes.txt",what="char")

one.iso.txs <- txdf$tx_id[txdf$gene_id %in% intersect(one.iso.genes, selected.genes)]
length(one.iso.txs)

ebt0 <- exonsBy(txdb, by="tx")
ebt <- ebt0[one.iso.txs]

# more than 1 exon
table(elementNROWS(ebt))
ebt <- ebt[elementNROWS(ebt) > 1]

# filter small genes and long genes
min.bp <- 800
max.bp <- 5000
gene.lengths <- sum(width(ebt))
summary(gene.lengths)
ebt <- ebt[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt)

set.seed(1)
ebt <- ebt[sample(length(ebt),8)]

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

library(devtools)
load_all("../../../alpine")

library(BSgenome.Hsapiens.NCBI.GRCh38)
library(Rsamtools)

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

fitpar <- mclapply(bam.files, function(bf) {
                     fitBiasModels(genes=ebt,
                                   bamfile=bf,
                                   fragtypes=fragtypes,
                                   genome=Hsapiens,
                                   models=models,
                                   readlength=readlength,
                                   minsize=minsize,
                                   maxsize=maxsize)
                   })

save(fitpar, file="fitpar.rda")
