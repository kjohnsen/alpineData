library(speedglm) # this first or AnnotationDbi's select() will be masked
library(GenomicAlignments)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(splines)
library(Matrix)
library(devtools)
library(graph)
library(RBGL)
load_all("../../alpine")

if (FALSE) {
  gtf <- "genesStandardTagValue.gtf"
  system.time({ txdb <- makeTxDbFromGFF(gtf) })
  saveDb(txdb, file="genesStandard.sqlite")
}

txdb <- loadDb("genesStandard.sqlite")

## 

samps <- read.delim("samples.txt", stringsAsFactors=FALSE)
samps$id <- samps$Comment.ENA_RUN.
bamfiles <- paste0("/n/irizarryfs01/mlove/geuvadis/star/",samps$id,
                   ".Aligned.sortedByCoord.out.bam")
names(bamfiles) <- samps$id
stopifnot(all(file.exists(bamfiles)))

gene.ids <- keys(txdb, keytype="GENEID")
txdf <- select(txdb, keys=gene.ids,
               columns=c("TXID","TXNAME","TXCHROM"), keytype="GENEID")
txdf <- txdf[txdf$GENEID != "",]

tab.tx <- table(txdf$GENEID)
single.tx.genes <- names(tab.tx[tab.tx == 1])
single.tx.txs <- txdf$TXNAME[txdf$GENEID %in% single.tx.genes]

ebt <- exonsBy(txdb, by="tx")
txname <- txdf$TXNAME[ match(names(ebt), txdf$TXID) ]
names(ebt) <- txname

ebt <- ebt[single.tx.txs]
# more than 1 exons
ebt <- ebt[elementLengths(ebt) > 1]

# filter small genes and long genes
min.bp <- 800
max.bp <- 5000
gene.lengths <- sum(width(ebt))
summary(gene.lengths)
ebt <- ebt[gene.lengths > min.bp & gene.lengths < max.bp]

# gene counts
min.count <- 200
max.count <- 15000
load("featureCounts.rda")
k <- rowMeans(fc$count)
txdf$count <- numeric(nrow(txdf))
txdf$count[txdf$GENEID %in% names(k)] <- k[match(txdf$GENEID, names(k))]

mcols(ebt)$count <- txdf$count[match(names(ebt), txdf$TXNAME)]
summary(mcols(ebt)$count)

ebt <- ebt[ mcols(ebt)$count > min.count & mcols(ebt)$count < max.count ]
length(ebt)

set.seed(1)
ebt <- ebt[sample(names(ebt), 100)]

# switch
extraStuff <- TRUE

if (!extraStuff) {
models <- list(
  "GC_str" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  GC40.80 + GC40.90 + GC20.80 + GC20.90 +
  gene",
  offset=c("fraglen"))
  )
} else {
models <- list(
  "readstart" = list(formula = "count ~ ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen","vlmm")),
  "all" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  GC40.80 + GC40.90 + GC20.80 + GC20.90 +
  gene",
  offset=c("fraglen","vlmm"))
  )
}

models

getReadlength(bamfiles)

options(mc.cores=4)

minsize <- 80
maxsize <- 350
readlength <- 75 

gene.names <- names(ebt)
names(gene.names) <- gene.names
fragtypes <- mclapply(gene.names, function(gene.name) {
               buildFragtypesFromExons(exons=ebt[[gene.name]], genome=Hsapiens,
               readlength=readlength, minsize=minsize, maxsize=maxsize, vlmm=extraStuff)
             })
fitpar <- mclapply(bamfiles, function(bamfile) {
            fitModelOverGenes(genes=ebt, bamfile=bamfile, fragtypes=fragtypes,
            genome=Hsapiens, models=models, readlength=readlength,
            minsize=minsize, maxsize=maxsize)
          })

sapply(fitpar, names)

if (!extraStuff) {
  save(fitpar, file="second_submit/fitpar_gc_str.rda")
} else {
  save(fitpar, file="second_submit/fitpar_all.rda")
}

