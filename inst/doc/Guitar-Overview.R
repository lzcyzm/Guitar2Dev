### R code from vignette source 'Guitar-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: Genomic Features
###################################################
library(Guitar)
narrowPeak <- system.file(
  "extdata", "m6A_hg19_1000peaks_macs2.narrowPeak", 
  package="Guitar")
# genomic features imported into named list
m6A_Bcell <- narrowPeaktoGRanges(narrowPeak) 
m6A_Bcell_1 <- m6A_Bcell[1:300]
m6A_Bcell_2 <- m6A_Bcell[301:600]
m6A_Bcell_3 <- m6A_Bcell[601:900]
feature_hg19 <- list(m6A_Bcell_1, m6A_Bcell_2, m6A_Bcell_3) 
names(feature_hg19) <- c("Group_1","Group_2","Group_3")


###################################################
### code chunk number 3: rep_dendro_1_plot (eval = FALSE)
###################################################
## count <- Guitar2Plot(feature_hg19, 
##            genome="hg19", 
##            saveToPDFprefix = "example")


###################################################
### code chunk number 4: Txdb
###################################################
txdb_file <- system.file("extdata", "hg19_toy.sqlite", 
                         package="Guitar")
txdb <- loadDb(txdb_file)
# Or use makeTxDbFromUCSC() to download TxDb from internet
# txdb <- makeTxDbFromUCSC(genome="hg19")


###################################################
### code chunk number 5: guitar coordiantes
###################################################
gc_txdb <- makeGuitarCoordsFromTxDb(txdb, noBins=100)


###################################################
### code chunk number 6: Guitarcoordiantes
###################################################
Guitar2Plot(gfeatures = feature_hg19, 
           GuitarCoordsFromTxDb = gc_txdb,
           saveToPDFprefix = "example2")


###################################################
### code chunk number 7: withDNA2
###################################################
Guitar2Plot(gfeatures = feature_hg19, 
           GuitarCoordsFromTxDb = gc_txdb,
           rescaleComponent=FALSE)


###################################################
### code chunk number 8: withDNA3
###################################################
Guitar2Plot(gfeatures = feature_hg19, 
           GuitarCoordsFromTxDb = gc_txdb,
           includeNeighborDNA =TRUE,
           rescaleComponent=FALSE,
           fill=TRUE)


###################################################
### code chunk number 9: mm10_example
###################################################
# import different data formats into a named list object.
# These genomic features are using mm10 genome assembly
bed3=system.file("extdata", "H3K4me3_mm10_1000peaks.bed", package="Guitar")
bed12=system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed", package="Guitar")
bam1=system.file("extdata", "866991_part.bam", package="Guitar")
bam2=system.file("extdata", "866992_part.bam", package="Guitar")
H3K4me3 <- import.bed(bed3) # bed3 imported as GRanges
m6A_HepG2 <- BED12toGRangesList(bed12) # bed12 imported as GRangesList
MeRIP_IP <- readGAlignments(bam1) # bam imported as GAlignments
MeRIP_Input <- readGAlignments(bam2) # bam imported as GAlignments
feature_mm10 <- list(H3K4me3,m6A_HepG2,MeRIP_IP,MeRIP_Input) 
names(feature_mm10) <- c("H3K4me3", "m6A_HepG2", "IP","Input")

# Build Guitar Coordinates
txdb_file <- system.file("extdata", "mm10_toy.sqlite", 
                         package="Guitar")
mm10_toy_txdb <- loadDb(txdb_file)
gc_mm10_txdb <- makeGuitarCoordsFromTxDb(mm10_toy_txdb, noBins=30)

# Guitar Plot
Guitar2Plot(gfeatures = feature_mm10, 
           GuitarCoordsFromTxDb = gc_mm10_txdb,
           rescaleComponent=FALSE)


###################################################
### code chunk number 10: bsseq (eval = FALSE)
###################################################
## f <-system.file("extdata", "DNAm5C_mm10_10000sites_bismark.cov", package="Guitar") 
## a <- read.table(f,header=FALSE,sep="\t")
## q <- a[[5]]
## size <- a[[5]]+a[[6]]
## mean_methy <- sum(q)/sum(size)
## d0 <- (pbinom(q, size, mean_methy, lower.tail = TRUE, log.p = FALSE) < 0.05)
## d1 <- (pbinom(q, size, mean_methy, lower.tail = FALSE, log.p = FALSE) < 0.05)  
## GR <- GRanges(seqnames = a[,1], ranges = IRanges(start=a[,2], width=1))
## peak <- list(GR[d1],GR[d1+d0 == 0],GR[d0])
## names(peak) <- c( "higher", "unsure", "lower")
## 
## TxDb.Mmusculus.UCSC.mm10.knownGene <- makeTxDbFromUCSC(genome="mm10")
## gc_mm10_txdb <- makeGuitarCoordsFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene)
## GuitarPlot(gfeatures = peak,
##            GuitarCoordsFromTxDb = gc_mm10_txdb,
##            includeNeighborDNA =TRUE,
##            fill=TRUE)


###################################################
### code chunk number 11: session
###################################################
head(gc_txdb)
mcols(gc_txdb)


###################################################
### code chunk number 12: lcp3
###################################################
GuitarCoords <- reduce(gc_txdb)  # extract the coordinates
gcl <- list(GuitarCoords) # put into a list
GuitarPlot(gfeatures = gcl,
           GuitarCoordsFromTxDb = gc_txdb,
           rescaleComponent=TRUE)


###################################################
### code chunk number 13: lcp2
###################################################
GuitarCoords <- reduce(gc_txdb)  # extract the coordinates
gcl <- list(GuitarCoords) # put into a list
GuitarPlot(gfeatures = gcl,
           GuitarCoordsFromTxDb = gc_txdb,
           rescaleComponent=FALSE)


###################################################
### code chunk number 14: comp
###################################################
GuitarCoords <- gc_txdb
type <- paste(mcols(GuitarCoords)$comp,mcols(GuitarCoords)$category)
key <- unique(type)
landmark <- list(1,2,3,4,5,6,7,8)
names(landmark) <- key  
for (i in 1:length(key)) {
  landmark[[i]] <- GuitarCoords[type==key[i]]
}


###################################################
### code chunk number 15: mcp
###################################################
GuitarPlot(gfeatures=landmark[1:5], 
           GuitarCoordsFromTxDb = gc_txdb,
           includeNeighborDNA =TRUE,
           rescaleComponent=FALSE)


###################################################
### code chunk number 16: lcp
###################################################
GuitarPlot(gfeatures=landmark[6:8], 
           GuitarCoordsFromTxDb = gc_txdb,
           includeNeighborDNA =TRUE,
           rescaleComponent=FALSE)


###################################################
### code chunk number 17: session
###################################################
sessionInfo()


