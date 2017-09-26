


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
## count <- GuitarPlot(feature_hg19, 
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
# GuitarPlot(gfeatures = feature_hg19, 
#            GuitarCoordsFromTxDb = gc_txdb,
#            saveToPDFprefix = "example2")
gfeatures = feature_hg19
GuitarCoordsFromTxDb=gc_txdb

Guitar2Plot(gfeatures= feature_hg19,
           GuitarCoordsFromTxDb=gc_txdb,
           txdb=txdb,
           genome=NA,
           noBins=10,
           transcript = 'mRNA',
           saveToPDFprefix='example',
           returnCount=FALSE,
           includeNeighborDNA=FALSE,
           maximalFeatureAmbiguity=5,
           rescaleComponent=TRUE,
           CI_interval=c(0.025,0.975),
           adjust=1)