\name{Guitar2Plot}
\alias{Guitar2Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Guitar2Plot
}
\description{
Plot the transcriptomic distribution of genomic features
}
\usage{
Guitar2Plot(gfeatures, GuitarCoordsFromTxDb = NA, 
  txdb = NA, genome = NA, 
  noBins = 10, saveToPDFprefix = NA, 
  returnCount = FALSE, includeNeighborDNA = FALSE,
  maximalFeatureAmbiguity=5, rescaleComponent=TRUE, adjust=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gfeatures}{
genomic features, which can be GRanges, GRangesList, GAlignment.
}
  \item{GuitarCoordsFromTxDb}{
Guitar coordiantes generated from makeGuitarCoordsFromTxDb object
}
  \item{txdb}{
If GuitarCoordsFromTxDb is not provided, you may optionally provide a TranscriptDb object, and the Guitar coordiante will be generated accordingly from it.
}
  \item{genome}{
If neither GuitarCoordsFromTxDb nortxdb  is not provided, you may optionally provide genome assembly number, gene annotation will be automatically downloaded from UCSC, and the Guitar coordiante will be generated accordingly from it.
}
  \item{noBins}{
The resolution of the analysis. Default: 100
}
  \item{saveToPDFprefix}{
Whether the figure generated should be saved into a PDF file.
}
  \item{returnCount}{
Whether the count information should be returned. This information is from countOverlaps between the provided gfeatures and Guitar coordiantes. It is not recommended to use this information.
}
  \item{includeNeighborDNA}{
Whether neighborhood DNA regions should be included in the Guitar2Plot.
}
  \item{maximalFeatureAmbiguity}{
Features that overlap with more than this number of transcripts will not be used in the analysis. If a genomic feature overlaps with multiple transcripts but the number doesn't exceed this threshold, its weight will be evenly divided. Default: 5
}
  \item{rescaleComponent}{
Whether calculate the relative location-standardized density or not. If set "TRUE", the density of tracks will be standardized with respect to each location. Default: 1
}
  \item{CI_interval}{
The confidence level used to plot the confidence interval. Default: c(0.025,0.975).
}
  \item{adjust}{
Adjustment of the bandwidth. Default: 1.
}
}
\details{
This function plots the transcriptomic distribution of genomic features. It is designed for a fast usage of the Guitar package without the needs to go into the details.
}
\value{
A figure showing the transcriptomic distribution of the genomic features will be generated. Post-editing with Adobe Illustrator or other graphic software is recommended.
}

\examples{

# read genomic features
narrowPeak <- system.file(
  "extdata", "m6A_hg19_1000peaks_macs2.narrowPeak", 
  package="Guitar")

# bam imported as GAlignments
m6A_Bcell <- narrowPeaktoGRanges(narrowPeak) 

# generate a list of genomic features
m6A_Bcell_1 <- m6A_Bcell[1:300]
m6A_Bcell_2 <- m6A_Bcell[301:600]
m6A_Bcell_3 <- m6A_Bcell[601:900]
feature_hg19 <- list(m6A_Bcell_1, m6A_Bcell_2, m6A_Bcell_3) 
names(feature_hg19) <- c("m6A_1","m6A_2","m6A_3")

# Make Guitar coordiantes
txdb_file <- system.file("extdata", "hg19_toy.sqlite", 
                         package="Guitar")
txdb <- loadDb(txdb_file)
gc_txdb <- makeGuitarCoordsFromTxDb(txdb,noBins =10)


# Plot
Guitar2Plot(feature_hg19, 
           GuitarCoordsFromTxDb = gc_txdb)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
