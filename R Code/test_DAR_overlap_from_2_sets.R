##### Test for overlap in ATAC peaks from two sets ####
## NOTE: this will only work well to test overlap given the same 
## starting set of peaks. Otherwise nearly impossible

## set working directory
setwd("/Users/jchap12/Desktop")

## specify filenames
DARfile1 <- "together.txt"
DARfile2 <- "LS_vs_dF.txt"

## read in files
require(data.table)
DARs1 <- as.data.table(read.delim(DARfile1))
DARs2 <- as.data.table(read.delim(DARfile2))

## combine columns to make a SYMBOL_distToTss column
DARs1$SymDist <- paste(DARs1$SYMBOL, DARs1$distanceToTSS, sep= "_")
DARs2$SymDist <- paste(DARs2$SYMBOL, DARs2$distanceToTSS, sep= "_")

## merge the two to assess overlap
DARs.merged <- merge(DARs1, DARs2, by.x= "SymDist", by.y= "SymDist")
