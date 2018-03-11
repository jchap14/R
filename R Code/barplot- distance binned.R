##### Create distance binned barplots for distance to nearest TSS for each DAR set ####
library(ggplot2)
Name <- "LS_vs_ST.DARs"
counts  <- LS_vs_ST.DARs$distanceToTSS
## count the # of DMCs in distance bins & make a df
Bin3kb   <- sum(abs(counts) < 3000)
Bin10kb  <- sum(abs(counts) > 3000 & abs(counts) < 10000)
Bin100kb <- sum(abs(counts) > 10000 & abs(counts) < 100000)
Bin1MB   <- sum(abs(counts) > 100000 & abs(counts) < 1000000)
BinUpper <- sum(abs(counts) > 1000000)

CountTable <- as.data.frame(c("<3kb"= Bin3kb, "3-10kb"= Bin10kb,
                              "10-100kb"= Bin100kb, ">100kb"= Bin1MB,
                              ">1MB"= BinUpper))
colnames(CountTable) [1] <- "count"
CountTable$bins <- rownames(CountTable)
rownames(CountTable) <- NULL
## set the bins as a factor, so that ggplot won't plot genenames alphabetically
CountTable$bins <- factor(CountTable$bins, levels = CountTable$bins)
## plot the distance to nearest TSS for count
ggplot() + geom_bar(data= CountTable, aes(y= count, x= bins), stat="identity",
                    position="dodge") + ggtitle(paste(Name," Distance from TSS",sep='')) +
  theme(plot.title= element_text(size= 14, face= "bold"),
        axis.text= element_text(size= 14), legend.text= element_text(size= 14),
        legend.title= element_text(size= 14), axis.title= element_text(size= 14))
## export
graph2ppt(file=paste(NAME,".distribution.ppt",sep=''), width=5, height=5, append=T)
