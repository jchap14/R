##### code to make a line plot of min overlap testing (improves the DiffBind built in) #####

##### in the script that is ran on SCG to get a data frame of overlap rate 
olap.rate <- dba.overlap(treatment,mode=DBA_OLAP_RATE)
overlapRate.df <- data.frame("sample_overlap" = c(1:length(olap.rate)), "peaks" = olap.rate)
## write out csv
write.csv(overlapRate.df, file="overlapRate.csv", row.names=F)

##### locally load & make overlap rate plot
overlapRate.df <- read.delim("overlapRate.csv", quote="\"'", sep = ",")
require(ggplot2)
textSize  <- element_text(size= 14)
oLapPlot <- ggplot(overlapRate.df, aes(y = peaks, x = sample_overlap)) + 
  geom_point(shape=16, fill="blue", color="darkred", size=3) + xlab("Overlap at least this many peaksets") +
  ylab("peak #") + geom_line() + scale_x_continuous(breaks=1:length(overlapRate.df$peaks)) +
  scale_y_continuous(breaks=c(250000, 500000, 750000, 1000000, 1250000, 1500000)) +
  theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
        legend.text= textSize, legend.title= textSize,
        axis.title= textSize, strip.text= textSize)
plot(oLapPlot)

## export to powerpoint
require(export)
graph2ppt(file=paste("Sample_peaks_overlap_plot.pptx",sep=''), width=10, height=7, append=T)
