##### Histogram ####
df <- DMCs.PAH_vs_con
qplot(df$dist.to.feature,
      geom="histogram", binwidth = 10000,
      main = "Distribution of DMCs relative to TSS", col="Green", 
      xlab="distToTSS") +
  theme(plot.title= element_text(size= 14, face= "bold"),
        axis.text= element_text(size= 14), legend.text= element_text(size= 14),
        legend.title= element_text(size= 14), axis.title= element_text(size= 14))

## Export ppt
graph2ppt(file=paste(title,".distribution.ppt",sep=''), width=7, height=7, append=T)