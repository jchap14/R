##### Create a function & plot rlog normalized counts for single genes on log scale ####
singleGenesofInterest <- as.character(FDR10.EBSeq_results$Symbol) #make boxplots for all DEGs
# singleGenesofInterest <- c("gene", "gene2") # specify genes here
# Gene <- "blah" #test with single Gene here

## make a boxplot (or violin plot) function to plot with ggplot2
PlotGeneTPMs <- function(Gene) {
  title <- deparse(substitute(Gene)) #substitute is weird (must use it before altering df)
  df2 <- FDR10.EBSeq_results
  df <- data.frame(t(df2[df2$Symbol==Gene,c(colnames(TPM_matrix))]))
  df$condition <- as.character(metadata$condition)
  colnames(df)[1] <- "count"
  df$names <- rownames(df)
  textSize  <- element_text(size= 14)
  plota <- ggplot(df, aes(x= factor(condition), y= count, fill= factor(condition))) + ggtitle(Gene) +
    geom_boxplot(width = 0.2, outlier.shape= NA) + xlab("Cell Type") + ylab("TPM") +
    geom_point(shape=21, size=1, color="black",
               aes(y= count, colour= factor(condition)), position = "jitter") +
    theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
          legend.text= textSize, legend.title= textSize,
          axis.title= textSize, strip.text= textSize) #+ geom_text(aes(label=names),hjust=0, vjust=0)
  print(plota)
  graph2ppt(file="mRNA.boxplots.pptx", width=4, height=3, append=T)
  dev.off()
}

PlotGeneReadCounts <- function(Gene) {
  title <- deparse(substitute(Gene)) #substitute is weird (must use it before altering df)
  df2 <- FDR10.EBSeq_results
  df <- data.frame(t(df2[df2$Symbol==Gene,c(colnames(TPM_matrix))]))
  df$condition <- as.character(metadata$condition)
  colnames(df)[1] <- "count"
  df$names <- rownames(df)
  textSize  <- element_text(size= 14)
  plota <- ggplot(df, aes(x= factor(condition), y= count, fill= factor(condition))) + ggtitle(Gene) +
    geom_boxplot(width = 0.2, outlier.shape= NA) + xlab("Cell Type") + ylab("Gene Counts (log2)") +
    geom_point(shape=21, size=1, color="black",
               aes(y= count, colour= factor(condition)), position = "jitter") +
    theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
          legend.text= textSize, legend.title= textSize,
          axis.title= textSize, strip.text= textSize) #+ geom_text(aes(label=names),hjust=0, vjust=0)
  print(plota)
  graph2ppt(file="mRNA.boxplots.pptx", width=4, height=3, append=T)
  dev.off()
}

for(i in singleGenesofInterest) {
  PlotGeneTPMs(i)
  PlotGeneReadCounts(i)
}

