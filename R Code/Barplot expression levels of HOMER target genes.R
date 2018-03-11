##### Set experiment name
Title <- "Rev_vs_Irr.DEGs.FDR0.1.FC2.Irr."

##### Read in
## Target genes
require(data.table)
print(filenames <- list.files(pattern="*.target_genes$"))
Target_gene.list <- lapply(filenames, function(files) {
  read.delim(files, header= T, sep = '\t', quote="\"'")
})

## change the names of each df in list
names(Target_gene.list) <- gsub("Rev_vs_Irr.DEGs.FDR0.1.FC2.Irr.", "", filenames)

## mRNA expression data
mRNA_expression.File <- "/Users/jchap12/Google Drive/PAH/Diederik/Rev_vs_Irr.full.csv"
mRNA_expression <- read.delim(file= mRNA_expression.File, header= T, sep = ',', quote="\"'")

##### Create a function & plot rlog normalized counts for single genes on log scale ####
singleGenesofInterest <- as.character(df$Name) #make boxplots for all DEGs

## make a barplot function to plot with ggplot2
df <- Target_gene.list$ID2.target_genes #for testing only

PlotTargetFCs <- function(df) {
  title <- deparse(substitute(df)) #substitute is weird (must use it before altering df)
  ## remove the duplicated rows
  df_dedup <- df[!duplicated(as.character(df$Ensembl)),]
  ## join the df with mRNA expression levels
  df2 <- merge(x= df_dedup, y= mRNA_expression, by.x = "Ensembl", by.y= "Ensembl_gene")
  ## change duplicates by appending ".num" string to end
  df2$Name = make.names(df2$Name, unique=T)
  ## sort by fold change
  df2 <- df2[order(df2$log2FoldChange, decreasing = T), ]
  ## set the genes as factors, so that ggplot won't plot genenames alphabetically
  df2$Name <- factor(df2$Name, levels = df2$Name)
  ## make barplot of log2FC vs condition
  df <- df2 
  require(ggplot2)
  require(export)
  textSize  <- element_text(size= 14)
  plota <- ggplot(data= df, aes(y= (log2FoldChange * -1), x= factor(Name))) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle("DMC Distance from TSS") +
    theme(plot.title= element_text(size= 14, face= "bold"),
          axis.text= textSize, legend.text= textSize,
          legend.title= textSize, axis.title= textSize) +
    geom_text(aes(label=round(padj, digits = 2)), position=position_dodge(width=1), hjust=-0.5) +
    labs(title=title, x="Gene", y = "Irr vs Rev mRNA (log2FC)") + coord_flip()
  print(plota)
  graph2ppt(file="Homer_Target_gene.barplots.pptx", width=7, height=7, append=T)
  dev.off()
}

### list all target_genes
list2env(Target_gene.list,.GlobalEnv)
print(testfiles <- names(Target_gene.list))

### generate plots
PlotTargetFCs(Target_gene.list$ID2.target_genes)
PlotTargetFCs(Target_gene.list$NFIC.target_genes)
PlotTargetFCs(Target_gene.list$P53.target_genes)
PlotTargetFCs(Target_gene.list$RUNX.target_genes)
PlotTargetFCs(Target_gene.list$SPIB.target_genes)
PlotTargetFCs(Target_gene.list$TBX1.target_genes)
