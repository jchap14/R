##### PCA: PC2 vs PC1 ####
## assign a numeric matrix to mtrx
mtrx <- x

## specify plot aesthetics
color     <- "condition"
label     <- "condition"
shape     <- "Date"
mainTitle <- "Transcriptome PCA"
textSize  <- element_text(size= 14)

## Calculations
rv <- rowVars(mtrx)
## select # of genes to consider for PCA (change max 2 min to specify top 500)
select <- order(rv, decreasing= T)[seq_len(max(500, length(rv)))]
pca <- prcomp(t(mtrx[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], metadata)

## Plot
ggplot(d, aes_string(x= "PC1", y= "PC2", color=color, shape=shape, label=label)) +
  geom_point(size= 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() +
  geom_text(hjust=0, vjust=0, size= 5) + ggtitle(mainTitle) +
  theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
        legend.text= textSize, legend.title= textSize, axis.title= textSize, strip.text= textSize)
## Export to powerpoint
graph2ppt(file="Transcriptome.PCA.pptx", width=10, height=7, append=T)

##### PCA: PC3 vs PC2 ####
## Plot
ggplot(d, aes_string(x= "PC2", y= "PC3", color=color, shape=shape, label=label)) +
  geom_point(size= 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() +
  geom_text(hjust=0, vjust=0, size= 5) + ggtitle(mainTitle) +
  theme(plot.title= element_text(size= 14, face= "bold"), axis.text= textSize,
        legend.text= textSize, legend.title= textSize, axis.title=textSize, strip.text=textSize)

## Export to powerpoint
graph2ppt(file="Transcriptome.PCA.pptx", width=10, height=7, append=T)