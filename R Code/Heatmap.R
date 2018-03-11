##### Transcriptome Heatmap ####
mtrx <- x #set mtrx again & ensure that it isn't the batch effect removed version
## select # of genes to consider for heatmap (change max 2 min to specify top 500)
rv <- rowVars(mtrx)
select <- order(rv, decreasing= T)[seq_len(min(500, length(rv)))]
mtrx <- mtrx[select, ]
## cluster rows by pearson, complete
hr <- hclust(as.dist(1-cor(t(mtrx), method='pearson')), method='complete')
hc <- hclust(as.dist(1-cor(mtrx, method='pearson')), method='complete')
## Heatmap2 w/ color bar. h & k modify cutree (k overrides h). Specify margins here as well.
mycl <- cutree(hr, h=max(hr$height)/3, k = 5);
mycol <- sample(rainbow(256)); mycol <- mycol[as.vector(mycl)]
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299)
## generate heatmap with rows clustered, but not columns
png('Transcriptome.heatmap_r.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(mtrx,Rowv=as.dendrogram(hr), Colv=NA, dendrogram= c("row"), col=my_palette, scale="row",
          density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## generate a heatmap with rows & columns clustered
png('Transcriptome.heatmap_rc.png', width= 7, height= 7, units= "in", res= 300, pointsize= 14)
heatmap.2(mtrx,Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), dendrogram= c("both"), col=my_palette,
          scale="row", density.info="none", trace="none", RowSideColors=mycol, margins= c(20, 5)) #margins c(height,width)
dev.off()
## plot & export column dendrograms
hc.dend <- hc %>% as.dendrogram #convert to dendrogram
plot(hc.dend, main = "Column Dendrogram")
graph2ppt(file="Transcriptome.dendrograms.pptx", width=10, height=7, append=T)
## plot & export row dendrograms
hr.dend <- hr %>% as.dendrogram #convert to dendrogram
hr.dend.cut <- cut(hr.dend, h= 1.9) #cut at desired height
plot(hr.dend.cut$upper, horiz = T, main = "Row Dendrogram")
graph2ppt(file="Transcriptome.dendrograms.pptx", width=10, height=7, append=T)