##### Test/plot correlation between DEGs: LSS vs ST & LSS vs OSS
Title  <- "Fold Change of Union DEGs" #set title for exported filenames
df     <- DEGs.union.df          #set dataframe of interest to df
dfy    <- df$log2FC.LS_vs_ST #set y feature
y_axis <- "log2FC.LS_vs_ST"    #set y axis
dfx    <- df$log2FC.LS_vs_OS #set x feature
x_axis <- "log2FC.LS_vs_OS"  #set x axis
Size   <- element_text(size= 14) #set text size for plot
## scatterplot with best fit line
ggplot(df, aes(x=dfx, y=dfy)) + geom_point(shape=1) + geom_smooth(method=lm) + xlab(x_axis)+
  ylab(y_axis) + ggtitle(Title) + theme(plot.title= element_text(size= 14, face= "bold"),
                                        axis.text= Size, legend.text= Size, legend.title= Size,
                                        axis.title= Size, strip.text= Size)
## fit linear model, get r2 (r is Pearson's correlation coefficient)
summary(lm(dfy ~ dfx))
print(paste("The Pearson correlation coefficient is ",round(cor(dfy, dfx),digits = 2), sep=''))
graph2ppt(file=paste(Title,".correlations.pptx",sep=''), width=3, height=3.15, append=T)
