###########################################################################################
########################## Take HOMER GO output and condense and graph it #################

## set Title to the "enriched terms" file in directory
enriched_terms <- list.files(pattern="*.enriched_terms.txt$")
Title <- gsub("\\.enriched_terms.txt", "", enriched_terms, perl=T)

##### Read in the functional_enrichment files from HOMER w/ header junk removed (made manually from HTML) ####
Terms.df <- read.delim(enriched_terms, quote="\"", header=T)

## all functional terms here
GO_Tree_terms <- as.vector(unique(Terms.df$GO.Tree))

## subset by 1 functional term at a time into a new df
GO_bp.df <- subset(Terms.df, Terms.df$GO.Tree == "biological process")
GO_bp.df <- GO_bp.df[,c("GO.ID","P.value","Term","ln.P.","GO.Tree","X..of.Genes.in.Term",
                        "X..of.Target.Genes.in.Term","X..of.Total.Genes",
                        "X..of.Target.Genes","Common.Genes")]

## write out df for use w Revigo online tool
write.table(GO_bp.df, paste(Title, ".RevigoIDs", sep=''), col.names= T, row.names=F, sep='\t')

##### Run Revigo online (column 1+2, "Small" option), save output ####

##### Read in Revigo output & filter ####
ReviGO <- read.delim("REVIGO.csv", quote="\"", header=T, sep=',')
## remove NULL (redundant) terms
ReviGO <- subset(ReviGO, plot_X != "null")
## merge w genes in GO terms
ReviGO.df <- merge(ReviGO, GO_bp.df, by.x= "term_ID" , by.y= "GO.ID")

##### remove terms who target >20% of tested terms 
number_tested_terms <- unique(ReviGO.df$X..of.Target.Genes) * 0.2
##or specify a number based on inspection
# number_tested_terms <- 71
ReviGO.df <- subset(ReviGO.df,X..of.Target.Genes.in.Term < number_tested_terms &
                           X..of.Target.Genes.in.Term > 2)

## subset interesting columns
ReviGO.df2 <- ReviGO.df[,c("term_ID","description","log10.p.value",
                           "GO.Tree","X..of.Genes.in.Term","X..of.Target.Genes.in.Term",
                           "Common.Genes")]
## sort by p-value
require("dplyr")
ReviGO.df2 <- arrange(ReviGO.df2, log10.p.value)
write.table(ReviGO.df2, paste(Title, ".ReviGO.GO_BP", sep=''), col.names= T, row.names=F, sep='\t')

##### Graph the results (plot more terms if going to remove manually) ####
require("ggplot2")
termNum <- nrow(ReviGO.df2) #"20"
df <- ReviGO.df2
size <- element_text(size= 14) #font size on plot

## make p-val positive, sort by p-value
df$log10.p.value <- df$log10.p.value * -1
df <- arrange(df, desc(log10.p.value))

##### Graph results (plot more terms if going to remove manually)
df$description <- factor(df$description, levels= df$description) #set X as factor preventing ABC order
a <- ggplot() + geom_bar(aes(y= log10.p.value, x= description), data= df, stat="identity") +
  coord_flip() + ggtitle("GO BP") + theme(plot.title= element_text(size= 14, face= "bold"),
                                          axis.text= size, legend.text= size,
                                          legend.title= size, axis.title= size) +
  geom_text(data=df, aes(x=description, y=log10.p.value, label=as.factor(X..of.Target.Genes.in.Term)),hjust=-0.5)
plot(a)
## export to powerpoint
require("export")
graph2ppt(file=paste(Title,".GO_and_Pathways.ppt",sep=''), width=10, height=9, append=T)

##### Subset & make graph for Reactome, KEGG, WikiPathways, BIOCYC, Pathway Interaction DB ####
Pathways.df <- subset(Terms.df, GO.Tree== "REACTOME pathways" | GO.Tree== "KEGG pathways" |
                        GO.Tree== "WikiPathways" | GO.Tree== "Pathway Interaction DB"| GO.Tree== "BIOCYC pathways") 
df <- Pathways.df[,c("GO.ID","P.value","Term", "ln.P.","GO.Tree","X..of.Genes.in.Term",
                     "X..of.Target.Genes.in.Term","X..of.Total.Genes","X..of.Target.Genes",
                     "Common.Genes")]

## remove non-unique terms here by keeping version w/ most genes
df <- arrange(df, desc(X..of.Target.Genes.in.Term), Term)
df <- subset(df, !duplicated(Term))

## remove terms who target >20% of tested terms
number_tested_terms <- max(unique(df$X..of.Target.Genes)) * 0.2
df <- subset(df, X..of.Target.Genes.in.Term < number_tested_terms &
               X..of.Target.Genes.in.Term > 2) #specify min here
## calculate -log10pVal & subset/reorg columns
df$log10pVal <- log10(df$P.value) * -1
df <- df[,c("Term", "log10pVal", "Common.Genes", "GO.Tree",
            "X..of.Target.Genes.in.Term", "X..of.Genes.in.Term",
            "GO.ID", "X..of.Total.Genes", "X..of.Target.Genes")]
## sort by # of Target genes in Term
df <- arrange(df, desc(X..of.Target.Genes.in.Term))
write.table(df, paste(Title, ".pathways", sep=''), col.names= T, row.names=F, sep='\t')

##### Plot Pathways results ####
df <- df[1:44,]  #plot more terms if going to remove manually
df <- arrange(df, desc(log10pVal)) ## sort by p-value
df$Term <- factor(df$Term, levels= df$Term) #set X as factor preventing ABC order
a <- ggplot() + geom_bar(aes(y= log10pVal, x= Term), data= df, stat="identity") +
  coord_flip() + ggtitle("Pathways") + theme(plot.title= element_text(size= 14, face= "bold"),
                                             axis.text= size, legend.text= size,
                                             legend.title= size, axis.title= size) +
  geom_text(data=df, aes(x=Term, y=log10pVal, label=as.factor(X..of.Target.Genes.in.Term)),hjust=-0.5)
plot(a)
graph2ppt(file=paste(Title,".GO_and_Pathways.ppt",sep=''), width=10, height=9, append=T)

