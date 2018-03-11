##### try to convert Rat to Human homolog names using BiomaRt (eventually make a function) #### 
## select genes of interest by Ensembl term
genesOfInterest <- as.character(Rev_vs_Irr.DEGs.FDR0.1.FC2.Irr$Ensembl_gene)
## load biomaRt and setup
library("biomaRt")
ensembl=useMart("ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat   <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
## make df of attributes for output table
humanAttribute.DF <- listAttributes(human)
ratAttribute.DF <- listAttributes(rat)
## do the conversion
converted_genes <- getLDS(attributes = c("ensembl_gene_id","rgd_symbol", "description"),
                          filters= "ensembl_gene_id", values= genesOfInterest , mart = rat,
                          attributesL= c("ensembl_gene_id","hgnc_symbol","description"), martL= human)

