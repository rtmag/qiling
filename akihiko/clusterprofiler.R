library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("org.Mm.eg.db")


data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")


#######################################################################################

intergenic_h2_plus_myc<-read.table("ENTREZ_intergenic_ach2az_+_mycTFREGULOME_peaks.genes.txt")
intergenic_h2_plus_myc<-as.character(intergenic_h2_plus_myc[,1])

edo <- enrichDGN(intergenic_h2_plus_myc)

intergenic_h2_plus_myc %in% de

enrichGO(de)

gene.df <- bitr(intergenic_h2_plus_myc, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
        
ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))
