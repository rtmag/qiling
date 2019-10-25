library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("org.Mm.eg.db")
options(scipen=999)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats) ## for reordering the factor

#######################################################################################
# INT H2 + MYC
intergenic_h2_plus_myc<-read.table("ENTREZ_intergenic_ach2az_+_mycTFREGULOME_peaks.genes.txt")
intergenic_h2_plus_myc<-as.character(intergenic_h2_plus_myc[,1])

gene.df <- bitr(intergenic_h2_plus_myc, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
        
int_h2_and_myc <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

##########################################
# INT H2 - MYC
intergenic_h2_neg_myc<-read.table("ENTREZ_intergenic_ach2az_-_mycTFREGULOME_peaks.genes.txt")
intergenic_h2_neg_myc<-as.character(intergenic_h2_neg_myc[,1])

gene.df <- bitr(intergenic_h2_neg_myc, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
        
int_h2_neg_myc <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
#######################################################################################
# PROMO H2 + MYC
promoter_h2_plus_myc<-read.table("ENTREZ_promoter_ach2az_+_mycTFREGULOME_peaks.genes.txt")
promoter_h2_plus_myc<-as.character(promoter_h2_plus_myc[,1])

gene.df <- bitr(promoter_h2_plus_myc, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
        
prom_h2_and_myc <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

##########################################
# PROMO H2 - MYC
promoter_h2_neg_myc<-read.table("ENTREZ_promoter_ach2az_-_mycTFREGULOME_peaks.genes.txt")
promoter_h2_neg_myc<-as.character(promoter_h2_neg_myc[,1])

gene.df <- bitr(promoter_h2_neg_myc, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
        
prom_h2_neg_myc <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)



head(summary(ego2))rowPercentage
dotplot(ego2,type="bar", showCategory=20) + ggtitle("dotplot for ORA")
ddf<-as.data.frame(ego2)


gcSample


kk <- enrichGO(gene         = gcSample,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

x=compareCluster(gcSample, fun='enrichGO')
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# compareCluster
intergenic_acH2AZ_plus_MYC<-read.table("ENTREZ_intergenic_ach2az_+_mycTFREGULOME_peaks.genes.txt")
intergenic_acH2AZ_plus_MYC<-as.character(intergenic_acH2AZ_plus_MYC[,1])

promoter_acH2AZ_plus_MYC<-read.table("ENTREZ_promoter_ach2az_+_mycTFREGULOME_peaks.genes.txt")
promoter_acH2AZ_plus_MYC<-as.character(promoter_acH2AZ_plus_MYC[,1])

intergenic_acH2AZ_neg_MYC<-read.table("ENTREZ_intergenic_ach2az_-_mycTFREGULOME_peaks.genes.txt")
intergenic_acH2AZ_neg_MYC<-as.character(intergenic_acH2AZ_neg_MYC[,1])

promoter_acH2AZ_neg_MYC<-read.table("ENTREZ_promoter_ach2az_-_mycTFREGULOME_peaks.genes.txt")
promoter_acH2AZ_neg_MYC<-as.character(promoter_acH2AZ_neg_MYC[,1])

geneEntrez <- list(intergenic_acH2AZ_plus_MYC = intergenic_acH2AZ_plus_MYC,
    promoter_acH2AZ_plus_MYC = promoter_acH2AZ_plus_MYC,
    intergenic_acH2AZ_neg_MYC = intergenic_acH2AZ_neg_MYC,
    promoter_acH2AZ_neg_MYC = promoter_acH2AZ_neg_MYC)

names(geneEntrez) <- c("Dis.H2+Myc","Pro.H2+Myc",
                       "Dis.H2-Myc","Pro.H2-Myc")

x=compareCluster(geneEntrez, fun='enrichGO',
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP")
pdf("dotplot.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()










