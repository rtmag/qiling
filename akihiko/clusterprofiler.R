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
##############################################################################################################################
# pie charts
#  H2 + MYC

pdf("pie_dist_acH2AZ+MYC.pdf",width=12)
res=read.table(pipe("more mycTFREGULOME_h2az_peaks.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>40]
pie(sort(tdown), main=,cex=1.7)
title("acH2az peaks overlapping Myc peaks", cex.main=2)
dev.off()


pdf("enrichment_acH2AZ+MYC.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more mycTFREGULOME_h2az_peaks.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,3]))
names(tdown) = res[,1]
tdown = tdown[as.numeric(as.character(res[,2]))>20]
barplot(sort(tdown),las=2,ylim=c(-4,6),ylab="Log2 Enrichment over random genomic background",col="lightblue3",
       cex.axis=2,cex.names=2)
abline(h=0)
dev.off()

#  H2 - MYC

pdf("pie_dist_acH2AZ-MYC.pdf",width=12)
res=read.table(pipe("more ach2az_-_mycTFREGULOME_peaks.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>250]
pie(sort(tdown), main=,cex=1.7)
title("acH2az peaks not overlapping Myc peaks", cex.main=2)
dev.off()


pdf("enrichment_acH2AZ-MYC.pdf.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more ach2az_-_mycTFREGULOME_peaks.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,3]))
names(tdown) = res[,1]
tdown = tdown[as.numeric(as.character(res[,2]))>250]
barplot(sort(tdown),las=2,ylim=c(-4,6),ylab="Log2 Enrichment over random genomic background",col="lightblue3",
        cex.axis=2,cex.names=2)
abline(h=0)
dev.off()





