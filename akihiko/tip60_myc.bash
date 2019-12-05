# getting log2 ratio

scp -P 60035 root@172.18.149.78:/root/akihiko/oriFiles/*bw ./

/home/rtm/anaconda2/bin/conda install -c bioconda deeptools

# 
/home/rtm/anaconda2/bin/bigwigCompare -b1 actH2AZ_neg_R1.bw -b2 H2AZ_neg_R1.bw --operation log2 -bs 20 -p max -o NEG_ac_h2_log2ratio.bw -of bigwig
/home/rtm/anaconda2/bin/bigwigCompare -b1 actH2AZ_pos_R1.bw -b2 H2AZ_pos_R1.bw --operation log2 -bs 20 -p max -o POS_ac_h2_log2ratio.bw -of bigwig

/home/rtm/anaconda2/bin/bigwigCompare -b1 actH2AZ_neg_R1.bw -b2 H2AZ_neg_R1.bw --operation ratio -bs 20 -p max -o NEG_ac_h2_ratio.bw -of bigwig
/home/rtm/anaconda2/bin/bigwigCompare -b1 actH2AZ_pos_R1.bw -b2 H2AZ_pos_R1.bw --operation ratio -bs 20 -p max -o POS_ac_h2_ratio.bw -of bigwig

#bigwigCompare -b1 ../oriFiles/actH2AZ_neg_R1.bw -b2 ../oriFiles/H2AZ_neg_R1.bw --operation log2 -bs 20 -p max -o NEG_ac_h2_log2ratio.bw -of bigwig
#bigwigCompare -b1 ../oriFiles/actH2AZ_pos_R1.bw -b2 ../oriFiles/H2AZ_pos_R1.bw --operation log2 -bs 20 -p max -o POS_ac_h2_log2ratio.bw

#bigwigCompare -b1 ../oriFiles/actH2AZ_neg_R1.bw -b2 ../oriFiles/H2AZ_neg_R1.bw --operation ratio -bs 20 -p max -o NEG_ac_h2_ratio.bw
#bigwigCompare -b1 ../oriFiles/actH2AZ_pos_R1.bw -b2 ../oriFiles/H2AZ_pos_R1.bw --operation ratio -bs 20 -p max -o POS_ac_h2_ratio.bw


############################################################################################################################
bedtools intersect -a /root/akihiko/oriFiles/TIP60_peaks.bed -b /root/akihiko/tfregulome/tfregulomeMYC.bed > tip60_myc.bed # 3712
bedtools intersect -a /root/akihiko/oriFiles/TIP60_peaks.bed -b /root/akihiko/tfregulome/tfregulomeMYC.bed  -v > tip60_-_myc.bed # 659

annotatePeaks.pl tip60_myc.bed mm10 -annStats tip60_myc.annStats > tip60_myc.anno 
annotatePeaks.pl tip60_-_myc.bed mm10 -annStats tip60_-_myc.annStats > tip60_-_myc.anno 

 grep "Intergenic" tip60_myc.anno|cut -f2,3,4,10,16 > tip60_myc_intergenic.bed
 grep "promoter-TSS" tip60_myc.anno|cut -f2,3,4,10,16 > tip60_myc_promoter.bed

 grep "Intergenic" tip60_-_myc.anno|cut -f2,3,4,10,16 > tip60_-_myc_intergenic.bed
 grep "promoter-TSS" tip60_-_myc.anno|cut -f2,3,4,10,16 > tip60_-_myc_promoter.bed
############################################################################################################################



 ############################################################################################################################
 
 # pie charts
#  TIP60 + MYC
pdf("pie_dist_tip60_+_myc.pdf",width=12)
res=read.table(pipe("more tip60_myc.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>20]
pie(sort(tdown), main=,cex=1.7)
title("TIP60 peaks overlapping Myc peaks", cex.main=2)
dev.off()

#  TIP60 - MYC
pdf("pie_dist_tip60_-_myc.pdf",width=12)
res=read.table(pipe("more tip60_-_myc.annStats|cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>100]
pie(sort(tdown), main=,cex=1.7)
title("TIP60 peaks not overlapping Myc peaks", cex.main=2)
dev.off()
############################################################################################################################
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library("org.Hs.eg.db")
library(ReactomePA)
library(reactome.db) 

intergenic_acH2AZ_plus_MYC<-read.table("tip60_myc_intergenic.bed",sep="\t")
intergenic_acH2AZ_plus_MYC<-as.character(intergenic_acH2AZ_plus_MYC[,5])

promoter_acH2AZ_plus_MYC<-read.table("tip60_myc_promoter.bed",sep="\t")
promoter_acH2AZ_plus_MYC<-as.character(promoter_acH2AZ_plus_MYC[,5])

intergenic_acH2AZ_neg_MYC<-read.table("tip60_-_myc_intergenic.bed",sep="\t")
intergenic_acH2AZ_neg_MYC<-as.character(intergenic_acH2AZ_neg_MYC[,5])

promoter_acH2AZ_neg_MYC<-read.table("tip60_-_myc_promoter.bed",sep="\t")
promoter_acH2AZ_neg_MYC<-as.character(promoter_acH2AZ_neg_MYC[,5])

geneEntrez <- list(intergenic_acH2AZ_plus_MYC = intergenic_acH2AZ_plus_MYC,
    promoter_acH2AZ_plus_MYC = promoter_acH2AZ_plus_MYC,
    intergenic_acH2AZ_neg_MYC = intergenic_acH2AZ_neg_MYC,
    promoter_acH2AZ_neg_MYC = promoter_acH2AZ_neg_MYC)

names(geneEntrez) <- c("Dis.TIP60+Myc","Pro.TIP60+Myc",
                       "Dis.TIP60-Myc","Pro.TIP60-Myc")

x=compareCluster(geneEntrez, fun='enrichGO',
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP")
pdf("dotplot.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun='enrichPathway',
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP")
pdf("dotplot_enrichPathway.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun='enrichKEGG',
                 OrgDb         = org.Mm.eg.db,
                 ont           = "BP")
pdf("dotplot_enrichKEGG.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()


