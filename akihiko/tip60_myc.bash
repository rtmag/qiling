# getting log2 ratio

scp -P 60035 root@172.18.149.78:/root/akihiko/oriFiles/*bw ./
scp -P 60035 root@172.18.149.78:/root/akihiko/tip60_myc/*intergenic.bed ./
scp -P 60035 root@172.18.149.78:/root/akihiko/tip60_myc/*promoter.bed ./

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
# AcH2
computeMatrix reference-point \
-S \
actH2AZ_neg_R1.bw \
actH2AZ_pos_R1.bw \
-R tip60_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ach2_int_+.mat

computeMatrix reference-point \
-S \
actH2AZ_neg_R1.bw \
actH2AZ_pos_R1.bw \
-R tip60_-_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ach2_int_-.mat

computeMatrix reference-point \
-S \
actH2AZ_neg_R1.bw \
actH2AZ_pos_R1.bw \
-R tip60_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ach2_pro_+.mat

computeMatrix reference-point \
-S \
actH2AZ_neg_R1.bw \
actH2AZ_pos_R1.bw \
-R tip60_-_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ach2_pro_-.mat

# H2

computeMatrix reference-point \
-S \
H2AZ_neg_R1.bw \
H2AZ_pos_R1.bw \
-R tip60_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out h2_int_+.mat

computeMatrix reference-point \
-S \
H2AZ_neg_R1.bw \
H2AZ_pos_R1.bw \
-R tip60_-_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out h2_int_-.mat

computeMatrix reference-point \
-S \
H2AZ_neg_R1.bw \
H2AZ_pos_R1.bw \
-R tip60_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out h2_pro_+.mat

computeMatrix reference-point \
-S \
H2AZ_neg_R1.bw \
H2AZ_pos_R1.bw \
-R tip60_-_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out h2_pro_-.mat

# log2ratio

computeMatrix reference-point \
-S \
NEG_ac_h2_log2ratio.bw \
POS_ac_h2_log2ratio.bw \
-R tip60_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out lr_int_+.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_log2ratio.bw \
POS_ac_h2_log2ratio.bw \
-R tip60_-_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out lr_int_-.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_log2ratio.bw \
POS_ac_h2_log2ratio.bw \
-R tip60_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out lr_pro_+.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_log2ratio.bw \
POS_ac_h2_log2ratio.bw \
-R tip60_-_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out lr_pro_-.mat

# ratio

computeMatrix reference-point \
-S \
NEG_ac_h2_ratio.bw \
POS_ac_h2_ratio.bw \
-R tip60_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out rt_int_+.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_ratio.bw \
POS_ac_h2_ratio.bw \
-R tip60_-_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out rt_int_-.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_ratio.bw \
POS_ac_h2_ratio.bw \
-R tip60_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out rt_pro_+.mat

computeMatrix reference-point \
-S \
NEG_ac_h2_ratio.bw \
POS_ac_h2_ratio.bw \
-R tip60_-_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out rt_pro_-.mat

# H3K27ac
computeMatrix reference-point \
-S \
H3K27ac_neg_R1.bw \
H3K27ac_pos_R1.bw \
-R tip60_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ac_int_+.mat

computeMatrix reference-point \
-S \
H3K27ac_neg_R1.bw \
H3K27ac_pos_R1.bw \
-R tip60_-_myc_intergenic.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ac_int_-.mat

computeMatrix reference-point \
-S \
H3K27ac_neg_R1.bw \
H3K27ac_pos_R1.bw \
-R tip60_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ac_pro_+.mat

computeMatrix reference-point \
-S \
H3K27ac_neg_R1.bw \
H3K27ac_pos_R1.bw \
-R tip60_-_myc_promoter.bed --referencePoint center \
--sortRegions descend -bs 20 -a 3000 -b 3000 -p max -out ac_pro_-.mat
 ############################################################################################################################
ach2_int_-.mat
ach2_int_+.mat
ach2_pro_-.mat
ach2_pro_+.mat

ac_int_-.mat
ac_int_+.mat
ac_pro_-.mat
ac_pro_+.mat

h2_int_-.mat
h2_int_+.mat
h2_pro_-.mat
h2_pro_+.mat

lr_int_-.mat
lr_int_+.mat
lr_pro_-.mat
lr_pro_+.mat

rt_int_-.mat
rt_int_+.mat
rt_pro_-.mat
rt_pro_+.mat
#
plotProfile  --colors 'Blue' 'Red' --regionsLabel "AcH2AZ" \
-m ach2_int_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Intergenic" \
-out ach2_int_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH2AZ" \
-m ach2_int_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Intergenic" \
-out ach2_int_-.pdf
#
plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH2AZ" \
-m ach2_pro_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Promoter" \
-out ach2_pro_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH2AZ" \
-m ach2_pro_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Promoter" \
-out ach2_pro_-.pdf
#######################
plotProfile  --colors 'Blue' 'Red' --regionsLabel "H2AZ" \
-m h2_int_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Intergenic" \
-out h2_int_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "H2AZ" \
-m h2_int_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Intergenic" \
-out h2_int_-.pdf
#
plotProfile --colors 'Blue' 'Red' --regionsLabel "H2AZ" \
-m h2_pro_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Promoter" \
-out h2_pro_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "H2AZ" \
-m h2_pro_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Promoter" \
-out h2_pro_-.pdf
#######################
plotProfile  --colors 'Blue' 'Red' --regionsLabel "AcH3K27" \
-m ac_int_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Intergenic" \
-out ac_int_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH3K27" \
-m ac_int_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Intergenic" \
-out ac_int_-.pdf
#
plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH3K27" \
-m ac_pro_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Promoter" \
-out ac_pro_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "AcH3K27" \
-m ac_pro_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Promoter" \
-out ac_pro_-.pdf
#
#######################
plotProfile  --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ log2Ratio" \
-m lr_int_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Intergenic" \
-out lr_int_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ log2Ratio" \
-m lr_int_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Intergenic" \
-out lr_int_-.pdf
#
plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ log2Ratio" \
-m lr_pro_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Promoter" \
-out lr_pro_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ log2Ratio" \
-m lr_pro_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Promoter" \
-out lr_pro_-.pdf
#######################
plotProfile  --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ Ratio" \
-m rt_int_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Intergenic" \
-out rt_int_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ Ratio" \
-m rt_int_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Intergenic" \
-out rt_int_-.pdf
#
plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ Ratio" \
-m rt_pro_+.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60+Myc Promoter" \
-out rt_pro_+.pdf

plotProfile --colors 'Blue' 'Red' --regionsLabel "acH2AZ/H2AZ Ratio" \
-m rt_pro_-.mat --perGroup --plotHeight 9 --plotWidth 9 --plotTitle "" \
 --samplesLabel "TIP60-FF" "TIP60-KO" --refPointLabel "Tip60-Myc Promoter" \
-out rt_pro_-.pdf

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
                 ont           = "BP",qvalueCutoff = 0.1,pvalueCutoff = 1)
pdf("dotplot.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun="enrichPathway", organism = "mouse",qvalueCutoff = 0.1,pvalueCutoff = 1)
pdf("dotplot_enrichPathway.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()

x=compareCluster(geneEntrez, fun="enrichKEGG", organism = "mouse", qvalueCutoff = 0.1,pvalueCutoff = 1)
pdf("dotplot_enrichKEGG.pdf",height=10,width=10)
dotplot(x, showCategory=15, includeAll=FALSE)
dev.off()


###############################################

genes <- read.table("45_distant_tip60_myc_genes.txt")
genes <-as.character(genes[,1])
exp <- read.csv("Tip60_RNAseq_triplicates.csv")
exp <-exp[exp$Gene %in% genes,]
exp.sig <- exp[,8:13]
rownames(exp.sig) <- as.character(exp[,1])
TIP60KO<-exp[,8:10] #tip60KO
TIP60WT<-exp[,11:13] #wt

x<-log2(cbind(rowMeans(TIP60WT)+1,rowMeans(TIP60KO)+1))
colnames(x)<-c("TIP60WT","TIP60KO")

pdf("tip60+myc_distal_genes_exp.pdf")
boxplot( x)
dev.off()


library(ComplexHeatmap)
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))
options(bitmapType="cairo")

exp.sig<-as.matrix(exp.sig)
exp.sig<-exp.sig[apply(exp.sig,1,sd)!=0,]

exp.sig.rm<-exp.sig-rowMeans(exp.sig)

pdf("TIP60+MYC_intergenic_associated_genes_expression.pdf",width=9)
Heatmap(exp.sig.rm,
show_row_names = TRUE,show_column_names = TRUE,name = "Expression",row_dend_reorder = TRUE, column_dend_reorder = TRUE,
clustering_distance_columns = "pearson", clustering_distance_rows = "pearson")
dev.off()
