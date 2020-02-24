# SPI1 mm9 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58128
annotatePeaks.pl pu1_GSE58128_liftover2mm10.bed mm10 -annStats pu1_GSE58128_liftover2mm10.annStats > pu1_GSE58128_liftover2mm10.anno 

# TAL1 mm10 https://www.encodeproject.org/search/?searchTerm=ENCSR000DIA
annotatePeaks.pl TAL1_MEL.bed mm10 -annStats TAL1_MEL.annStats > TAL1_MEL.anno 

more pu1_GSE58128_liftover2mm10.anno|grep "promoter-TSS"|cut -f16|more|sort|uniq > pu1_GSE58128_Targetgenes.txt
more TAL1_MEL.anno|grep "promoter-TSS"|cut -f16|more|sort|uniq > TAL1_MEL.txt
#####################################################################################################################################
#####################################################################################################################################
library(gplots)
library(factoextra)
library(RColorBrewer)

exp <- read.csv("Tip60_RNAseq_triplicates.csv")
exp$Gene <- toupper(exp$Gene)
exp <- exp[,c(1,3)]

act <- read.table("actRatio.txt",sep="\t",header=T)
#act$Gene <- toupper(act$Gene)

tal1<-read.table("TAL1_MEL.txt")
pu1<-read.table("pu1_GSE58128_Targetgenes.txt")
tal1 <- toupper(tal1[,1])
pu1 <- toupper(pu1[,1])

ix = match(act$Gene,exp$Gene)
act.m <- act[!is.na(ix),]
exp.m <- exp[ix[!is.na(ix)],]

dataf <- data.frame(act = act.m$log2RatioChange, exp = exp.m$log2FoldChange)

pdf("TAL1_targets.pdf")
plot(dataf$act,dataf$exp,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Tal1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )

  points(dataf$act[act$Gene %in% tal1],
       dataf$exp[act$Gene %in% tal1],
      col="red",pch=20)
dev.off()


pdf("PU1_targets.pdf")
plot(dataf$act,dataf$exp,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="PU1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )

  points(dataf$act[act$Gene %in% pu1],
       dataf$exp[act$Gene %in% pu1],
      col="red",pch=20)
dev.off()

#####################################################################################################################################
library(gplots)
library(factoextra)
library(RColorBrewer)
exp <- read.table("Tip60_actRatio_log2FC.txt",sep="\t",header = T)

tal1<-read.table("TAL1_MEL.txt")
pu1<-read.table("pu1_GSE58128_Targetgenes.txt")
tal1 <- toupper(tal1[,1])
pu1 <- toupper(pu1[,1])

pdf("TAL1_targets_chipseq.pdf")
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Tal1 targets (ChIP-Seq)",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% tal1],
       exp$log2FoldChange[exp$Gene %in% tal1],
      col="red",pch=20)

abline(v=0)
abline(h=0)
dev.off()

pdf("PU1_targets_chipseq.pdf")
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="PU1 targets (ChIP-Seq)",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% pu1],
       exp$log2FoldChange[exp$Gene %in% pu1],
      col="red",pch=20)

abline(v=0)
abline(h=0)
dev.off()

spi1<-read.table("Spi1_targets.mouse.tsv")
gata2<-read.table("Gata2_targets.mouse.tsv")
spi1 <- toupper(spi1[,2])
gata2 <- toupper(gata2[,2])

pdf("spi1_targets_trrust.pdf")
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="PU1 targets (TRRUST)",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% spi1],
       exp$log2FoldChange[exp$Gene %in% spi1],
      col="red",pch=20)

abline(v=0)
abline(h=0)
dev.off()

pdf("GATA2_targets_trrust.pdf")
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="GATA2 targets (TRRUST)",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% gata2],
       exp$log2FoldChange[exp$Gene %in% gata2],
      col="red",pch=20)

abline(v=0)
abline(h=0)
dev.off()


pdf("tf_targets.pdf",height=2.5)
par(mfrow=c(1,3))
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Myc targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$group=="myc"],
       exp$log2FoldChange[exp$group=="myc"],
      col=alpha("red",.5),pch=20)

abline(v=0)
abline(h=0)

plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Tal1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% tal1],
       exp$log2FoldChange[exp$Gene %in% tal1],
      col=alpha("darkgreen",.5),pch=20)

abline(v=0)
abline(h=0)

plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Pu.1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% pu1],
       exp$log2FoldChange[exp$Gene %in% pu1],
      col=alpha("blue",.5),pch=20)

abline(v=0)
abline(h=0)
dev.off()

pdf("legends_targets.pdf")
plot.new()
legend("center",  legend=c("Myc targets","Tal1 targets","Pu.1 targets"), fill=c('red','darkgreen','blue'), bty = "n")
dev.off()
############################################

legend("topright", paste(length(which(dds_res$log2FoldChange>1 & dds_res$padj<0.05))), bty="n") 




pdf("tf_targets.pdf",height=2.5)
par(mfrow=c(1,3))
plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Myc targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$group=="myc"],
       exp$log2FoldChange[exp$group=="myc"],
      col=alpha("red",.5),pch=20)
myc_n = 
legend("topright",length(which(exp$log2RatioChange[exp$group=="myc"]>0 & exp$log2FoldChange[exp$group=="myc"]>0))/myc_n, bty="n") 
legend("topright",length(which(exp$log2RatioChange[exp$group=="myc"]<0 & exp$log2FoldChange[exp$group=="myc"]<0))/myc_n, bty="n") 
legend("topright",length(which(exp$log2RatioChange[exp$group=="myc"]<0 & exp$log2FoldChange[exp$group=="myc"]>0))/myc_n, bty="n") 
legend("topright",length(which(exp$log2RatioChange[exp$group=="myc"]>0 & exp$log2FoldChange[exp$group=="myc"]<0))/myc_n, bty="n") 

abline(v=0)
abline(h=0)

plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Tal1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% tal1],
       exp$log2FoldChange[exp$Gene %in% tal1],
      col=alpha("darkgreen",.5),pch=20)

abline(v=0)
abline(h=0)

plot(exp$log2RatioChange,exp$log2FoldChange,xlab=expression('Log'[2]*paste('acH2AZ/H2AZ Fold Change')),main="Pu.1 targets",
              ylab=expression('Log'[2]*paste('RNA expression Fold Change')),col=alpha("grey",.5),pch=20 )
  points(exp$log2RatioChange[exp$Gene %in% pu1],
       exp$log2FoldChange[exp$Gene %in% pu1],
      col=alpha("blue",.5),pch=20)

abline(v=0)
abline(h=0)
dev.off()
