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
