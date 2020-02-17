# SPI1 mm9 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58128
annotatePeaks.pl pu1_GSE58128_liftover2mm10.bed mm10 -annStats pu1_GSE58128_liftover2mm10.annStats > pu1_GSE58128_liftover2mm10.anno 

# TAL1 mm10 https://www.encodeproject.org/search/?searchTerm=ENCSR000DIA
annotatePeaks.pl TAL1_MEL.bed mm10 -annStats TAL1_MEL.annStats > TAL1_MEL.anno 

more pu1_GSE58128_liftover2mm10.anno|grep "promoter-TSS"|cut -f16|more|sort|uniq > pu1_GSE58128_Targetgenes.txt
more TAL1_MEL.anno|grep "promoter-TSS"|cut -f16|more|sort|uniq > TAL1_MEL.txt
#####################################################################################################################################
#####################################################################################################################################
exp <- read.csv("Tip60_RNAseq_triplicates.csv")
#exp$Gene <- toupper(exp$Gene)
exp <- exp[,c(1,3)]

act <- read.table("actRatio.txt",sep="\t",header=T)
#act$Gene <- toupper(act$Gene)

tal1<-read.table("TAL1_MEL.txt")
pu1<-read.table("pu1_GSE58128_Targetgenes.txt")

ix = match(act$Gene,exp$Gene)
