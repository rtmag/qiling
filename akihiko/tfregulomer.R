library(TFregulomeR)

myc_mouse = dataBrowser(species="mouse", tf="myc")

 for (i in  myc_mouse$ID){
    peak_i = loadPeaks(i, includeMotifOnly=T)
    write.table(peak_i, paste0(i,"_peaks_with_motif.bed", sep=""), col.names= FALSE,row.names=FALSE, quote=FALSE)
 }
 
 
system("cat GTRD*MYC_peaks_with_motif.bed||perl -pe 's/\n/HOLAMUNDO/g'|perl -pe 's/\s/\t/g'|
        perl -pe 's/HOLAMUNDO/\n/g'|cut -f1,2,3|awk -F\"\t\" {\'print $1\"\t\"$2\'}")
        
cat GTRD*MYC_peaks_with_motif.bed|perl -pe 's/\n/HOLAMUNDO/g'|perl -pe 's/\s/\t/g'| \
perl -pe 's/HOLAMUNDO/\n/g'|cut -f1,2,3|awk -F"\t" {'print $1"\t"$2-50"\t"$3+50'}|sort -k1,1 -k2,2n| \
bedtools merge -i - > tfregulomeMYC.bed

# 54992 MYC peaks total, 20,746 peaks,  34,246 myc peaks in 
###
bedtools intersect -a tfregulomeMYC.bed -b ../oriFiles/Ql_actH2AZ_neg_R1.bed > mycTFREGULOME_h2az_peaks.bed

annotatePeaks.pl mycTFREGULOME_h2az_peaks.bed mm10 -annStats mycTFREGULOME_h2az_peaks.annStats > mycTFREGULOME_h2az_peaks.anno 

 grep "Intergenic" mycTFREGULOME_h2az_peaks.anno|cut -f2,3,4,10,16 > Intergenic_mycTFREGULOME_h2az_peaks.bed
 grep "promoter-TSS" mycTFREGULOME_h2az_peaks.anno|cut -f2,3,4,10,16 > promoter_mycTFREGULOME_h2az_peaks.bed
###

bedtools intersect -a ../oriFiles/Ql_actH2AZ_neg_R1.bed -b tfregulomeMYC.bed -v  > ach2az_-_mycTFREGULOME_peaks.bed

annotatePeaks.pl ach2az_-_mycTFREGULOME_peaks.bed mm10 -annStats ach2az_-_mycTFREGULOME_peaks.annStats > ach2az_-_mycTFREGULOME_peaks.anno 
annotatePeaks.pl ../oriFiles/Ql_actH2AZ_neg_R1.bed mm10 -annStats ach2az_peaks.annStats > ach2az_peaks.anno 
##
 grep "Intergenic" ach2az_-_mycTFREGULOME_peaks.anno|cut -f2,3,4,10,16 > Intergenic_ach2az_-_mycTFREGULOME_peaks.bed
 grep "promoter-TSS" ach2az_-_mycTFREGULOME_peaks.anno|cut -f2,3,4,10,16 > promoter_ach2az_-_mycTFREGULOME_peaks.bed
###
 grep "Intergenic" ach2az_peaks.anno|cut -f2,3,4,10,16 > Intergenic_ach2az_peaks.bed
 grep "promoter-TSS" ach2az_peaks.anno|cut -f2,3,4,10,16 > promoter_ach2az_peaks.bed

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
library(TFregulomeR)
my_peak <- read.delim("Intergenic_mycTFREGULOME_h2az_peaks.bed", sep = "\t", header = FALSE)

mouse_TFBS = dataBrowser(species = "mouse") # or TFBSBrowser() before v1.2.0

intersectMatrix_common <- intersectPeakMatrix(user_peak_list_x = list(my_peak),
                                              peak_id_y = mouse_TFBS$ID, 
                                              motif_only_for_id_y = TRUE)

intersectMatrix_common <- readRDS("TFREGULOMER_MYC+H2AZ_intersectMatrix_common.rds")

intersectPeakMatrixResult.test<-intersectPeakMatrixResult(intersectMatrix_common,angle_of_matrix="x",return_intersection_matrix=TRUE)

intersectPeakMatrixResult.result <- intersectPeakMatrixResult.test$intersection_matrix

cobind <- rev(sort(intersectPeakMatrixResult.result[1,]))

cobind <- data.frame(dataset=names(cobind[1,]),cobinding=as.numeric(cobind[1,]))

op <- par(mar=c(4,17,4,2)) # the 10 allows the names.arg below the barplot
barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=TRUE,col="skyblue",xlim=c(0,40),xlab="Cobinding(%) of distal Myc peaks overlapping acH2AZ peaks")
abline(h=0)

pdf("akihiko_cobinding.pdf")
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylim=c(0,40),ylab="Cobinding factors(%) of distal Myc + acH2AZ peaks")
abline(h=0)
dev.off()
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
library(TFregulomeR)

cobinding_mypeak <- function(bedfile,outname){
   my_peak <- read.delim(bedfile, sep = "\t", header = FALSE)
   mouse_TFBS = dataBrowser(species = "mouse") # or TFBSBrowser() before v1.2.0

   intersectMatrix_common <- intersectPeakMatrix(user_peak_list_x = list(my_peak),
                                              peak_id_y = mouse_TFBS$ID, 
                                              motif_only_for_id_y = TRUE)

   intersectPeakMatrixResult.test<-intersectPeakMatrixResult(intersectMatrix_common,angle_of_matrix="x",return_intersection_matrix=TRUE)
   intersectPeakMatrixResult.result <- intersectPeakMatrixResult.test$intersection_matrix
   cobind <- rev(sort(intersectPeakMatrixResult.result[1,]))
   cobind <- data.frame(dataset=names(cobind[1,]),cobinding=as.numeric(cobind[1,]))
   write.csv(cobind,paste0(outname,".csv"))
   
   pdf(paste0(outname,".pdf"))
   op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for ",dim(my_peak)[1]," peaks") )
   abline(h=0)
   dev.off()
}

cobinding_mypeak("promoter_ach2az_peaks.bed","promoter_ach2az_cobinding")
cobinding_mypeak("Intergenic_ach2az_peaks.bed","distal_ach2az_cobinding")

cobinding_mypeak("promoter_ach2az_-_mycTFREGULOME_peaks.bed","promoter_ach2az_-_myc_cobinding")
cobinding_mypeak("Intergenic_ach2az_-_mycTFREGULOME_peaks.bed","distal_ach2az_-_myc_cobinding")

cobinding_mypeak("promoter_mycTFREGULOME_h2az_peaks.bed","promoter_ach2az_+_myc_cobinding")
cobinding_mypeak("Intergenic_mycTFREGULOME_h2az_peaks.bed","distal_ach2az_+_myc_cobinding")
####################


   pdf("distal_ach2az_+_myc_cobinding.pdf")
cobind<- read.csv("distal_ach2az_+_myc_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for distal acH2AZ+Myc",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()

   pdf("distal_ach2az_-_myc_cobinding.pdf")
cobind<- read.csv("distal_ach2az_-_myc_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for distal acH2AZ-Myc",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()

   pdf("distal_ach2az_cobinding.pdf")
cobind<- read.csv("distal_ach2az_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for distal acH2AZ",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()
########

   pdf("promoter_ach2az_+_myc_cobinding.pdf")
cobind<- read.csv("distal_ach2az_+_myc_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for promoter acH2AZ+Myc",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()

   pdf("promoter_ach2az_-_myc_cobinding.pdf")
cobind<- read.csv("distal_ach2az_-_myc_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for promoter acH2AZ-Myc",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()

   pdf("promoter_ach2az_cobinding.pdf")
cobind<- read.csv("distal_ach2az_cobinding.csv",stringsAsFactors=FALSE,row.names=1)
op <- par(mar=c(17,4,4,2)) # the 10 allows the names.arg below the barplot
   barplot( height=cobind[,2][1:20], names.arg=gsub("GTRD-EXP.+_MMU_","",cobind$dataset[1:20],perl=TRUE),border=NA,las=2,
        horiz=FALSE,col="skyblue",ylab=paste0("Cobinding factors(%) for promoter acH2AZ",dim(cobind)[1]," peaks"),
           ylim=c(0,60),yaxt='n' )
   abline(h=0)
axis(2, at=seq(0, 60, by=10), labels = FALSE)
text(y = seq(0, 60, by=10)+2,-1, labels = seq(0, 60, by=10), srt = 90, pos = 2, xpd = TRUE)
   dev.off()
