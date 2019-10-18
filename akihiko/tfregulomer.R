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

bedtools intersect -a tfregulomeMYC.bed -b ../oriFiles/Ql_actH2AZ_neg_R1.bed > mycTFREGULOME_h2az_peaks.bed

annotatePeaks.pl mycTFREGULOME_h2az_peaks.bed mm10 -annStats mycTFREGULOME_h2az_peaks.annStats > mycTFREGULOME_h2az_peaks.anno 

 grep "Intergenic" mycTFREGULOME_h2az_peaks.anno|cut -f2,3,4,10,16 > mycTFREGULOME_h2az_peaks.bed
 grep "promoter-TSS" mycTFREGULOME_h2az_peaks.anno|cut -f2,3,4,10,16 > mycTFREGULOME_h2az_peaks.bed
 