
bedtools intersect -a oriFiles/cMyc_peaks.bed -b oriFiles/Ql_actH2AZ_neg_R1.bed > myc_h2az_peaks.bed
# 3650 and 30009 > 3039
annotatePeaks.pl myc_h2az_peaks.bed mm10 -annStats myc_h2az_peaks.annStats > myc_h2az_peaks.anno 

 grep "Intergenic" myc_h2az_peaks.anno|cut -f2,3,4,10,16 > myc_h2az_peaks_intergenic.bed
 grep "promoter-TSS" myc_h2az_peaks.anno|cut -f2,3,4,10,16 > myc_h2az_peaks_promoter.bed
 
cat myc_h2az_peaks_intergenic.bed > myc_h2az_peaks.deeptoolsbed
echo "#Intergenic" >> myc_h2az_peaks.deeptoolsbed
 cat myc_h2az_peaks_promoter.bed >> myc_h2az_peaks.deeptoolsbed
echo "#Promoter" >> myc_h2az_peaks.deeptoolsbed

computeMatrix reference-point \
-S \
/root/akihiko/oriFiles/actH2AZ_neg_R1.bw \
/root/akihiko/oriFiles/actH2AZ_pos_R1.bw \
-R myc_h2az_peaks.deeptoolsbed --referencePoint center \
--sortRegions descend -bs 20 -a 5000 -b 5000 -p max -out actH2AZ_myc_h2az_peaks.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "Center" --colorMap Blues \
-m actH2AZ_myc_h2az_peaks.mat \
-out actH2AZ_myc_h2az_peaks.pdf 

computeMatrix reference-point \
-S \
/root/akihiko/oriFiles/H2AZ_neg_R1.bw \
/root/akihiko/oriFiles/H2AZ_pos_R1.bw \
-R myc_h2az_peaks.deeptoolsbed --referencePoint center \
--sortRegions descend -bs 20 -a 5000 -b 5000 -p max -out H2AZ_myc_h2az_peaks.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "Center" --colorMap Blues \
-m H2AZ_myc_h2az_peaks.mat \
-out H2AZ_myc_h2az_peaks.pdf 

computeMatrix reference-point \
-S \
/root/akihiko/oriFiles/H3K27ac_neg_R1.bw \
/root/akihiko/oriFiles/H3K27ac_pos_R1.bw \
-R myc_h2az_peaks.deeptoolsbed --referencePoint center \
--sortRegions descend -bs 20 -a 5000 -b 5000 -p max -out H3K27ac_myc_h2az_peaks.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "Center" --colorMap Blues \
-m H3K27ac_myc_h2az_peaks.mat \
-out H3K27ac_myc_h2az_peaks.pdf 

computeMatrix reference-point \
-S \
/root/akihiko/oriFiles/Tip60_neg_R1.bw \
/root/akihiko/oriFiles/Tip60_pos_R1.bw \
-R myc_h2az_peaks.deeptoolsbed --referencePoint center \
--sortRegions descend -bs 20 -a 5000 -b 5000 -p max -out Tip60_myc_h2az_peaks.mat

plotHeatmap --xAxisLabel "" --yAxisLabel "" --refPointLabel "Center" --colorMap Blues \
-m Tip60_myc_h2az_peaks.mat \
-out Tip60_myc_h2az_peaks.pdf 
