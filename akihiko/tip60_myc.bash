# getting log2 ratio

bigwigCompare -b1 ../oriFiles/actH2AZ_neg_R1.bw -b2 ../oriFiles/H2AZ_neg_R1.bw --operation log2 -bs 20 -p max -o NEG_ac_h2_log2ratio.bw
bigwigCompare -b1 ../oriFiles/actH2AZ_pos_R1.bw -b2 ../oriFiles/H2AZ_pos_R1.bw --operation log2 -bs 20 -p max -o POS_ac_h2_log2ratio.bw

bigwigCompare -b1 ../oriFiles/actH2AZ_neg_R1.bw -b2 ../oriFiles/H2AZ_neg_R1.bw --operation ratio -bs 20 -p max -o NEG_ac_h2_ratio.bw
bigwigCompare -b1 ../oriFiles/actH2AZ_pos_R1.bw -b2 ../oriFiles/H2AZ_pos_R1.bw --operation ratio -bs 20 -p max -o POS_ac_h2_ratio.bw
############################################################################################################################

bedtools intersect -a oriFiles/cMyc_peaks.bed -b oriFiles/TIP60_peaks.bed > cmyc_tip60_peaks.bed
# 3650 and 30009 > 3039
annotatePeaks.pl cmyc_tip60_peaks.bed mm10 -annStats myc_h2az_peaks.annStats > cmyc_tip60_peaks.anno 

 grep "Intergenic" cmyc_tip60_peaks.anno|cut -f2,3,4,10,16 > cmyc_tip60_peaks_intergenic.bed
 grep "promoter-TSS" cmyc_tip60_peaks.anno|cut -f2,3,4,10,16 > cmyc_tip60_peaks_promoter.bed
