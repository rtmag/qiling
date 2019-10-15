##################################################################################################################
# INSTALL
# Crane 2015 code 11/oct/2019
git clone https://github.com/dekkerlab/crane-nature-2015

cd crane-nature-2015/

#update cpam build lib
sudo cpan App::cpanminus
cpanm Module::Build

# build perl script
perl Build.pl
./Build
./Build install

##################################################################################################################
# MODIFY juicer .hic dump in R
# R
# juicer_to_hicPro_denseMatrix.Rscript
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)

denseMatrix <- as.character(args[1])
resolution <- as.numeric(args[2])
chr <- as.character(args[3])
sampleID <- as.character(args[4])
matrix <- read.table(denseMatrix,header=FALSE)
mat_size = dim(matrix)[1]
ix1 <- seq(0, (mat_size-1)*resolution, by=resolution)
ix2 <- seq(resolution, ((mat_size-1)*resolution)+resolution, by=resolution)
labels<-paste0("bin",1:dim(matrix)[1],"|mm10|",chr,":",ix1,"-",ix2)
rownames(matrix) <- labels
colnames(matrix) <- labels
outfile <- paste0(sampleID,"_",chr,"_observed_NONE_25kb_D_hicPROformat.txt")
write.table(matrix,outfile,quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)

##################################################################################################################
# Runing perl script
perl /root/qiling/crane-nature-2015/scripts/matrix2insulation.pl -i chr19_observed_NONE_25kb_D_hicPROformat.txt
perl /root/qiling/crane-nature-2015/scripts/matrix2insulation.pl -i /root/qiling/hicProMatrix/KO_19_observed_NONE_25kb_D_hicPROformat.txt

##################################################################################################################


for matrixFile in /root/qiling/dense_matrix/*/*;
do ls -lh $matrixFile; 
sampleOut=$(basename $matrixFile) ;
chrN=${sampleOut//_observed_NONE_25kb_D.txt/} ;
chrN=${chrN//WT_/} ;
chrN=${chrN//KO_/} ;
res=25000 ;
sampleID="$(echo $sampleOut|perl -pe 's/\_.+//g')" ;
Rscript juicer_to_hicPro_denseMatrix.R $matrixFile $res $chrN $sampleID;
done

for matrixFile in /root/qiling/hicProMatrix/*_hicPROformat.txt;
do ls -lh $matrixFile; 
perl /root/qiling/crane-nature-2015/scripts/matrix2insulation.pl -i $matrixFile ;
done

more mm10.chrom.sizes|perl -pe 's/chr//g' > mm10.chrom.sizes_noChr 

cat KO_*.insulation.bedGraph|grep -v "track type=bedGraph name"|grep -v -w "NA" > ../KO_insulationScore_observed_NONE_25kb.bedGraph
cat WT_*.insulation.bedGraph|grep -v "track type=bedGraph name"|grep -v -w "NA" > ../WT_insulationScore_observed_NONE_25kb.bedGraph

more KO_insulationScore_observed_NONE_25kb.bedGraph|perl -pe "s/^/chr/g" > KO_insulationScore_observed_NONE_25kb_chrNames.bedGraph
more WT_insulationScore_observed_NONE_25kb.bedGraph|perl -pe "s/^/chr/g" > WT_insulationScore_observed_NONE_25kb_chrNames.bedGraph

/root/myPrograms/kentUtils/bin/bedGraphToBigWig KO_insulationScore_observed_NONE_25kb.bedGraph mm10.chrom.sizes_noChr KO_insulationScore_observed_NONE_25kb.bw &         #
/root/myPrograms/kentUtils/bin/bedGraphToBigWig WT_insulationScore_observed_NONE_25kb.bedGraph mm10.chrom.sizes_noChr WT_insulationScore_observed_NONE_25kb.bw &         #

/root/myPrograms/kentUtils/bin/bedGraphToBigWig KO_insulationScore_observed_NONE_25kb_chrNames.bedGraph mm10.chrom.sizes KO_insulationScore_observed_NONE_25kb_chrNames.bw &          #
/root/myPrograms/kentUtils/bin/bedGraphToBigWig WT_insulationScore_observed_NONE_25kb_chrNames.bedGraph mm10.chrom.sizes WT_insulationScore_observed_NONE_25kb_chrNames.bw &          #


