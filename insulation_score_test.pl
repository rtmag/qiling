##################################################################################################################
# INSTALL
# Crane 2015 code 11/oct/2019
git clone https://github.com/dekkerlab/crane-nature-2015

cd crane-nature-2015/
perl Build.PL
./Build
./Build install

##################################################################################################################
# MODIFY juicer .hic dump in R
# R
options(scipen=999)
matrix <- read.table("chr19_observed_NONE_25kb_D.txt",header=FALSE)
chr <- "chr19"
mat_size = dim(matrix)[1]
ix1 <- seq(0, (mat_size-1)*25000, by=25000)
ix2 <- seq(25000, ((mat_size-1)*25000)+25000, by=25000)
labels<-paste0("bin",1:dim(matrix)[1],"|mm10|",chr,":",ix1,"-",ix2)
rownames(matrix) <- labels
colnames(matrix) <- labels
write.table(matrix,"chr19_observed_NONE_25kb_D_hicPROformat.txt",quote=FALSE)

