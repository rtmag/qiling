library(ComplexHeatmap)
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))
options(bitmapType="cairo")

genes <- read.table("20200211_myc_targets_gene_list.txt")
genes<-as.character(genes[,1])

exp <- read.csv("Tip60_RNAseq_triplicates.csv")
exp$Gene <- toupper(exp$Gene)
exp <-exp[exp$Gene %in% genes,]
exp.sig <- exp[,8:13]
rownames(exp.sig) <- as.character(exp[,1])

exp.log2 <- as.matrix(log2(exp.sig))
exp.log2 = exp.log2[apply(exp.log2,1,sd)!=0,]

Heatmap(exp.log2,
show_row_names = TRUE,show_column_names = TRUE,name = "Expression",row_dend_reorder = TRUE, column_dend_reorder = TRUE,
clustering_distance_columns = "pearson", clustering_distance_rows = "pearson")
