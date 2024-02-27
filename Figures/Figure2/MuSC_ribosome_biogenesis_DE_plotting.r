# MuSC ribosome biogenesis DE plotting (Figure 2e)
rm(list = ls())
library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

DE_genes_all = fread("/public/LiuTL/Rstudio/skeletal_muscle/DE_analysis_LMM_v3/de_celltype_DE_human_scell_granular_v2_2023May.txt")
DE_genes <- DE_genes_all[which(DE_genes_all$celltype %in% c('MYOG+MuSCs','TNFRSF12A+MuSCs','ICAM1+MuSCs','MuSCs')),]

colnames(DE_genes)[2] = "SYMBOL"
DE_genes[, log2fc:=beta_old-beta_young]
DE_genes[, cell_type_common:=celltype]
DE_genes[, REGULATION:=case_when(
  (log2fc>0 & ltsr > 0.9) ~ "UP",
  (log2fc<0 & ltsr > 0.9) ~ "DW",
  TRUE ~ 'none')]

gene_list = c('POLR1A','POLR1B','POLR1C','POLR1D','POLR1E','WDR74','WDR3','MRTO4','EBNA1BP2','RRS1','RRP1','RRP9','BRIX1','IMP4','UTP3','NOP56','NOP58','DCAF13','NOB1','MPHOSPH10','GAR1','DDX21','NIP7','SURF6','SNRPB')
fraction <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/Veronika_ICM_data/write/ICM_MuSC_scvi_v2_ribosome_biogenesis_fraction.csv',row.names = 1)

unique(DE_genes$celltype)

ctypes = c("MYOG+MuSCs","TNFRSF12A+MuSCs","ICAM1+MuSCs","MuSCs")

gene_subset =DE_genes[(SYMBOL %in% gene_list) & (celltype %in% ctypes)]
gene_subset$celltype= factor(gene_subset$celltype, levels = ctypes)
gene_subset$SYMBOL = factor(gene_subset$SYMBOL, levels = rev(gene_list))
gene_subset$REGULATION = factor(gene_subset$REGULATION, levels = c("none", "UP", "DW"))

gene_subset$idenify <- paste(gene_subset$SYMBOL,gene_subset$celltype,sep = '_')
fraction$idenify <- paste(fraction$gene,fraction$cell_type,sep = '_')
gene_subset <- merge(gene_subset,fraction,by='idenify')

pdf( "./figure/MuSC_gene_ribosome_biogenesis_v2.pdf", width = 4, height = 8)
ggplot(gene_subset, aes(x= SYMBOL, y=celltype, size=fraction, fill=log2fc, color = REGULATION, group = celltype)) + 
  geom_point(pch=21) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) + 
  scale_fill_gradient2(midpoint = 0, low = "royalblue", mid = "white",
                       high = "red", space = "Lab" ) + 
  scale_color_manual(values = c("gray70","blue",  "red")) + 
  scale_radius(limits = range(0, 1)) + coord_flip()
dev.off()
