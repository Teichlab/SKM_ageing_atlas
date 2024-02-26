# 画RNA聚合酶I和ribosome biogenesis的基因差异表达 Figure 2D
rm(list = ls())
library(tidyverse)
library(readxl)
library(data.table)
library (ggplot2)
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


# gene_list = c('DLK1','PAX7','POLR1A','POLR1B','POLR1C','POLR1D','POLR1E','WDR12','EBNA1BP2','MRTO4','WDR74','RRS1','RRP1','RSL1D1','RPL7L1','BRIX1','NIFK','IMP4','UTP3','NOP58','UTP23','WDR3','RRP9','DCAF13','NOP56','NOB1','MPHOSPH10','GAR1','DKC1','NOLC1','DDX21','XRN2','NIP7','SURF6','PAK1IP1','NOP16','URB2','SRSF6','MRPL36','SNRPD1','SNRPB','ACTB','MRPS26','MRPL14','MRPL13','PUS7','METTL1','FARSB')
gene_list = c('POLR1A','POLR1B','POLR1C','POLR1D','POLR1E','WDR74','WDR3','MRTO4','EBNA1BP2','RRS1','RRP1','RRP9','BRIX1','IMP4','UTP3','NOP56','NOP58','DCAF13','NOB1','MPHOSPH10','GAR1','DDX21','NIP7','SURF6','SNRPB')
# 'POLR1A','POLR1B','POLR1C','POLR1D',  'RPA3','POLR1E',
# gene_list = c('BOP1','BRIX1','DDX28','DHX30','FASTKD2','MDN1','MRTO4','NLE1','NOP2','PPAN','RPF2','RRS1','TRAF7','ABT1','ERAL1','PRKDC','PWP2','RRP7A','XRCC5','C1QBP','DDX3X','EIF6','SBDS')

# ribosome assembly gene set
# gene_list = c('BOP1','BRIX1','DDX28','DHX30','FASTKD2','MDN1','MRPL20','MRTO4','NLE1','NOP2','PPAN','RPF2','RRS1','TRAF7','ABT1','ERAL1','MRPS11','MRPS7','PRKDC','PWP2','RRP7A','XRCC5','C1QBP','DDX3X','EIF6','SBDS')

# ribosome protein gene set
gene_list = c('RPLP0','RPL3','RPL4','RPL5','RPL6','RPL7','RPL7A','RPL8','RPL9','RPL10','RPL11','RPL12','RPL13','RPL13A','RPL14','RPL15','RPL17','RPL18','RPL18A','RPL19','RPL21','RPL22','RPL22L1',
              'RPL23','RPL23A','RPL24','RPL26','RPL27','RPL27A','RPL28','RPL29','RPL30','RPL31','RPL32','RPL34',
              'RPL35','RPL35A','RPL36','RPL36A','RPL37','RPL37A','RPL38','RPL39','RPL41','RPSA','RPS2','RPS3','RPS3A','RPS4Y1','RPS4Y2','RPS4X','RPS5','RPS6','RPS7','RPS8','RPS9','RPS10','RPS11','RPS12','RPS13','RPS14',
              'RPS15','RPS16','RPS17','RPS18','RPS19','RPS20','RPS21','RPS23','RPS24','RPS25','RPS26','RPS27','RPS27A','RPS28','RPS29')


fraction <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/Veronika_ICM_data/write/ICM_MuSC_scvi_v2_ribosome_biogenesis_fraction.csv',row.names = 1)

# ribosome assembly
# fraction <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/human_SKM_subcluster_scvi/write/humsn_SKM_MuSC_scVI_v2.0-analysis3_ribosome_assembly_fraction.csv',row.names = 1)

# ribosome protein
fraction <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/human_SKM_subcluster_scvi/write/humsn_SKM_MuSC_scVI_v2.0-analysis3_ribosome_protein_fraction.csv',row.names = 1)



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
  # scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) + 
  scale_fill_gradient2(midpoint = 0, low = "royalblue", mid = "white",
                       high = "red", space = "Lab" ) + 
  scale_color_manual(values = c("gray70","blue",  "red")) + 
  scale_radius(limits = range(0, 1)) + coord_flip()
dev.off()
