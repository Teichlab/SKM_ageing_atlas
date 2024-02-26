# 画细胞因子的基因差异表达 Figure 4
rm(list = ls())
library(tidyverse)
library(readxl)
library(data.table)
library (ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

DE_genes_all = fread("/public/LiuTL/Rstudio/skeletal_muscle/DE_analysis_LMM_v3/data/de_celltype_DE_human_scell_granular_v3_2023May.txt")
DE_genes <- DE_genes_all[which(DE_genes_all$celltype %in% c('B_naive','B_memory','B-plasma','CD16-NK','CD16+NK','CD4+CD8+T-naive','CD4+T','CD8+T','CD8+T-CRTAM+','CD14+Mono','CD16+Mono','Mϕ_CD16hi','cDC1','cDC2','Cyc-cDC2','Mast_cell',
                                                            'Adv_FB','Inter_FB','Par_FB','Tenocyte','Perineural_FB','Endoneural_FB-CDH19+','Endoneural_FB-TAC1+','mSchwann_cell','nmSchwann_cell',
                                                            'Artery','Arteriole','Cap','Cap-CCL2+','Cap-Ven','Vein','Vein-CCL2+','Lymphatic','SMC','SMC-RAMP1+','SMC-CCL2+','SMC-PC','Pericyte','Pericyte-CCL2+','Pericyte-CCL26+')),]

c('B_naive','B_memory','B-plasma','CD16-NK','CD16+NK','CD4+CD8+T-naive','CD4+T','CD8+T','CD8+T-CRTAM+','CD14+Mono','CD16+Mono','Mϕ_CD16hi','cDC1','cDC2','Cyc-cDC2','Mast_cell',
  'Adv_FB','Inter_FB','Par_FB','Tenocyte','Perineural_FB','Endoneural_FB-CDH19+','Endoneural_FB-TAC1+','mSchwann_cell','nmSchwann_cell',
  'Artery','Arteriole','Cap','Cap-CCL2+','Cap-Ven','Vein','Vein-CCL2+','Lymphatic','SMC','SMC-RAMP1+','SMC-CCL2+','SMC-PC','Pericyte','Pericyte-CCL2+','Pericyte-CCL26+')


colnames(DE_genes)[2] = "SYMBOL"
DE_genes[, log2fc:=beta_old-beta_young]
DE_genes[, cell_type_common:=celltype]
DE_genes[, REGULATION:=case_when(
  (log2fc>0 & ltsr > 0.9) ~ "UP",
  (log2fc<0 & ltsr > 0.9) ~ "DW",
  TRUE ~ 'none')]


gene_list = c('CCL2','CCL3','CCL4','CXCL8','IL6','IL10')


fraction1 <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/Veronika_ICM_data/write/ICM_Immune_cell_scvi_cytokines_fraction.csv',row.names = 1)
fraction2 <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/Veronika_ICM_data/write/ICM_FB_Schwann_scvi_cytokines_fraction.csv',row.names = 1)
fraction3 <- read.csv('/public/LiuTL/jupyter_notebook/skeletal_muscle/Veronika_ICM_data/write/ICM_Endo_SMC_scvi-v2_cytokines_fraction.csv',row.names = 1)

fraction <- rbind(fraction1,fraction2,fraction3)

unique(DE_genes$celltype)

ctypes = c('B_naive','B_memory','B-plasma','CD16-NK','CD16+NK','CD4+CD8+T-naive','CD4+T','CD8+T','CD8+T-CRTAM+','CD14+Mono','CD16+Mono','Mϕ_CD16hi','cDC1','cDC2','Cyc-cDC2','Mast_cell',
           'Adv_FB','Inter_FB','Par_FB','Tenocyte','Perineural_FB','Endoneural_FB-CDH19+','Endoneural_FB-TAC1+','mSchwann_cell','nmSchwann_cell',
           'Artery','Arteriole','Cap','Cap-CCL2+','Cap-Ven','Vein','Vein-CCL2+','Lymphatic','SMC','SMC-RAMP1+','SMC-CCL2+','SMC-PC','Pericyte','Pericyte-CCL2+','Pericyte-CCL26+')


gene_subset =DE_genes[(SYMBOL %in% gene_list) & (celltype %in% ctypes)]
gene_subset$celltype= factor(gene_subset$celltype, levels = ctypes)
gene_subset$SYMBOL = factor(gene_subset$SYMBOL, levels = rev(gene_list))
gene_subset$REGULATION = factor(gene_subset$REGULATION, levels = c("none", "UP", "DW"))

gene_subset$idenify <- paste(gene_subset$SYMBOL,gene_subset$celltype,sep = '_')
fraction$idenify <- paste(fraction$gene,fraction$cell_type,sep = '_')
gene_subset <- merge(gene_subset,fraction,by='idenify')

pdf( "./figure/Im_Fb_Ec_SMC_cytokines_v4.pdf", width = 9, height = 3)
ggplot(gene_subset, aes(x= SYMBOL, y=celltype, size=fraction, fill=log2fc, color = REGULATION, group = celltype)) + 
  geom_point(pch=21) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) + 
  # scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu"))) + 
  scale_fill_gradient2(midpoint = 0, low = "royalblue", mid = "white",
                       high = "red", space = "Lab" ) + 
  scale_color_manual(values = c("gray70",  "red","blue")) + 
  scale_radius(limits = range(0.0000000000000000000000001, 1)) + coord_flip()
dev.off()
