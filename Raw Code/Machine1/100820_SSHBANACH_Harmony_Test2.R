  
# Load required Libraries

library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)
library(data.table)

#Set paths

mA5.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Lung-Bleo/"
BDL2.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Liver-BDL/"
CCl4.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Liver-CCl4/"
# Bleo2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
# BDL2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
# CCl4.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"

banach.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/"
mA5Output.path <- paste0(banach.path,"Kidney-UUO/Outputs/")
Bleo2Output.path <-  paste0(banach.path,"Lung-Bleo/Outputs/")
BDL2Output.path <-  paste0(banach.path,"Liver-BDL/Outputs/")
CCl4Output.path  <-  paste0(banach.path,"Liver-CCl4/Outputs/")

#Set date

Today.date <- "070820" 

# Set colour palettes

mA5.colpal <- c("#CA9858",	"#FF7F00",	"#C8828D",	"#FFCA62")

Bleo2.colpal <- c("#9683B0",	"#6A3D9A",	"#C8828D",	"#DA98D6")

BDL2.colpal <- c("#8EB26E",	"#88BCB2",	"#33A02C",	"#00CDB3",	"#C8828D",	"#F4A1AF",	"#B2DF8A",	"#AFF1E1")

CCl4.colpal <- c("#74909E",	"#167EDB",	"#C8828D",	"#A6CEE3")

 
#Load datasets
mA5 <- readRDS(paste0(mA5.path, Today.date, "_mA5_UMAP.RDS"))
Bleo2 <- readRDS(paste0(Bleo2.path, Today.date, "_Bleo2_UMAP.RDS"))
BDL2 <- readRDS(paste0(BDL2.path, Today.date, "_BDL2_UMAP.RDS"))
CCl4 <- readRDS(paste0(CCl4.path, Today.date, "_CCl4_UMAP.RDS"))

# mA5.markers <- FindAllMarkers(mA5, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# # mA5.markersOUTPUT <- mA5.markers %>% group_by(cluster) %>% arrange(cluster, pct.2)
# # write.csv(mA5.markersOUTPUT, file = paste0(mA5Output.path, Today.date, "_mA5_biomarkers.csv"))
# 
# Bleo2.markers <- FindAllMarkers(Bleo2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# # Bleo2.markersOUTPUT <- Bleo2.markers %>% group_by(cluster) %>%  arrange(cluster, pct.2)
# # write.csv(Bleo2.markersOUTPUT, file = paste0(Bleo2Output.path, Today.date, "_Bleo2_biomarkers.csv"))
# 
# BDL2.markers <- FindAllMarkers(BDL2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# # BDL2.markersOUTPUT <- BDL2.markers %>% group_by(cluster) %>%  arrange(cluster, pct.2)
# # write.csv(BDL2.markersOUTPUT, file = paste0(BDL2Output.path, Today.date, "_BDL2_biomarkers.csv"))
# 
# CCl4.markers <- FindAllMarkers(CCl4, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# # CCl4.markersOUTPUT <- CCl4.markers %>% group_by(cluster) %>%  arrange(cluster, pct.2)
# write.csv(CCl4.markersOUTPUT, file = paste0(CCl4Output.path, Today.date, "_CCl4_biomarkers.csv"))


# # Cluster and Lineage Annotation + Save
# ### cell type annotation ###
# mA5_cluster.annotation <- c("Blood B cell 1", "Kidney macrophage 1", "Kidney macrophage 2", "UUO blood Ly6Chi monocyte", "Blood Ly6Clo monocyte", "CD8+ T cell 1", "Kidney macrophage 3", "Kidney macrophage 4", "Cd4+ T cell 1", "NK cell 1", "Kidney monocyte", "Healthy blood Ly6Chi monoyte", "Kidney cDC2 1", "Cd4+ T cell 2", "Kidney cDC2 2", "Kidney B cell", "Kidney cDC1", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney endothelial cells", "Blood neutrophils", "NK cell 2", "pDC","Neutrophil/B cell doublet","Basophil","Proliferating","Neutrophil/T cell doublet", "ILC2", "CCR7+ cDC", "Cd8+ T cell 2", "Monocyte/T cell doublet", "Platelet")
# 
# names(mA5_cluster.annotation) <- levels(mA5)
# mA5 <- RenameIdents(mA5, mA5_cluster.annotation)
# 
# mA5_cell.data <- data.table(barcode = colnames(mA5),
#                             Cell.type = Idents(mA5))
# 
# mA5_cell.data <- data.frame(mA5_cell.data, row.names = mA5_cell.data$barcode)
# mA5_cell.data$barcode <- NULL
# mA5 <- AddMetaData(mA5, mA5_cell.data, col.name = "Cell.type")
# 
# ### save annotated UMAP ###
# annotated.umap.plot.split <- DimPlot(mA5, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_aUMAP_cluster_OGIDsplit.png"), plot = annotated.umap.plot.split, width = 45, height = 10, units = "cm")
# 
# annotated.umap.plot.merge <- DimPlot(mA5, reduction = "umap", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_aUMAP_cluster.png"), plot = annotated.umap.plot.merge, width = 30, units = "cm")
# 
# 
# ### cell lineage annotation ###
# mA5_cell.data <- data.table(barcode = colnames(mA5),
#                             Cell.type = Idents(mA5))
# 
# mA5_lineage.annotation <- c("B cell",	"MP",	"MP",	"MP",	"MP",	"T cell/ILC",	"MP",	"MP",	"T cell/ILC",	"T cell/ILC",	
#                             "MP",	"MP",	"MP",	"T cell/ILC",	"MP",	"B cell",	"MP",	"MP",	"MP",	"MP",	"Endothelia",	
#                             "Neutrophil",	"T cell/ILC",	"MP",	"Doublet",	"Basophil",	"Proliferating",	"Doublet",	
#                             "T cell/ILC",	"MP",	"T cell/ILC",	"Doublet",	"Platelet")
# 
# 
# mA5_lineage.data <- data.table(Cell.type = mA5_cluster.annotation, lineage = mA5_lineage.annotation)
# meta.data <- merge(mA5_cell.data, mA5_lineage.data, by = "Cell.type")
# meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
# meta.data$barcode <- NULL
# meta.data$Cell.type <- NULL
# mA5 <- AddMetaData(mA5, meta.data, col.name = "lineage")
# # save annotated UMAP
# lineage.umap.plot <- DimPlot(mA5, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 45, height = 10, units = "cm")
# 
# lineage.umap.plot <- DimPlot(mA5, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")
# 
# saveRDS(mA5, file = paste0(mA5.path, Today.date, "_mA5_annotated.RDS"))
# 
# #- Bleo2 -
# 
# ### cell type annotation ###
# Bleo2_cluster.annotation <- c("Blood B cell 1","Lung alveolar macrophage 1","Blood Ly6Chi Monocyte","Blood Ly6Clo Monocyte 1","Lung B Cell 1","Cd4+ T cell 1","NK cell 1", "Blood CD8+ T cell","Blood CD4+ T cell 2","Lung CD4+ T cell 2","Lung CD8+ T cell","Lcn2+ Monocyte","cDC1","cDC2","Lung interstitial macrophage 1","Lung Monocyte","Lung macrophage 3","Lung interstitial macrophage 2","Blood B cell 2","Blood Ly6Clo Monocyte 2","Proliferating T cell","Proliferating cDC","ILC2","CD8+ T cell 2","Proliferating Lung Alveolar macrophage","Blood IFN-primed Ly6Chi Monocyte","T cell/B cell doublet","Basophil","pDC","Lung macrophage 4","Lung endothelial cells","CCR7+ cDC","Lung alveolar macrophage 2")
# 
# names(Bleo2_cluster.annotation) <- levels(Bleo2)
# Bleo2 <- RenameIdents(Bleo2, Bleo2_cluster.annotation)
# 
# Bleo2_cell.data <- data.table(barcode = colnames(Bleo2),
#                               Cell.type = Idents(Bleo2))
# 
# Bleo2_cell.data <- data.frame(Bleo2_cell.data, row.names = Bleo2_cell.data$barcode)
# Bleo2_cell.data$barcode <- NULL
# Bleo2 <- AddMetaData(Bleo2, Bleo2_cell.data, col.name = "Cell.type")
# 
# ### save annotated UMAP ###
# annotated.umap.plot.split <- DimPlot(Bleo2, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2_aUMAP_cluster_OGIDsplit.png"), plot = annotated.umap.plot.split, width = 45, height = 10, units = "cm")
# 
# annotated.umap.plot.merge <- DimPlot(Bleo2, reduction = "umap", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2_aUMAP_cluster.png"), plot = annotated.umap.plot.merge, width = 30, units = "cm")
# 
# 
# ### cell lineage annotation ###
# Bleo2_cell.data <- data.table(barcode = colnames(Bleo2),
#                               Cell.type = Idents(Bleo2))
# 
# Bleo2_lineage.annotation <- c("B cell","MP","MP","MP","B cell","T cell/ILC","T cell/ILC","T cell/ILC","T cell/ILC","T cell/ILC", "T cell/ILC","MP","MP","MP","MP","MP","MP","MP","B cell","MP","T cell/ILC","MP","T cell/ILC","T cell/ILC","MP","MP","Doublet","Basophil","MP","MP","Endothelia","MP","MP")
# 
# 
# Bleo2_lineage.data <- data.table(Cell.type = Bleo2_cluster.annotation, lineage = Bleo2_lineage.annotation)
# Bleo2_meta.data <- merge(Bleo2_cell.data, Bleo2_lineage.data, by = "Cell.type")
# Bleo2_meta.data <- data.frame(Bleo2_meta.data, row.names = Bleo2_meta.data$barcode)
# Bleo2_meta.data$barcode <- NULL
# Bleo2_meta.data$Cell.type <- NULL
# Bleo2 <- AddMetaData(Bleo2, Bleo2_meta.data, col.name = "lineage")
# # save annotated UMAP
# lineage.umap.plot <- DimPlot(Bleo2, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename =  paste0(Bleo2Output.path, Today.date, "_Bleo2_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 45, height = 10, units = "cm")
# 
# lineage.umap.plot <- DimPlot(Bleo2, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename =  paste0(Bleo2Output.path, Today.date, "_Bleo2_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")
# 
# 
# saveRDS(Bleo2, file = paste0(Bleo2.path, Today.date, "_Bleo2_annotated.RDS"))
# 
# #- CCl4 -
# 
# ### cell type annotation ###
# CCl4_cluster.annotation <- c("Red blood cell/myeloid doublet","Blood Ly6Clo monocyte", "CCl4 Blood Ly6Clo monocyte","Healthy Blood Ly6Chi monocyte", "B cell 1", "CD8+ T cell 1", "CCl4 Blood Ly6Chi monocyte 1", "Liver macrophage 1", "NK cell 1","Cd4+ T cell 1","pDC", "cDC2","Liver Kupffer Cell 1","cDC1","CCl4 Liver Kupffer Cell 2","Blood IFN-primed Ly6Chi monocyte","Liver monocyte","CCl4 Blood Ly6Chi monocyte 2","Neutrophil","CCl4 Blood Ly6Chi monocyte 3","NK cell 2","Proliferating","Liver Endothelial Cells","B cell 2","B cell 3","Liver CCR7+ cDC","Macrophage/endothelial Doublet","Basophil","CD8+ T cell 2","Monocyte/T cell doublet","Eosinophil")
# 
# names(CCl4_cluster.annotation) <- levels(CCl4)
# CCl4 <- RenameIdents(CCl4, CCl4_cluster.annotation)
# 
# CCl4_cell.data <- data.table(barcode = colnames(CCl4),
#                              Cell.type = Idents(CCl4))
# 
# CCl4_cell.data <- data.frame(CCl4_cell.data, row.names = CCl4_cell.data$barcode)
# CCl4_cell.data$barcode <- NULL
# CCl4 <- AddMetaData(CCl4, CCl4_cell.data, col.name = "Cell.type")
# 
# ### save annotated UMAP ###
# annotated.umap.plot.split <- DimPlot(CCl4, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4_aUMAP_cluster_OGIDsplit.png"), plot = annotated.umap.plot.split, width = 45, height = 10, units = "cm")
# 
# annotated.umap.plot.merge <- DimPlot(CCl4, reduction = "umap", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4_aUMAP_cluster.png"), plot = annotated.umap.plot.merge, width = 30, units = "cm")
# 
# 
# ### cell lineage annotation ###
# CCl4_cell.data <- data.table(barcode = colnames(CCl4),
#                              Cell.type = Idents(CCl4))
# 
# CCl4_lineage.annotation <- c("Doublet","MP","MP","MP","B cell","T cell/ILC","MP","MP","T cell/ILC","T cell/ILC","MP","MP","MP","MP","MP","MP","MP","MP","Neutrophil","MP","T cell/ILC","Proliferating","Endothelia","B cell","B cell","MP","Doublet","Basophil","T cell/ILC","Doublet","Eosinophil")
# 
# 
# CCl4_lineage.data <- data.table(Cell.type = CCl4_cluster.annotation, lineage = CCl4_lineage.annotation)
# meta.data <- merge(CCl4_cell.data, CCl4_lineage.data, by = "Cell.type")
# meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
# meta.data$barcode <- NULL
# meta.data$Cell.type <- NULL
# CCl4 <- AddMetaData(CCl4, meta.data, col.name = "lineage")
# # save annotated UMAP
# lineage.umap.plot <- DimPlot(CCl4, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename =  paste0(CCl4Output.path, Today.date, "_CCl4_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 45, height = 10, units = "cm")
# 
# lineage.umap.plot <- DimPlot(CCl4, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename =  paste0(CCl4Output.path, Today.date, "_CCl4_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")
# 
# saveRDS(CCl4, file = paste0(CCl4.path, Today.date, "_CCl4_annotated.RDS"))
# # # Remove Doublet populations + Save
# # 
# # 
# mA5_cleaned <- subset(mA5, subset = lineage != "Doublet")
# saveRDS(mA5_cleaned, file = paste0(mA5.path, Today.date, "_mA5_CLEANED.RDS"))
# 
# Bleo2_cleaned <- subset(Bleo2, subset = lineage != "Doublet")
# saveRDS(Bleo2_cleaned, file = paste0(Bleo2.path, Today.date, "_Bleo2_CLEANED.RDS"))
# 
# CCl4_cleaned <- subset(CCl4, subset = lineage != "Doublet")
# saveRDS(CCl4_cleaned, file = paste0(CCl4.path, Today.date, "_CCl4_CLEANED.RDS"))
# 
# 
# 
#   

#Load stuff
mA5 <- readRDS(paste0(mA5.path, Today.date, "_mA5_annotated.RDS"))
Bleo2 <- readRDS(paste0(Bleo2.path, Today.date, "_Bleo2_annotated.RDS"))
CCl4 <- readRDS(paste0(CCl4.path, Today.date, "_CCl4_annotated.RDS"))

mA5_cleaned <- readRDS(paste0(mA5.path, Today.date, "_mA5_CLEANED.RDS"))

Bleo2_cleaned <- readRDS(paste0(Bleo2.path, Today.date, "_Bleo2_CLEANED.RDS"))

CCl4_cleaned  <- readRDS(paste0(CCl4.path, Today.date, "_CCl4_CLEANED.RDS"))

# # Post-doublet visualisationss
# 
# mA5cleaned_UMAP <- DimPlot(mA5_cleaned, reduction = "umap",label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_aUMAP_cluster.png"), plot = mA5cleaned_UMAP, width = 30, units = "cm")
# 
# mA5cleaned_UMAP_split <- DimPlot(mA5_cleaned, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_aUMAP_cluster_OGIDsplit.png"), plot = mA5cleaned_UMAP_split, width = 45, height = 10, units = "cm")
# 
# mA5cleaned_lineage_split <- DimPlot(mA5_cleaned, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_aUMAP_lineage_OGIDsplit.png"), plot = mA5cleaned_lineage_split, width = 45, height = 10, units = "cm")
# 
# mA5cleaned_lineage <- DimPlot(mA5_cleaned, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.02) + NoLegend()
# ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_aUMAP_lineage.png"), plot = mA5cleaned_lineage, width = 30, units = "cm")
# 
# 
# #-
# 
# Bleo2cleaned_UMAP <- DimPlot(Bleo2_cleaned, reduction = "umap",label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_aUMAP_cluster.png"), plot = Bleo2cleaned_UMAP, width = 30, units = "cm")
# 
# Bleo2cleaned_UMAP_split <- DimPlot(Bleo2_cleaned, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_aUMAP_cluster_OGIDsplit.png"), plot = Bleo2cleaned_UMAP_split, width = 45, height = 10, units = "cm")
# 
# Bleo2cleaned_lineage_split <- DimPlot(Bleo2_cleaned, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_aUMAP_lineage_OGIDsplit.png"), plot = Bleo2cleaned_lineage_split, width = 45, height = 10, units = "cm")
# 
# Bleo2cleaned_lineage <- DimPlot(Bleo2_cleaned, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.02) + NoLegend()
# ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_aUMAP_lineage.png"), plot = Bleo2cleaned_lineage, width = 30, units = "cm")
# 
# #- 
# 
# CCl4cleaned_UMAP <- DimPlot(CCl4_cleaned, reduction = "umap",label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_aUMAP_cluster.png"), plot = CCl4cleaned_UMAP, width = 30, units = "cm")
# 
# CCl4cleaned_UMAP_split <- DimPlot(CCl4_cleaned, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_aUMAP_cluster_OGIDsplit.png"), plot = CCl4cleaned_UMAP_split, width = 45, height = 10, units = "cm")
# 
# CCl4cleaned_lineage_split <- DimPlot(CCl4_cleaned, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_aUMAP_lineage_OGIDsplit.png"), plot = CCl4cleaned_lineage_split, width = 45, height = 10, units = "cm")
# 
# CCl4cleaned_lineage <- DimPlot(CCl4_cleaned, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.02) + NoLegend()
# ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_aUMAP_lineage.png"), plot = CCl4cleaned_lineage, width = 30, units = "cm")

mA5.markers <- FindAllMarkers(mA5, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# mA5.markersOUTPUT <- mA5.markers %>% group_by(cluster) %>% arrange(cluster, pct.2)
# write.csv(mA5.markersOUTPUT, file = paste0(mA5Output.path, Today.date, "_mA5_biomarkers.csv"))

Bleo2.markers <- FindAllMarkers(Bleo2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# Bleo2.markersOUTPUT <- Bleo2.markers %>% group_by(cluster) %>%  arrange(cluster, pct.2)
# write.csv(Bleo2.markersOUTPUT, file = paste0(Bleo2Output.path, Today.date, "_Bleo2_biomarkers.csv"))

BDL2.markers <- FindAllMarkers(BDL2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# BDL2.markersOUTPUT <- BDL2.markers %>% group_by(cluster) %>%  arrange(cluster, pct.2)
# write.csv(BDL2.markersOUTPUT, file = paste0(BDL2Output.path, Today.date, "_BDL2_biomarkers.csv"))

CCl4.markers <- FindAllMarkers(CCl4, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

# Individual Visualisations GENES -> Heatmap of cluster/marker genes
mA5_top10 <- mA5.markers  %>% group_by(cluster) %>% top_n(n = 10, wt = pct.2)
mA5_clusterHM <- DoHeatmap(mA5, features = mA5_top10$genes) + NoLegend()
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_clusterHM.png"), plot = mA5_clusterHM, width = 45, units = "cm")

mA5CLEANED_top10 <- mA5.markers% >% filter(cluster != 24 & cluster != 27 & cluster != 31) %>% group_by(Cell.type) %>% top_n(n = 10, wt = pct.2)
mA5_cellHM <- DoHeatmap(mA5_cleaned, group.by = "Cell.type", features = mA5CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_celltypeHM.png"), plot = mA5_celltypeHM, width = 45, units = "cm")

mA5_lineageHM <- DoHeatmap(mA5_cleaned, groub.by = "lineage", features = mA5CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5_lineage_HM.png"), plot = mA5_celltypeHM, width = 45, units = "cm")

#-
Bleo2_top10 <- Bleo2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = pct.2)
Bleo2_clusterHM <- DoHeatmap(Bleo2, features = Bleo2_top10$genes) + NoLegend()
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2_clusterHM.png"), plot = Bleo2_clusterHM, width = 45, units = "cm")

Bleo2CLEANED_top10 <- Bleo2.markers %>% filter(cluster != 26)  %>% group_by(Cell.type) %>% top_n(n = 10, wt = pct.2)

Bleo2_cellHM <- DoHeatmap(Bleo2_cleaned, group.by = "Cell.type", features = Bleo2CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2_celltypeHM.png"), plot = Bleo2_celltypeHM, width = 45, units = "cm")

Bleo2_lineageHM <- DoHeatmap(Bleo2_cleaned, group.by = "lineage", features = Bleo2CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2_lineage_HM.png"), plot = Bleo2_celltypeHM, width = 45, units = "cm")

#-

BDL2_top10 <- BDL2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = pct.2)
BDL2_clusterHM <- DoHeatmap(BDL2, features = BDL2_top10$genes) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_clusterHM.png"), plot = BDL2_clusterHM, width = 45, units = "cm")

#-

CCl4_top10 <- CCl4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = pct.2)
CCl4_clusterHM <- DoHeatmap(CCl4, features = CCl4_top10$genes) + NoLegend()
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4_clusterHM.png"), plot = CCl4_clusterHM, width = 45, units = "cm")

CCl4CLEANED_top10 <- CCl4.markers %>% filter(cluster != 0 & cluster != 26 & cluster != 29)  %>% group_by(Cell.type) %>% top_n(n = 10, wt = pct.2)
CCl4_cellHM <- DoHeatmap(CCl4_cleaned, group.by = "Cell.type", features = CCl4CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4_celltypeHM.png"), plot = CCl4_celltypeHM, width = 45, units = "cm")

CCl4_lineageHM <- DoHeatmap(CCl4_cleaned, group.by = "lineage", features = CCl4CLEANED_top10$genes) + NoLegend()
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4_lineage_HM.png"), plot = CCl4_celltypeHM, width = 45, units = "cm")


# Individual Visualisations CELL -> Proportion
### Dataset membership vary by cell type ###
mA5$Cell.type <- factor(mA5$Cell.type, levels = c("Basophil",	"Blood B cell 1",	"Kidney B cell",	"Blood Ly6Clo monocyte",	"Healthy blood Ly6Chi monoyte",	"UUO blood Ly6Chi monocyte",	"NK cell 1",	"NK cell 2",	"Blood neutrophils",	"Cd4+ T cell 1",	"Cd4+ T cell 2",	"CD8+ T cell 1",	"Cd8+ T cell 2",	"ILC2",	"pDC",	"CCR7+ cDC",	"Kidney cDC1",	"Kidney cDC2 1",	"Kidney cDC2 2",	"Kidney endothelial cells",	"Kidney IFN-primed macrophage",	"Kidney macrophage 1",	"Kidney macrophage 2",	"Kidney macrophage 3",	"Kidney macrophage 4",	"Kidney monocyte",	"Kidney proliferating macrophage 1",	"Kidney proliferating macrophage 2",	"Monocyte/T cell doublet",	"Neutrophil/B cell doublet",	"Neutrophil/T cell doublet",	"Platelet",	"Proliferating"))


Idents(mA5) <- "Names"
table(Idents(mA5))
mA5_NameCellProp <- prop.table(table(Idents(mA5), mA5$Cell.type))
mA5_NameCellProp
mA5_CellProp <- as.data.frame(mA5_NameCellProp)

names(mA5_CellProp)[names(mA5_CellProp) == "Var1"] <- "Dataset"
names(mA5_CellProp)[names(mA5_CellProp) == "Var2"] <- "Cell_type"
names(mA5_CellProp)[names(mA5_CellProp) == "Freq"] <- "Proportion"

mA5_CellProp <- mA5_CellProp %>% filter(Cell_type != "Monocyte/T cell doublet" & Cell_type != "Neutrophil/T cell doublet" &Cell_type != "Neutrophil/B cell doublet")
mA5_CellProp

plot_mA5_CellPropBAR <- ggplot(mA5_CellProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = mA5.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_CellPropBAR_Celltype.png"), plot = plot_mA5_CellPropBAR, width = 30, units = "cm")

plot_mA5_CellPropBAR_Dataset <- ggplot(mA5_CellProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_CellPropBAR_Dataset.png"), plot = plot_mA5_CellPropBAR_Dataset, width = 30, units = "cm")

#-
Bleo2$Cell.type <- factor(Bleo2$Cell.type, levels = c("Basophil",	"Blood B cell 1",	"Blood B cell 2",	"Lung B Cell 1",	"Blood CD4+ T cell 2",	"Blood CD8+ T cell",	"Cd4+ T cell 1",	"CD8+ T cell 2",	"Lung CD4+ T cell 2",	"Lung CD8+ T cell",	"Proliferating T cell",	"ILC2",	"NK cell 1",	"Blood IFN-primed Ly6Chi Monocyte",	"Blood Ly6Chi Monocyte",	"Blood Ly6Clo Monocyte 1",	"Blood Ly6Clo Monocyte 2",	"Lcn2+ Monocyte",	"CCR7+ cDC",	"cDC1",	"cDC2",	"pDC",	"Proliferating cDC",	"Lung alveolar macrophage 1",	"Lung alveolar macrophage 1",	"Proliferating Lung Alveolar macrophage",	"Lung interstitial macrophage 1",	"Lung interstitial macrophage 2",	"Lung macrophage 3",	"Lung macrophage 4",	"Lung Monocyte",	"Lung endothelial cells",	"T cell/B cell doublet"))

Idents(Bleo2) <- "Names"
table(Idents(Bleo2))
Bleo2_NameCellProp <- prop.table(table(Idents(Bleo2), Bleo2$Cell.type))
Bleo2_NameCellProp
Bleo2_CellProp <- as.data.frame(Bleo2_NameCellProp)

names(Bleo2_CellProp)[names(Bleo2_CellProp) == "Var1"] <- "Dataset"
names(Bleo2_CellProp)[names(Bleo2_CellProp) == "Var2"] <- "Cell_type"
names(Bleo2_CellProp)[names(Bleo2_CellProp) == "Freq"] <- "Proportion"

Bleo2_CellProp <- Bleo2_CellProp %>% filter(Cell_type != "T cell/B cell doublet")
Bleo2_CellProp

plot_Bleo2_CellPropBAR <- ggplot(Bleo2_CellProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = Bleo2.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_CellPropBAR_Celltype.png"), plot = plot_Bleo2_CellPropBAR, width = 30, units = "cm")

plot_Bleo2_CellPropBAR_Dataset <- ggplot(Bleo2_CellProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_CellPropBAR_Dataset.png"), plot = plot_Bleo2_CellPropBAR_Dataset, width = 30, units = "cm")

#-
CCl4$Cell.type <- factor(CCl4$Cell.type, levels = c("B cell 1","B cell 2","B cell 3","Basophil","Blood IFN-primed Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","CCl4 Blood Ly6Clo monocyte","Healthy Blood Ly6Chi monocyte","Liver monocyte","Cd4+ T cell 1","CD8+ T cell 1","CD8+ T cell 2","Eosinophil","Neutrophil","NK cell 1","NK cell 2","Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Liver CCR7+ cDC","cDC1","cDC2","pDC","Liver Endothelial Cells","Proliferating","Macrophage/endothelial Doublet","Monocyte/T cell doublet","Red blood cell/myeloid doublet"))

Idents(CCl4) <- "Names"
table(Idents(CCl4))
CCl4_NameCellProp <- prop.table(table(Idents(CCl4), CCl4$Cell.type))
CCl4_NameCellProp
CCl4_CellProp <- as.data.frame(CCl4_NameCellProp)

names(CCl4_CellProp)[names(CCl4_CellProp) == "Var1"] <- "Dataset"
names(CCl4_CellProp)[names(CCl4_CellProp) == "Var2"] <- "Cell_type"
names(CCl4_CellProp)[names(CCl4_CellProp) == "Freq"] <- "Proportion"

CCl4_CellProp <- CCl4_CellProp %>% filter(Cell_type != "Monocyte/T cell doublet" & Cell_type != "Red blood cell/myeloid doublet" &Cell_type != "Macrophage/endothelial Doublet")
CCl4_CellProp

plot_CCl4_CellPropBAR <- ggplot(CCl4_CellProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = CCl4.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_CellPropBAR_Celltype.png"), plot = plot_CCl4_CellPropBAR, width = 30, units = "cm")

plot_CCl4_CellPropBAR_Dataset <- ggplot(CCl4_CellProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_CellPropBAR_Dataset.png"), plot = plot_CCl4_CellPropBAR_Dataset, width = 30, units = "cm")


### Cell type membership vary by dataset ###

Idents(mA5) <- "Cell.type"
table(Idents(mA5))
mA5_CellNameProp <- prop.table(table(Idents(mA5), mA5$Names))
mA5_CellNameProp
mA5_DataProp <- as.data.frame(mA5_CellNameProp)

names(mA5_DataProp)[names(mA5_DataProp) == "Var1"] <- "Cell_type"
names(mA5_DataProp)[names(mA5_DataProp) == "Var2"] <- "Dataset"
names(mA5_DataProp)[names(mA5_DataProp) == "Freq"] <- "Proportion"

mA5_DataProp <- mA5_DataProp %>% filter(Cell_type != "Monocyte/T cell doublet" & Cell_type != "Neutrophil/T cell doublet" &Cell_type != "Neutrophil/B cell doublet")

plot_mA5_DataPropBAR <- ggplot(mA5_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = mA5.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_DataPropBAR_Celltype.png"), plot = plot_mA5_DataPropBAR, width = 30, units = "cm")

plot_mA5_DataPropBAR_Dataset <- ggplot(mA5_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_DataPropBAR_Dataset.png"), plot = plot_mA5_DataPropBAR_Dataset, width = 30, units = "cm")

#- 

Idents(Bleo2) <- "Cell.type"
table(Idents(Bleo2))
Bleo2_CellNameProp <- prop.table(table(Idents(Bleo2), Bleo2$Names))
Bleo2_CellNameProp
Bleo2_DataProp <- as.data.frame(Bleo2_CellNameProp)

names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Var1"] <- "Cell_type"
names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Var2"] <- "Dataset"
names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Freq"] <- "Proportion"

Bleo2_DataProp <- Bleo2_DataProp %>% filter(Cell_type != "T cell/B cell doublet")

plot_Bleo2_DataPropBAR <- ggplot(Bleo2_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = Bleo2.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_DataPropBAR_Celltype.png"), plot = plot_Bleo2_DataPropBAR, width = 30, units = "cm")

plot_Bleo2_DataPropBAR_Dataset <- ggplot(Bleo2_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_DataPropBAR_Dataset.png"), plot = plot_Bleo2_DataPropBAR_Dataset, width = 30, units = "cm")

#- 

Idents(CCl4) <- "Cell.type"
table(Idents(CCl4))
CCl4_CellNameProp <- prop.table(table(Idents(CCl4), CCl4$Names))
CCl4_CellNameProp
CCl4_DataProp <- as.data.frame(CCl4_CellNameProp)

names(CCl4_DataProp)[names(CCl4_DataProp) == "Var1"] <- "Cell_type"
names(CCl4_DataProp)[names(CCl4_DataProp) == "Var2"] <- "Dataset"
names(CCl4_DataProp)[names(CCl4_DataProp) == "Freq"] <- "Proportion"

CCl4_DataProp <- CCl4_DataProp %>% filter(Cell_type != "Monocyte/T cell doublet" & Cell_type != "Red blood cell/myeloid doublet" &Cell_type != "Macrophage/endothelial Doublet")

plot_CCl4_DataPropBAR <- ggplot(CCl4_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = CCl4.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_DataPropBAR_Celltype.png"), plot = plot_CCl4_DataPropBAR, width = 30, units = "cm")

plot_CCl4_DataPropBAR_Dataset <- ggplot(CCl4_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_DataPropBAR_Dataset.png"), plot = plot_CCl4_DataPropBAR_Dataset, width = 30, units = "cm")

  

  
# Individual Visualisations LINEAGE -> Proportion

### Dataset membership vary by lineage ###
Idents(mA5) <- "Names"
table(Idents(mA5))
mA5_NameLinProp <- prop.table(table(Idents(mA5), mA5$lineage))
mA5_LinProp <- as.data.frame(mA5_NameLinProp)

names(mA5_LinProp)[names(mA5_LinProp) == "Var1"] <- "Dataset"
names(mA5_LinProp)[names(mA5_LinProp) == "Var2"] <- "Lineage"
names(mA5_LinProp)[names(mA5_LinProp) == "Freq"] <- "Proportion"

mA5_LinProp <- mA5_LinProp %>% filter(Lineage != "Doublet")

mA5_plot1 <- ggplot(mA5_LinProp, aes(fill=Dataset, y=Proportion, x=Lineage)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = mA5.colpal)

ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_DataLinPropBAR.png"), plot = mA5_plot1, width = 30, units = "cm")

mA5_plot2 <- mA5_plot1 + coord_polar(theta = "x")
mA5_plot3 <- mA5_plot1 + coord_polar(theta = "y")
mA5_plot1_pie <- mA5_plot2 + mA5_plot3
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_DataLinPropPIE.png"), plot = mA5_plot1_pie, width = 40, units = "cm")

#- 
Idents(Bleo2) <- "Names"
table(Idents(Bleo2))
Bleo2_NameLinProp <- prop.table(table(Idents(Bleo2), Bleo2$lineage))
Bleo2_LinProp <- as.data.frame(Bleo2_NameLinProp)

names(Bleo2_LinProp)[names(Bleo2_LinProp) == "Var1"] <- "Dataset"
names(Bleo2_LinProp)[names(Bleo2_LinProp) == "Var2"] <- "Lineage"
names(Bleo2_LinProp)[names(Bleo2_LinProp) == "Freq"] <- "Proportion"

Bleo2_LinProp <- Bleo2_LinProp %>% filter(Lineage != "Doublet")

Bleo2_plot1 <- ggplot(Bleo2_LinProp, aes(fill=Dataset, y=Proportion, x=Lineage)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = Bleo2.colpal)

ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_DataLinPropBAR.png"), plot = Bleo2_plot1, width = 30, units = "cm")

Bleo2_plot2 <- Bleo2_plot1 + coord_polar(theta = "x")
Bleo2_plot3 <- Bleo2_plot1 + coord_polar(theta = "y")
Bleo2_plot1_pie <- Bleo2_plot2 + Bleo2_plot3
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_DataLinPropPIE.png"), plot = Bleo2_plot1_pie, width = 40, units = "cm")

#-
Idents(CCl4) <- "Names"
table(Idents(CCl4))
CCl4_NameLinProp <- prop.table(table(Idents(CCl4), CCl4$lineage))
CCl4_LinProp <- as.data.frame(CCl4_NameLinProp)

names(CCl4_LinProp)[names(CCl4_LinProp) == "Var1"] <- "Dataset"
names(CCl4_LinProp)[names(CCl4_LinProp) == "Var2"] <- "Lineage"
names(CCl4_LinProp)[names(CCl4_LinProp) == "Freq"] <- "Proportion"

CCl4_LinProp <- CCl4_LinProp %>% filter(Lineage != "Doublet")

CCl4_plot1 <- ggplot(CCl4_LinProp, aes(fill=Dataset, y=Proportion, x=Lineage)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = CCl4.colpal)

ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_DataLinPropBAR.png"), plot = CCl4_plot1, width = 30, units = "cm")

CCl4_plot2 <- CCl4_plot1 + coord_polar(theta = "x")
CCl4_plot3 <- CCl4_plot1 + coord_polar(theta = "y")
CCl4_plot1_pie <- CCl4_plot2 + CCl4_plot3
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_DataLinPropPIE.png"), plot = CCl4_plot1_pie, width = 40, units = "cm")
### Lineage membership vary by dataset ###

Idents(mA5) <- "lineage"
table(Idents(mA5))
mA5_LinNameProp <- prop.table(table(Idents(mA5), mA5$Names))
mA5_LinNameProp <- as.data.frame(mA5_LinNameProp)

names(mA5_LinNameProp)[names(mA5_LinNameProp) == "Var1"] <- "Lineage"
names(mA5_LinNameProp)[names(mA5_LinNameProp) == "Var2"] <- "Dataset"
names(mA5_LinNameProp)[names(mA5_LinNameProp) == "Freq"] <- "Proportion"

mA5_LinNameProp <- mA5_LinNameProp %>% filter(Lineage != "Doublet")

plot_LinNameBAR <- ggplot(mA5_LinNameProp, aes(fill=Lineage, y=Proportion, x=Dataset)) + geom_bar(position ="fill", stat="identity")
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_LinDataPropBAR.png"), plot = plot_LinNameBAR, width = 30, units = "cm")

plotLinNamePIE1 <- plot_LinNameBAR + coord_polar(theta = "x")
plotLinNamePIE2 <- plot_LinNameBAR + coord_polar(theta = "y") + scale_fill_manual(values = mA5.colpal)
plot_LinNamePIE <- plotLinNamePIE1 + plotLinNamePIE2
ggsave(filename = paste0(mA5Output.path, Today.date, "_mA5CLEANED_LinDataPropPIE.png"), plot = plot_LinNamePIE, width = 45, units = "cm")

#- 
Idents(Bleo2) <- "lineage"
table(Idents(Bleo2))
Bleo2_LinNameProp <- prop.table(table(Idents(Bleo2), Bleo2$Names))
Bleo2_LinNameProp <- as.data.frame(Bleo2_LinNameProp)

names(Bleo2_LinNameProp)[names(Bleo2_LinNameProp) == "Var1"] <- "Lineage"
names(Bleo2_LinNameProp)[names(Bleo2_LinNameProp) == "Var2"] <- "Dataset"
names(Bleo2_LinNameProp)[names(Bleo2_LinNameProp) == "Freq"] <- "Proportion"

Bleo2_LinNameProp <- Bleo2_LinNameProp %>% filter(Lineage != "Doublet")

plot_LinNameBAR <- ggplot(Bleo2_LinNameProp, aes(fill=Lineage, y=Proportion, x=Dataset)) + geom_bar(position ="fill", stat="identity")
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_LinDataPropBAR.png"), plot = plot_LinNameBAR, width = 30, units = "cm")

plotLinNamePIE1 <- plot_LinNameBAR + coord_polar(theta = "x")
plotLinNamePIE2 <- plot_LinNameBAR + coord_polar(theta = "y") + scale_fill_manual(values = Bleo2.colpal)
plot_LinNamePIE <- plotLinNamePIE1 + plotLinNamePIE2
ggsave(filename = paste0(Bleo2Output.path, Today.date, "_Bleo2CLEANED_LinDataPropPIE.png"), plot = plot_LinNamePIE, width = 45, units = "cm")

#-
Idents(CCl4) <- "lineage"
table(Idents(CCl4))
CCl4_LinNameProp <- prop.table(table(Idents(CCl4), CCl4$Names))
CCl4_LinNameProp <- as.data.frame(CCl4_LinNameProp)

names(CCl4_LinNameProp)[names(CCl4_LinNameProp) == "Var1"] <- "Lineage"
names(CCl4_LinNameProp)[names(CCl4_LinNameProp) == "Var2"] <- "Dataset"
names(CCl4_LinNameProp)[names(CCl4_LinNameProp) == "Freq"] <- "Proportion"

CCl4_LinNameProp <- CCl4_LinNameProp %>% filter(Lineage != "Doublet")

plot_LinNameBAR <- ggplot(CCl4_LinNameProp, aes(fill=Lineage, y=Proportion, x=Dataset)) + geom_bar(position ="fill", stat="identity")
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_LinDataPropBAR.png"), plot = plot_LinNameBAR, width = 30, units = "cm")

plotLinNamePIE1 <- plot_LinNameBAR + coord_polar(theta = "x")
plotLinNamePIE2 <- plot_LinNameBAR + coord_polar(theta = "y") + scale_fill_manual(values = CCl4.colpal)
plot_LinNamePIE <- plotLinNamePIE1 + plotLinNamePIE2
ggsave(filename = paste0(CCl4Output.path, Today.date, "_CCl4CLEANED_LinDataPropPIE.png"), plot = plot_LinNamePIE, width = 45, units = "cm")


# Remove Duplicates

# Merge + Save

# HVG 

# pre-harmony plots

# Harmony + Save

# post-harmony plots

# UMAP by group.by and split.by

# Feature Plots

# Cluster + Annotate + Save

