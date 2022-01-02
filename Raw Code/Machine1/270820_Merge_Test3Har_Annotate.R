# Title     : TODO
# Objective : TODO
# Created by: mjdioli
# Created on: 28/08/2020

#Load libraries
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)
library(data.table)

#Set dates
Load.date <- "070820"
Save.date <- "160820"
date.1 <- "200820"
date.2 <- "210820"
date.3 <- "220820"
date.4 <- "240820"

#Set paths
merge.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/"
Experiment.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/RegressExp/"

#Load Harmonised data
Test3 <- readRDS(paste0(merge.path, date.3,"_Harmony_Test3.RDS"))

#Annotate Cell-type
Test3_cluster.annotation <- c("B cell 1", "Ly6C lo monocyte 1", "Ly6C hi monocyte", "Liver Kupffer cell 1",
                              "Fibrosis-associated Macrophage 1", "Fibrosis-associated Ly-6Chi monocyte",
                              "Macrophage 2", "CD8+ T cell 1", "NK cell 1", "Cd4+ T cell 1", "cDC2",
                              "Fibrosis-associated Macrophage 2", "Fibrosis-associated Macrophage 3", "Neutrophil",
                              "NKT cell", "Proliferating 1", "pDC", "cDC1", "Kidney resident macrophage",
                              "Ly6C lo monocyte 2", "IFN-primed monocyte", "Cd4+ T cell 2", "Lung alveolar macrophage",
                              "Endothelial cells ", "NK cell 2", "B cell 2", "Proliferating 2", "Platelet 1",
                              "Ly6C lo monocyte 3", "Liver Kupffer cell 2", "Fibrosis-associated Macrophage 4",
                              "Basophil", "CCR7+ cDC", "Platelet 2", "T cell/ B cell doublet", "B cell 3"
)

 names(Test3_cluster.annotation) <- levels(Test3)
 Test3 <- RenameIdents(Test3, Test3_cluster.annotation)

 Test3_cell.data <- data.table(barcode = colnames(Test3),
                             Cell.type = Idents(Test3))

Test3_cell.data <- data.frame(Test3_cell.data, row.names = Test3_cell.data$barcode)
Test3_cell.data$barcode <- NULL
Test3 <- AddMetaData(Test3, Test3_cell.data, col.name = "Cell.type")


#Annotate Lineage
Test3_cell.data <- data.table(barcode = colnames(Test3),
                             Cell.type = Idents(Test3))

Test3_lineage.annotation <- c("B cell", "MP", "MP", "MP",
                              "MP", "MP",
                              "MP", "T cell/ILC", "T cell/ILC", "T cell/ILC", "MP",
                              "MP", "MP", "Neutrophil",
                              "T cell/ILC", "Proliferating", "MP", "MP", "MP",
                              "MP", "MP", "T cell/ILC", "MP",
                              "Enodthelia", "T cell/ILC", "B cell", "Proliferating", "Platelet",
                              "MP", "MP", "MP",
                              "Basophil", "MP", "Platelet", "Doublet", "B cell"
)



Test3_lineage.data <- data.table(Cell.type = Test3_cluster.annotation, lineage = Test3_lineage.annotation)
meta.data <- merge(Test3_cell.data, Test3_lineage.data, by = "Cell.type")

meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$Cell.type <- NULL

Test3 <- AddMetaData(Test3, meta.data, col.name = "lineage")

#Save
saveRDS(Test3, file = paste0(merge.path, date.4, "_Harmony_aTest3.RDS"))

#Graphs
cluster.umap.plot <- DimPlot(Test3, reduction = "umap", group.by = "Cell.type", split.by = "Tissue", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_cl_Tissplit.png"), plot = cluster.umap.plot, width = 50, height = 10, units = "cm")

cluster.umap.plot <- DimPlot(Test3, reduction = "umap", group.by = "Cell.type", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_cl.png"), plot = cluster.umap.plot, width = 30, units = "cm")

lineage.umap.plot <- DimPlot(Test3, reduction = "umap", group.by = "lineage", split.by = "Tissue", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 50, height = 10, units = "cm")

lineage.umap.plot <- DimPlot(Test3, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")



Test3_cleaned <- subset(Test3, lineage != "Doublet")
saveRDS(Test3_cleaned, file = paste0(Experiment.path, date.4, "_Harmony_aTest3CLEANED.RDS"))



Idents(Test3_cleaned) <- "Cell.type"
ALL_Conserved <- c("Fibrosis-associated Macrophage 1", "Macrophage 2", "Fibrosis-associated Macrophage 3")

lineage.umap.plot <- DimPlot(Test3_cleaned, reduction = "umap", group.by = "lineage", split.by = "Condition", cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage_Condsplit.png"), plot = lineage.umap.plot, width = 40, units = "cm")

lineage.umap.plot <- DimPlot(Test3_cleaned, reduction = "umap", group.by = "lineage", split.by = "Tissue", cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659",  label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage_Tissplit.png"), plot = lineage.umap.plot, width = 50, height = 10, units = "cm")

lineage.umap.plot <- DimPlot(Test3_cleaned, reduction = "umap", group.by = "lineage", label = T, cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659",  label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")

Test3MP_cleaned <- subset(Test3, subset = lineage == "MP")

lineage.umap.plot <- DimPlot(Test3MP_cleaned, reduction = "umap", group.by = "Cell.type", split.by = "Condition", cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_cell_Condsplit_MP.png"), plot = lineage.umap.plot, width = 40, units = "cm")

lineage.umap.plot <- DimPlot(Test3MP_cleaned, reduction = "umap", group.by = "Cell.type", split.by = "Tissue", cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659",  label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage_Tissplit_MP.png"), plot = lineage.umap.plot, width = 50, height = 15, units = "cm")

lineage.umap.plot <- DimPlot(Test3MP_cleaned, reduction = "umap", group.by = "Cell.type", label = T, cells.highlight = WhichCells(Test3_cleaned, idents = ALL_Conserved), cols.highlight = "#FF5659",  label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(Experiment.path,date.4 , "_MERGE_aUMAP_lineage_MP.svg"), plot = lineage.umap.plot, width = 30, units = "cm")
