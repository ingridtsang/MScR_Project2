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
harmony.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Merge/"
harmonyOutput.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Merge/PresentationOutputs/"
BDL2Output.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Liver-BDL/Outputs/"
#Set date

Today.date <- "070820" 

# # # Annotate BDL
# BDL2 <- readRDS(paste0(BDL2.path, Today.date, "_BDL2_UMAP.RDS"))
# BDL2_cluster.annotation <- c("Liver Kupffer Cell 1",	"Blood Ly6Clo monocyte 1",	"BDL Blood Ly6Chi monocyte 1",	"Healthy Blood Ly6Chi monocyte 1",	"Liver macrophage 2",	"B cell 1",	"Liver macrophage 4",	"Liver macrophage 1",	"B cell 2",	"NK cell 1",	"BDL Blood Ly6Clo monocyte 1",	"CD8+ T cell 1",	"T cell 2",	"Cd4+ T cell 1",	"Proliferating",	"cDC2",	"NK cell 2",	"Neutrophil 1",	"Healthy Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Chi monocyte 2",	"Liver Endothelial Cells",	"Liver macrophage 3",	"Platelet 1",	"pDC",	"Blood Ly6Clo monocyte 2",	"B cell 3",	"Blood IFN-primed Ly6Chi monocyte",	"Liver Kupffer Cell 2",	"cDC1",	"Neutrophil 2",	"Liver macrophage 5",	"Monocyte/B cell doublet",	"Basophil",	"Monocyte/B cell doublet 1",	"Monocyte/T cell doublet 2",	"Platelet 2")

# names(BDL2_cluster.annotation) <- levels(BDL2)
# BDL2 <- RenameIdents(BDL2, BDL2_cluster.annotation)

# BDL2_cell.data <- data.table(barcode = colnames(BDL2),
                             # Cell.type = Idents(BDL2))

# BDL2_cell.data <- data.frame(BDL2_cell.data, row.names = BDL2_cell.data$barcode)
# BDL2_cell.data$barcode <- NULL
# BDL2 <- AddMetaData(BDL2, BDL2_cell.data, col.name = "Cell.type")

# ### save annotated UMAP ###
# annotated.umap.plot.split <- DimPlot(BDL2, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_cluster_OGIDsplit.png"), plot = annotated.umap.plot.split, width = 45, height = 10, units = "cm")

# annotated.umap.plot.merge <- DimPlot(BDL2,library(uwot)
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
harmony.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Merge/"
harmonyOutput.path <- "/home/shomea/m/mjdioli/Ingrid/Harmony/Harmony_Test2_SSH/Merge/PresentationOutputs/"
 reduction = "umap", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_cluster.png"), plot = annotated.umap.plot.merge, width = 30, units = "cm")


# ### cell lineage annotation ###
# BDL2_cell.data <- data.table(barcode = colnames(BDL2),
                             # Cell.type = Idents(BDL2))

# BDL2_lineage.annotation <- c("MP",	"MP",	"MP",	"MP",	"MP",	"B cell",	"MP",	"MP",	"B cell",	"T cell/ILC",	"MP",	"T cell/ILC",	"T cell/ILC",	"T cell/ILC",	"Prolferating",	"MP",	"MP",	"MP",	"MP",	"MP",	"Endothelia",	"MP",	"Platelet",	"MP",	"MP",	"B cell",	"MP",	"MP",	"MP",	"Neutrophil",	"MP",	"Doublet",	"Basophil",	"Doublet",	"Doublet",	"Platelet")


# BDL2_lineage.data <- data.table(Cell.type = BDL2_cluster.annotation, lineage = BDL2_lineage.annotation)

# meta.data <- merge(BDL2_cell.data, BDL2_lineage.data, by = "Cell.type")

# meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
# meta.data$barcode <- NULL
# meta.data$Cell.type <- NULL

# BDL2 <- AddMetaData(BDL2, meta.data, col.name = "lineage")

# # save annotated UMAP
# lineage.umap.plot <- DimPlot(BDL2, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 45, height = 10, units = "cm")

# lineage.umap.plot <- DimPlot(BDL2, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
# ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")

# saveRDS(BDL2, file = paste0(BDL2.path, Today.date, "_BDL2_annotated.RDS"))


# 
# #Load data
# 
# mA5_cleaned <- readRDS(paste0(mA5.path, Today.date, "_mA5_CLEANED.RDS"))
# 
# Bleo2_cleaned <- readRDS(paste0(Bleo2.path, Today.date, "_Bleo2_CLEANED.RDS"))
# 
# BDL2_cleaned <-readRDS(paste0(BDL2.path, Today.date, "_BDL2_CLEANED.RDS"))
# 
# CCl4_cleaned  <- readRDS(paste0(CCl4.path, Today.date, "_CCl4_CLEANED.RDS"))
# 
# mA5_harmony <- subset(mA5_cleaned, subset = Duplicate == 0)
# Bleo2_harmony <- subset(Bleo2_cleaned, subset = Duplicate == 0)
# BDL2_harmony <- subset(BDL2_cleaned, subset = Duplicate == 0)
# CCl4_harmony <- subset(CCl4_cleaned, subset = Duplicate == 0)
# 
# 
# # Merge + Save
# #harTest2 <- merge(x = mA5_harmony, y = c(Bleo2_harmony, BDL2_harmony, CCl4_harmony), add.cell.ids = c("Kidney_UUO", "Lung_Bleo", "Liver-BDL2", "Liver-CCl4"), project = "harTest2")
# #saveRDS(harTest2, file = paste0(harmony.path, "120820_preRH_Test2.RDS"))
# harTest2 <- readRDS(paste0(harmony.path, "120820_preRH_Test2.RDS"))
# 
# # HVG 
# mA5_harmony <- FindVariableFeatures(mA5_harmony, selection.method = "vst", nfeatures = 5000)
# Bleo2_harmony <- FindVariableFeatures(Bleo2_harmony, selection.method = "vst", nfeatures = 5000)
# BDL2_harmony <- FindVariableFeatures(BDL2_harmony, selection.method = "vst", nfeatures = 5000)
# CCl4_harmony <- FindVariableFeatures(CCl4_harmony, selection.method = "vst", nfeatures = 5000)
# 
# Test2_HVG <- intersect(VariableFeatures(mA5_harmony), intersect(VariableFeatures(Bleo2_harmony), intersect(VariableFeatures(BDL2_harmony), VariableFeatures(CCl4_harmony))))
# VariableFeatures(harTest2) <- Test2_HVG
# 
# 
# # pre-harmony + plots
# #harTest2 <- readRDS(paste0(harmony.path, "120820_preRH_Test2.RDS"))
# harTest2 <- ScaleData(harTest2)
# harTest2 <- RunPCA(harTest2, features = VariableFeatures(object = harTest2))
# saveRDS(harTest2, file = paste0(harmony.path, "120820_preRH_Test2_PCA.RDS"))
# 
# plot_preRH_Elbow <- ElbowPlot(harTest2, ndims = 50, reduction = "pca")
# ggsave(filename = paste0(harmony.path, "120820_Test2_Elbow_preRH.png"), plot = plot_preRH_Elbow, width = 30, units = "cm")
# 
# plot_preRH_Dim <- DimPlot(object = harTest2, reduction = "pca", group.by = "Names")
# ggsave(filename = paste0(harmony.path, "120820_Test2_DimPlot_preRH.png"), plot = plot_preRH_Dim, width = 30, units = "cm")
# 
# # Harmony + Save
# 
# harTest2 <- RunHarmony(harTest2, group.by.vars = "orig.ident")
# make.names("harTest2") 
# saveRDS(harTest2, file = paste0(harmony.path, "120820_HarmonyTest2.RDS"))
# 
# # post-harmony plots
# 
# # UMAP + plot group.by and split.by
# 
# harTest2  <- FindNeighbors(harTest2, reduction = "harmony", dims = 1:50)
# harTest2 <- FindClusters(harTest2, reduction = "harmony", resolution = 1.00)
# harTest2_UMAP <- RunUMAP(harTest2, reduction = "harmony", dims = 1:50)
# saveRDS(harTest2_UMAP, file = paste0(harmony.path, "120820_HarmonyTest2_UMAP.RDS"))

harTest2_UMAP <- readRDS(paste0(harmony.path, "120820_HarmonyTest2_UMAP.RDS"))
har_UMAP_OGID <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "Names")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_OGID.png"), plot = har_UMAP_OGID, width = 30, units = "cm")

har_UMAP_OGID_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "Names", split.by = "Experiment")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_OGIDsplitExp.png"), plot = har_UMAP_OGID_split, width = 30, units = "cm")

har_UMAP_cl <- DimPlot(harTest2_UMAP, reduction = "umap", label = T)
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_CL.png"), plot = har_UMAP_cl, width = 30, units = "cm")

har_UMAP_cl_split <- DimPlot(harTest2_UMAP, reduction = "umap", split.by = "Experiment")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_CLsplitExp.png"), plot = har_UMAP_cl_split, width = 30, units = "cm")

har_UMAP_cl_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "cluster", split.by = "Injury.type")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_CLsplitInjT.png"), plot = har_UMAP_cl_split, width = 30, units = "cm")

har_UMAP_cl_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "cluster", split.by = "Tissue")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_CLsplitTiss.png"), plot = har_UMAP_cl_split, width = 30, units = "cm")

har_UMAP_cl_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "cluster", split.by = "Injury")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_CLsplitInj.png"), plot = har_UMAP_cl_split, width = 30, units = "cm")

har_UMAP_lin <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "lineage")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_LIN.png"), plot = har_UMAP_lin, width = 30, units = "cm")

har_UMAP_lin_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "lineage", split.by = "Experiment")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_LINsplitExp.png"), plot = har_UMAP_lin_split, width = 30, units = "cm")

har_UMAP_lin_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "lineage", split.by = "Injury.type")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_LINsplitInjT.png"), plot = har_UMAP_lin_split, width = 30, units = "cm")

har_UMAP_lin_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "lineage", split.by = "Tissue")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_LINsplitTiss.png"), plot = har_UMAP_lin_split, width = 30, units = "cm")

har_UMAP_lin_split <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "lineage", split.by = "Injury")
ggsave(filename = paste0(harmonyOutput.path, "120820_Test2_UMAP_LINsplitInj.png"), plot = har_UMAP_lin_split, width = 30, units = "cm")


# # Marker Genes + Feature Plots
# 
# harmony.markers <- FindAllMarkers(harTest2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(harmony.markers, file = paste0(harmony.path, "120820_Test2_ALLmarkers.csv"))
# harmony.markersOUTPUT <- harmony.markers %>% group_by(cluster) %>% arrange(cluster, pct.2) 
# write.csv(harmony.markersOUTPUT, file = paste0(harmony.path, "120820_Test2_biomarkers.csv"))