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

mA5.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
BDL2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
CCl4.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"
harmony.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Merge/"

#Load data

mA5_cleaned <- readRDS(paste0(mA5.path, Today.date, "_mA5_CLEANED.RDS"))

Bleo2_cleaned <- readRDS(paste0(Bleo2.path, Today.date, "_Bleo2_CLEANED.RDS"))

BDL2_cleaned <- readRDS(paste0(BDL2Output.path, Today.date, "_BDL2_CLEANED.RDS"))

CCl4_cleaned  <- readRDS(paste0(CCl4.path, Today.date, "_CCl4_CLEANED.RDS"))

mA5_harmony <- subset(mA5_cleaned, subset = Duplicate == 0)
Bleo2_harmony <- subset(Bleo2_cleaned, subset = Duplicate == 0)
BDL2_harmony <- subset(BDL2_cleaned, subset = Duplicate == 0)
CCl4_harmony <- subset(CCl4_cleaned, subset = Duplicate == 0)

# HVG 
mA5_harmony <- FindVariableFeatures(mA5_harmony, selection.method = "vst", nfeatures = 5000)
Bleo2_harmony <- FindVariableFeatures(Bleo2_harmony, selection.method = "vst", nfeatures = 5000)
BDL2_harmony <- FindVariableFeatures(BDL2_harmony, selection.method = "vst", nfeatures = 5000)
CCl4_harmony <- FindVariableFeatures(CCl4_harmony, selection.method = "vst", nfeatures = 5000)

Test2_HVG <- intersect(VariableFeatures(mA5_harmony), intersect(VariableFeatures(Bleo2_harmony), intersect(VariableFeatures(BDL2_harmony), VariableFeatures(CCl4_harmony))))
VariableFeatures(harTest2) <- Test2_HVG


# pre-harmony + plots
harTest2 <- readRDS(paste0(harmony.path, "120820_preRH_Test2.RDS"))
harTest2 <- ScaleData(harTest2)
harTest2 <- RunPCA(harTest2, features = VariableFeatures(object = harTest2))
saveRDS(harTest2, file = paste0(harmony.path, "120820_preRH_Test2_PCA.RDS"))

plot_preRH_Elbow <- ElbowPlot(harTest2, ndims = 50, reduction = "pca")
ggsave(filename = paste0(harmony.path, "120820_Test2_Elbow_preRH.png"), plot = plot_preRH_Elbow, width = 30, units = "cm")

plot_preRH_Dim <- DimPlot(object = harTest2, reduction = "pca", group.by = "Names")
ggsave(filename = paste0(harmony.path, "120820_Test2_DimPlot_preRH.png"), plot = plot_preRH_Dim, width = 30, units = "cm")

# Harmony + Save
make.names("harTest2")
harTest2 <- RunHarmony(harTest2, group.by.vars = "orig.ident", dims = 1:50)
saveRDS(harTest2, file = paste0(harmony.path, "120820_HarmonyTest2.RDS"))

# post-harmony plots

# UMAP + plot group.by and split.by

harTest2  <- FindNeighbors(harTest2, reduction = "harmony", dims = 1:50)
harTest2 <- FindClusters(harTest2, reduction = "harmony", resolution = 1.00)
harTest2_UMAP <- RunUMAP(harTest2, reduction = "harmony", dims = 1:50)
saveRDS(harTest2_UMAP, file = paste0(harmony.path, "120820_HarmonyTest2_UMAP.RDS"))

har_UMAP_OGID <- DimPlot(harTest2_UMAP, reduction = "UMAP", group.by = "Names")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_OGID.png"), plot = har_UMAP_cluster, width = 30, units = "cm")

har_UMAP_OGID_split <- DimPlot(harTest2_UMAP, reduction = "UMAP", group.by = "Names", split.by = "Names")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_OGIDsplit.png"), plot = har_UMAP_OGID_split, width = 30, units = "cm")

har_UMAP_cl <- DimPlot(harTest2_UMAP, reduction = "UMAP")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_CL.png"), plot = har_UMAP_cl, width = 30, units = "cm")

har_UMAP_cl_split <- DimPlot(harTest2_UMAP, reduction = "UMAP", split.by = "Names")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_CLsplit.png"), plot = har_UMAP_cl_split, width = 30, units = "cm")

har_UMAP_lin <- DimPlot(harTest2_UMAP, reduction = "UMAP", group.by = "lineage")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_LIN.png"), plot = har_UMAP_lin, width = 30, units = "cm")

har_UMAP_lin_split <- DimPlot(harTest2_UMAP, reduction = "UMAP", group.by = "lineage", split.by = "Names")
ggsave(filename = paste0(harmony.path, "120820_Test2_UMAP_LINsplit.png"), plot = har_UMAP_lin_split, width = 30, units = "cm")

# # Marker Genes + Feature Plots

harmony.markers <- FindAllMarkers(harTest2, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(harmony.markers, file = paste0(harmony.path, "120820_Test2_ALLmarkers.csv"))
harmony.markersOUTPUT <- harmony.markers %>% group_by(cluster) %>% arrange(cluster, pct.2) 
write.csv(harmony.markersOUTPUT, file = paste0(harmony.path, "120820_Test2_biomarkers.csv"))
