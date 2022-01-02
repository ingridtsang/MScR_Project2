#Load libraries
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)

#Set dates
Load.date <- "070820"
Save.date <- "160820"
date.1 <- "200820"
date.2 <- "210820"
date.3 <- "220820"

#Set paths
mA5origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
BDL2origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
CCl4origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"

mA5Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Kidney_UUO/"
Bleo2Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Lung_Bleo/"
BDL2Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_BDL/"
CCl4Output.path  <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_CCl4/"

merge.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/"
Experiment.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/RegressExp/"

#Load data
amA5_cleaned <- readRDS(paste0(mA5Output.path, date.1, "_amA5_cleaned.RDS"))

aBleo2_cleaned <- readRDS(paste0(Bleo2Output.path, date.1, "_aBleo2_cleaned.RDS"))

aBDL2_cleaned <- readRDS(paste0(BDL2Output.path, date.1, "_aBDL2_cleaned.RDS"))

aCCl4_cleaned <- readRDS(paste0(CCl4Output.path, date.1, "_aCCl4_cleaned.RDS"))


mA5_harmony <- subset(amA5_cleaned, subset = Duplicate == 0)
Bleo2_harmony <- subset(aBleo2_cleaned, subset = Duplicate == 0)
BDL2_harmony <- subset(aBDL2_cleaned, subset = Duplicate == 0)
CCl4_harmony <- subset(aCCl4_cleaned, subset = Duplicate == 0)
rm(amA5_cleaned, aBleo2_cleaned, aBDL2_cleaned, aCCl4_cleaned)

#
mA5_harmony <- FindVariableFeatures(mA5_harmony, selection.method = "vst", nfeatures = 4000)
Bleo2_harmony <- FindVariableFeatures(Bleo2_harmony, selection.method = "vst", nfeatures = 4000)
BDL2_harmony <- FindVariableFeatures(BDL2_harmony, selection.method = "vst", nfeatures = 4000)
CCl4_harmony <- FindVariableFeatures(CCl4_harmony, selection.method = "vst", nfeatures = 4000)

Test3 <- merge(x = mA5_harmony, y = c(Bleo2_harmony, BDL2_harmony, CCl4_harmony), add.cell.ids = c("Kidney_UUO", "Lung_Bleo", "Liver-BDL2", "Liver-CCl4"), project = "Test3")
saveRDS(Test3, file = paste0(merge.path, date.3,"_preRH_Test3.RDS"))

Test3_HVG <- intersect(VariableFeatures(mA5_harmony), intersect(VariableFeatures(Bleo2_harmony), intersect(VariableFeatures(BDL2_harmony), VariableFeatures(CCl4_harmony))))
VariableFeatures(Test3) <- Test3_HVG


rm(mA5_harmony, Bleo2_harmony, BDL2_harmony, CCl4_harmony)
# pre-harmony + plots

Test3 <- ScaleData(Test3)
Test3 <- RunPCA(Test3, features = VariableFeatures(object = Test3))
saveRDS(Test3, file = paste0(merge.path, date.3,"_preRH_Test3_PCA.RDS"))

plot_preRH_Elbow <- ElbowPlot(Test3, ndims = 50, reduction = "pca")
ggsave(filename = paste0(merge.path, date.3, "_preRH_Test3_Elbow.png"), plot = plot_preRH_Elbow, width = 30, units = "cm")

plot_preRH_Dim <- DimPlot(Test3, reduction = "pca", group.by = "Names")
ggsave(filename = paste0(merge.path, date.3, "_preRH_Test3_DimPlot.png"), plot = plot_preRH_Dim, width = 30, units = "cm")

Test3 <- readRDS(paste0(merge.path, date.3,"_preRH_Test3_PCA.RDS"))
Test3_UMAP  <- FindNeighbors(Test3, reduction = "pca", dims = 1:50)
Test3_UMAP <- FindClusters(Test3_UMAP, resolution = 1.0)
Test3_UMAP <- RunUMAP(Test3_UMAP, reduction = "pca", dims = 1:50)

plot_preRH_UMAP <- DimPlot(Test3_UMAP, reduction = "umap", label = T) + NoLegend()
ggsave(filename = paste0(merge.path, date.3, "_preRH_Test3_UMAP_cl.png"), width = 30, units = "cm")

plot_preRH_UMAP_Exp <- DimPlot(Test3_UMAP, reduction = "umap", group.by = "Experiment")
ggsave(filename = paste0(merge.path, date.3, "_preRH_Test3_UMAP_OGID.png"), width = 30, units = "cm")




#v <- c("orig.ident", "Tissue", "Injury", "Experiment")
#group.by.vars <- v
#for (v in group.by.vars){
Test3 <- readRDS(paste0(merge.path, date.3,"_preRH_Test3_PCA.RDS"))

v <- "Experiment"

Test3 <- RunHarmony(Test3, group.by.vars = v)

Test3  <- FindNeighbors(Test3, reduction = "harmony", dims = 1:50)
Test3 <- FindClusters(Test3, resolution = 1.0)
Test3 <- RunUMAP(Test3, reduction = "harmony", dims = 1:50)
saveRDS(Test3, file = paste0(merge.path, date.3,"_Harmony_Test3.RDS"))

  har_UMAP_cl <- DimPlot(Test3, reduction = "umap", group.by = "Experiment", label = T)
  ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPcl.png"), plot = har_UMAP_cl, width = 30, units = "cm")

Test3.markers <- FindAllMarkers(Test3, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Test3.markers, file = paste0(Experiment.path, date.3, "_", v, "_ALLmarkers.csv"))
Test3OUTPUT <- Test3.markers %>% group_by(cluster) %>% arrange(desc(p_val_adj), pct.2, .by_group = T)
write.csv(Test3OUTPUT, file = paste0(Experiment.path, date.3, "_", v, "_biomarkers.csv"))


Test3 <- readRDS(paste0(merge.path, date.3,"_Harmony_Test3.RDS"))

v = "Experiment"
harMP_UMAP <- DimPlot(Test3, split.by = "Experiment", reduction = "umap", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPclExp.png"), plot = harMP_UMAP, width = 50, height= 10, units = "cm")

harMP_UMAP <- DimPlot(Test3, group.by  = "Experiment", reduction = "umap")
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPExp.png"), plot = harMP_UMAP, width = 30, units = "cm")


har_UMAP_cl <- DimPlot(Test3, reduction = "umap", label = T)
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPcl.png"), plot = har_UMAP_cl, width = 40, units = "cm")

har_UMAP_cond <- DimPlot(Test3, reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPCond.png"), plot = har_UMAP_cond, width = 40, units = "cm")

har_UMAP_condmA5 <- DimPlot(subset(Test3, subset = Experiment == "Kidney.UUO"), reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPCond_UUO.png"), plot = har_UMAP_condmA5, width = 40, units = "cm")

har_UMAP_condBleo <- DimPlot(subset(Test3, subset = Experiment == "Lung.Bleo"), reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPCond_Bleo.png"), plot = har_UMAP_condBleo, width = 40, units = "cm")

har_UMAP_condBDL <- DimPlot(subset(Test3, subset = Experiment == "Liver.BDL"), reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_", v, "_UMAPCond_BDL.png"), plot = har_UMAP_condBDL, width = 40, units = "cm")

har_UMAP_condCCl4 <- DimPlot(subset(Test3, subset = Experiment == "Liver.CCl4"), reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPCond_CCl4.png"), plot = har_UMAP_condCCl4, width = 40, units = "cm")

harMP_UMAP <- DimPlot(subset(Test3, subset = lineage == "MP"), reduction = "umap", label = T, label.size = 2)  + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_MP_", v, "_UMAPcl.png"), plot = harMP_UMAP, width = 30, units = "cm")

harMP_UMAP <- DimPlot(subset(Test3, subset = lineage == "MP"), split.by = "Experiment", reduction = "umap", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_MP_", v, "_UMAPclExp.png"), plot = harMP_UMAP, width = 50, height= 10, units = "cm")

harMP_UMAP <- DimPlot(subset(Test3, subset = lineage == "MP"), group.by  = "Experiment", reduction = "umap")
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_MP_", v, "_UMAPExp.png"), plot = harMP_UMAP, width = 30, units = "cm")

harMP_UMAP <- DimPlot(subset(Test3, subset = lineage == "MP"), group.by = "Cell.type", reduction = "umap", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_MP_", v, "_UMAPCell.png"), plot = harMP_UMAP, width = 50, units = "cm")

harMP_UMAP <- DimPlot(subset(Test3, subset = lineage == "MP"), group.by = "Cell.type", split.by = "Experiment", reduction = "umap", label = T, label.size = 3) + NoLegend()
ggsave(filename = paste0(Experiment.path, date.3, "_harPara_MP_", v, "_UMAPCellExp.png"), plot = harMP_UMAP, width = 70, height = 20, units = "cm")


--
har_UMAP_cond <- DimPlot(Test3, reduction = "umap", split.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPclCond.png"), plot = har_UMAP_cond, width = 50, units = "cm")

har_UMAP_cond <- DimPlot(Test3, reduction = "umap", group.by = "Condition", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPCond.png"), plot = har_UMAP_cond, width = 40, units = "cm")

har_UMAP_cond <- DimPlot(Test3, reduction = "umap", split.by = "Tissue", label = T) + NoLegend()
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPclTiss.png"), plot = har_UMAP_cond, width =50, height = 10, units = "cm")

har_UMAP_cond <- DimPlot(Test3, reduction = "umap", group.by = "Tissue", label = F)
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPTiss.png"), plot = har_UMAP_cond, width = 40, units = "cm")

har_UMAP_cond <- DimPlot(Test3, reduction = "umap", split.by = "Tissue", group.by = "Condition", label = F)
ggsave(filename = paste0(Experiment.path,  date.3, "_harPara_", v, "_UMAPCondTiss.png"), plot = har_UMAP_cond, width =50, height = 10, units = "cm")
