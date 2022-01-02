# Load library
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

mA5origin.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2origin.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
BDL2origin.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
CCl4origin.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"

mA5Output.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test3/Kidney_UUO/"
Bleo2Output.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test3/Lung_Bleo/"
BDL2Output.path <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_BDL/"
CCl4Output.path  <- "/Users/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_CCl4/"

#Set date

Load.date <- "070820"
Save.date <- "160820"
date.1 <- "200820"

#Load annotated datasets
mA5 <- readRDS(paste0(mA5origin.path, Load.date, "_mA5_UMAP.RDS"))
Bleo2 <- readRDS(paste0(Bleo2origin.path, Load.date, "_Bleo2_UMAP.RDS"))
BDL2 <- readRDS(paste0(BDL2origin.path, Load.date, "_BDL2_UMAP.RDS"))
CCl4 <- readRDS(paste0(CCl4origin.path, Load.date, "_CCl4_UMAP.RDS"))

mA5_lin.colpal <- c("#ffe5cc",	"#ffcb99",	"#ffb266",	"#ff9832",	"#ff7f00",	"#cc6500",	"#994c00",	"#663200",	"#331900")
Bleo2_lin.colpal <-	c("#a58ac2",	"#9062c1",	"#7744ac",	"#5d3688",	"#442763",	"#3f245c")
BDL2_lin.colpal <- c("#d6ecd4",	"#add9aa",	"#84c680",	"#5bb356",	"#33a02c",	"#2d9027",	"#23701e",	"#195016",	"#0f300d")
CCl4_lin.colpal <- c("#d0e5f7",	"#a1cbf0",	"#73b1e9",	"#4497e2",	"#167edb",	"#1371c5",	"#0f5899",	"#0b3f6d",	"#062541")


#Load marker gene list
mA5.markers <- read.csv(paste0(mA5origin.path, "Outputs/" , Load.date, "_mA5_ALLmarkers.csv"))
Bleo2.markers <- read.csv(paste0(Bleo2origin.path,"Outputs/", Load.date, "_Bleo2_ALLmarkers.csv"))
BDL2.markers <- read.csv(paste0(BDL2origin.path,"Outputs/", Load.date, "_BDL2_ALLmarkers.csv"))
CCl4.markers <- read.csv(paste0(CCl4origin.path,"Outputs/", Load.date, "_CCl4_ALLmarkers.csv"))

#HM cluster
mA5.markersTop10 <- mA5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
mA5.markersTop3 <- mA5.markersTop10 %>% top_n(n = 3, wt = pct.2)
print(mA5.markersTop3$gene)
mA5_cellHM <- DoHeatmap(mA5, group.by = "seurat_clusters", features = mA5.markersTop3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(mA5Output.path, date.1, "_mA5_Supp_clusterHM.png"), plot = mA5_cellHM, width = 45, units = "cm")

Bleo2.markersTop10 <- Bleo2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
Bleo2.markersTop3 <- Bleo2.markersTop10 %>% top_n(n = 3, wt = pct.2)
print(Bleo2.markersTop3$gene)
Bleo2_cellHM <- DoHeatmap(Bleo2, group.by = "seurat_clusters", features = Bleo2.markersTop3$gene, assay = "RNA", slot = "data", label = T,size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(Bleo2Output.path, date.1, "_Bleo2_Supp_clusterHM.png"), plot = Bleo2_cellHM, width = 45, units = "cm")

BDL2.markersTop10 <- BDL2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
BDL2.markersTop3 <- BDL2.markersTop10 %>% top_n(n = 3, wt = pct.2)
print(BDL2.markersTop3$gene)
BDL2_cellHM <- DoHeatmap(BDL2, group.by = "seurat_clusters", features = BDL2.markersTop3$gene, assay = "RNA", slot = "data", label = T,size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, date.1, "_BDL2_Supp_clusterHM.png"), plot = BDL2_cellHM, width = 45, units = "cm")

CCl4.markersTop10 <- CCl4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
CCl4.markersTop3 <- CCl4.markersTop10 %>% top_n(n = 3, wt = pct.2)
print(CCl4.markersTop3$gene)
CCl4_cellHM <- DoHeatmap(CCl4, group.by = "seurat_clusters", features = CCl4.markersTop3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(CCl4Output.path, date.1, "_CCl4_Supp_clusterHM.png"), plot = CCl4_cellHM, width = 45, units = "cm")
