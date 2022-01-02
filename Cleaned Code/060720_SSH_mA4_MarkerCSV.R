library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(uwot)

# Load dataset
#mA4 <- readRDS("./060720/mA4.RDS")
#mA4

# Load and merge individual datasets
MKhealthyrmD <- readRDS(file = "./060720/mA4Indie_MKhealthy_rmD.RDS")
MBhealthyrmD <- readRDS(file = "./060720/mA4Indie_MBhealthy_rmD.RDS")

MKinjuredrmD <- readRDS(file = "./060720/mA4Indie_MKinjured_rmD.RDS")
MBinjuredrmD <- readRDS(file = "./060720/mA4Indie_MBinjured_rmD.RDS")

mA4 <- merge(MKhealthyrmD, y = c(MBhealthyrmD, MKinjuredrmD, MBinjuredrmD), add.cell.ids = c("H-kidney", "H-blood", "Inj-kidney", "Inj-blood"), project = "mA4")
unique(sapply(X = strsplit(colnames(mA4), split = "_"), FUN = "[", 1))

# QC, LogNorm and Scale data
mA4[["percent.mt"]] <- PercentageFeatureSet(mA4, pattern = "^mt-")
plot_QC_Vlnplot <- VlnPlot(mA4, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "./060720/mA4_QCviolin.png", plot = plot_QC_Vlnplot, width = 25, units = "cm")

plot_QC_Scatter1 <- FeatureScatter(mA4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_QC_Scatter2 <- FeatureScatter(mA4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2
ggsave(filename = "./060720/mA4_QCscatter.png", plot = plot_QC_Scatter, width = 25, units = "cm")

mA4 <- ScaleData(mA4)

# PCA + UMAP
mA4 <- FindVariableFeatures(object = mA4)
mA4 <- RunPCA(mA4, features = VariableFeatures(object = mA4))
mA4 <- FindNeighbors(mA4, dims = 1:43)
mA4 <- FindClusters(mA4, resolution = 1.0)
mA4 <- RunUMAP(mA4, dims = 1:43)

saveRDS(mA4, file = "./060720/mA4.RDS")
# UMAP plot of mA4 by sample type
plot_mA4_tissue <- DimPlot(mA4, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "./060720/mA4_UMAP_tissue.jpeg", plot = plot_mA4_tissue, width = 25, units = "cm")


# UMAP plot of mA4 by cluster

plot_mA4_cluster <- DimPlot(mA4, reduction = "umap")
ggsave(filename = "./060720/mA4_UMAP_cluster.jpeg", plot = plot_mA4_cluster, width = 25, units = "cm")

# Create marker gene CSV

mA4.markers <- FindAllMarkers(mA4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mA4.markersOUTPUT <- mA4.markers %>% group_by(cluster)
write.csv(mA4.markersOUTPUT, file = "./060720/K_mA4_biomarkers.csv")

