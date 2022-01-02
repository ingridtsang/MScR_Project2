

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(uwot)


# Load data
MliHealthy_A.data <- Read10X(data.dir = "./Filtered Matrices/H-liverA/")
MliHealthy_A <- CreateSeuratObject(counts = MliHealthy_A.data, project = "H-liver_A")
saveRDS(MliHealthy_A, file = "./100720/MLiH_OG_A.RDS")

MliHealthy_B.data <- Read10X(data.dir = "./Filtered Matrices/H-liverB/")
MliHealthy_B <- CreateSeuratObject(counts = MliHealthy_B.data, project = "H-liver_B")
saveRDS(MliHealthy_B, file = "./100720/MLiH_OG_B.RDS")

MBhealthy_A.data <- Read10X(data.dir = "./Filtered Matrices/H-bloodA/")
MBhealthy_A <- CreateSeuratObject(counts = MBhealthy_A.data, project = "H-blood_A")
saveRDS(MBhealthy_A, file = "./100720/MBHealthy_OG_A.RDS")

MBhealthy_B.data <- Read10X(data.dir = "./Filtered Matrices/H-bloodB/")
MBhealthy_B <- CreateSeuratObject(counts = MBhealthy_B.data, project = "H-blood_B")
saveRDS(MBhealthy_B, file = "./100720/MBhealthyOG_B.RDS")

MliBDL_A.data <- Read10X(data.dir = "./Filtered Matrices/BDL-liverA/")
MliBDL_A <- CreateSeuratObject(counts = MliBDL_A.data, project = "BDL-liver_A")
saveRDS(MliBDL_A, file = "./100720/MLiBDL_OG_A.RDS")

MliBDL_B.data <- Read10X(data.dir = "./Filtered Matrices/BDL-liverB/")
MliBDL_B <- CreateSeuratObject(counts = MliBDL_B.data, project = "BDL-liver_B")
saveRDS(MliBDL_B, file = "./100720/MLiBDL_OG_B.RDS")

MB_BDL_A.data <- Read10X(data.dir = "./Filtered Matrices/BDL-bloodA/")
MB_BDL_A <- CreateSeuratObject(counts = MB_BDL_A.data, project = "BDL-blood_A")
saveRDS(MB_BDL_A, file = "./100720/MBBDL_OG_A.RDS")

MB_BDL_B.data <- Read10X(data.dir = "./Filtered Matrices/BDL-bloodB/")
MB_BDL_B <- CreateSeuratObject(counts = MB_BDL_B.data, project = "BDL-blood_B")
saveRDS(MB_BDL_B, file = "./100720/MBBDL_OG_B.RDS")

BDL_1 <- merge(MliHealthy_A, y = c(MliHealthy_B, MBhealthy_A, MBhealthy_B, MliBDL_A, MliBDL_B, MB_BDL_A, MB_BDL_B), add.cell.ids = c("H-liver-A","H-liver-B", "H-blood-A", "H-blood-B", "BDL-liver-A", "BDL-liver-B", "BDL-blood-A", "BDL-blood-B"), project = "BDL_1") 
unique(sapply(X = strsplit(colnames(BDL_1), split = "_"), FUN = "[", 1))

# set %mt
BDL_1[["percent.mt"]] <- PercentageFeatureSet(BDL_1, pattern = "^mt-")

# QC
plot_QC_Vlnplot <- VlnPlot(BDL_1, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "./100720/BDL1_QCviolin.png", plot = plot_QC_Vlnplot, width = 25, units = "cm")

plot_QC_Scatter1 <- FeatureScatter(BDL_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_QC_Scatter2 <- FeatureScatter(BDL_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2
ggsave(filename = "./100720/BDL1_QCscatter.png", plot = plot_QC_Scatter, width = 25, units = "cm")

                                                                            


#Filter cells
BDL_1 <- subset(BDL_1, subset = nFeature_RNA > 417 & percent.mt < 20)

#Normalise
BDL_1 <- NormalizeData(BDL_1, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection

BDL_1 <- FindVariableFeatures(BDL_1, selection.method = "vst", nfeatures = 2000)
BDL_1_top15 <- head(VariableFeatures(BDL_1), 15)
Plot_FS_15 <- LabelPoints(plot = VariableFeaturePlot(BDL_1), points = BDL_1_top15, repel = TRUE)
ggsave(filename = "./100720/BDL1_FS15.png", plot = Plot_FS_15, width = 30, units = "cm")


#Scale 
BDL_1 <- ScaleData(BDL_1)


#PCA
BDL_1 <- RunPCA(BDL_1, features = VariableFeatures(object = BDL_1))
capture.output(print(BDL_1[["pca"]], dims = 1:10, nfeatures = 5), file = "./100720/BDL1_PCA10.txt")

saveRDS(BDL_1, file = "./100720/BDL1.RDS")

#n(PC) Plots
jpeg(file="./100720/BDL1_Heatmap18.jpeg", width = 480, height = 480, quality = 90)
DimHeatmap(BDL_1, dims = 1:18, cells = 500, balanced = TRUE)
dev.off()

jpeg(file="./100720/BDL1_Heatmap36.jpeg", width = 480, height = 480, quality = 90)
DimHeatmap(BDL_1, dims = 19:36, cells = 500, balanced = TRUE)
dev.off()

jpeg(file="./100720/BDL1_Heatmap50.jpeg", width = 480, height = 480, quality = 90)
DimHeatmap(BDL_1, dims = 37:50, cells = 500, balanced = TRUE)
dev.off()

plot_Elbow <- ElbowPlot(BDL_1, ndims = 50)
ggsave(filename = "./100720/BDL1_ElbowPlot.png", plot = plot_Elbow, width = 30, units = "cm")

BDL1_JSP <- JackStraw(BDL_1, dims = 50)
BDL1_JSP <- ScoreJackStraw(BDL1_JSP, dims = 1:50)
plot_Jackstraw <- JackStrawPlot(BDL1_JSP, dims = 1:50)
ggsave(filename = "./100720/BDL1_JSP.png", plot = plot_Jackstraw, width = 30, units = "cm")


# UMAP of Possible n(PC)s -> 24, 34, 47

#24
BDL1_PC24 <- FindNeighbors(BDL_1, dims = 1:24)
BDL1_PC24 <- FindClusters(BDL1_PC24, resolution = 1.0)
BDL1_PC24 <- RunUMAP(BDL1_PC24, dims = 1:24)

Plot_PC24 <- DimPlot(BDL1_PC24, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "./110720/BDL1_PC24_UMAP.jpeg", plot = Plot_PC24, width = 30, units = "cm")

#34
BDL1_PC34 <- FindNeighbors(BDL_1, dims = 1:34)
BDL1_PC34 <- FindClusters(BDL1_PC34, resolution = 1.0)
BDL1_PC34 <- RunUMAP(BDL1_PC34, dims = 1:34)

Plot_PC34 <- DimPlot(BDL1_PC34, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "./110720/BDL1_PC34_UMAP.jpeg", plot = Plot_PC34, width = 30, units = "cm")


#47
BDL1_PC47 <- FindNeighbors(BDL_1, dims = 1:47)
BDL1_PC47 <- FindClusters(BDL1_PC47, resolution = 1.0)
BDL1_PC47 <- RunUMAP(BDL1_PC47, dims = 1:47)

Plot_PC47 <- DimPlot(BDL1_PC47, reduction = "umap", group.by = "orig.ident")
ggsave(filename = "./110720/BDL1_PC47_UMAP.jpeg", plot = Plot_PC47, width = 30, units = "cm")

