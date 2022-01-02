
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
getwd()

# Load data
MLhealthy.data <- Read10X(data.dir = "./Filtered Matrices/H-lung/")
MLhealthy <- CreateSeuratObject(counts = MLhealthy.data, project = "H-lung")
#saveRDS(MLhealthy, file = "./100720/MLhealthyOG.RDS")

MBhealthy.data <- Read10X(data.dir = "./Filtered Matrices/H-blood/")
MBhealthy <- CreateSeuratObject(counts = MBhealthy.data, project = "H-blood")
#saveRDS(MBhealthy, file = "./100720/MBhealthyOG.RDS")

MLbleo.data <- Read10X(data.dir = "./Filtered Matrices/Bleo-lung/")
MLbleo <- CreateSeuratObject(counts = MLbleo.data, project = "Bleo-lung")
#saveRDS(MLbleo, file = "./100720/MLbleoOG.RDS")

MBbleo.data <- Read10X(data.dir = "./Filtered Matrices/Bleo-blood/")
MBbleo <- CreateSeuratObject(counts = MBbleo.data, project = "Bleo-blood")
#saveRDS(MBbleo, file = "./100720/MBbleoOG.RDS")

Bleo1 <- merge(MLhealthy, y = c(MBhealthy, MLbleo, MBbleo), add.cell.ids = c("H-lung", "H-blood", "Bleo-lung", "Bleo-blood"), project = "Bleo1") 
unique(sapply(X = strsplit(colnames(Bleo1), split = "_"), FUN = "[", 1))

# set %mt
Bleo1[["percent.mt"]] <- PercentageFeatureSet(Bleo1, pattern = "^mt-")

# QC
#plot_QC_Vlnplot <- VlnPlot(Bleo1, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#ggsave(filename = "./100720/Bleo1_QCviolin.png", plot = plot_QC_Vlnplot, width = 25, units = "cm")

#plot_QC_Scatter1 <- FeatureScatter(Bleo1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot_QC_Scatter2 <- FeatureScatter(Bleo1, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2
#ggsave(filename = "./100720/Bleo1_QCscatter.png", plot = plot_QC_Scatter, width = 25, units = "cm")




#Filter cells
Bleo1 <- subset(Bleo1, subset = nFeature_RNA > 550 & nFeature_RNA < 5760 & percent.mt < 21)

#Normalise
Bleo1 <-  NormalizeData(Bleo1, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection

Bleo1 <- FindVariableFeatures(Bleo1, selection.method = "vst", nfeatures = 2000)
#Bleo1_top15 <- head(VariableFeatures(Bleo1), 15)
#Plot_FS_15 <- LabelPoints(plot = VariableFeaturePlot(Bleo1), points = Bleo1_top15, repel = TRUE)
#ggsave(filename = "./100720/Bleo1_FS15.png", plot = Plot_FS_15, width = 30, units = "cm")


#Scale 
Bleo1 <- ScaleData(Bleo1)

#PCA
Bleo1 <- RunPCA(Bleo1, features = VariableFeatures(object = Bleo1))
#capture.output(print(Bleo1[["pca"]], dims = 1:10, nfeatures = 5), file = "./100720/Bleo1_PCA10.txt")

#saveRDS(Bleo1, file = "./100720/Bleo1.RDS")

#n(PC) Plots
#jpeg(file="./100720/Bleo1_Heatmap18.jpeg", width = 480, height = 480, quality = 100)
#DimHeatmap(Bleo1, dims = 1:18, cells = 500, balanced = TRUE)
#dev.off()

#jpeg(file="./100720/Bleo1_Heatmap36.jpeg", width = 480, height = 480, quality = 90)
#DimHeatmap(Bleo1, dims = 19:36, cells = 500, balanced = TRUE)
#dev.off()

#jpeg(file="./100720/Bleo1_Heatmap50.jpeg", width = 480, height = 480, quality = 90)
#DimHeatmap(Bleo1, dims = 37:50, cells = 500, balanced = TRUE)
#dev.off()

#plot_Elbow <- ElbowPlot(Bleo1, ndims = 50)
#ggsave(filename = "./100720/Bleo1_ElbowPlot.png", plot = plot_Elbow, width = 30, units = "cm")

#Bleo1_JSP <- JackStraw(Bleo1, dims = 50)
#Bleo1_JSP <- ScoreJackStraw(Bleo1_JSP, dims = 1:50)
#plot_Jackstraw <- JackStrawPlot(Bleo1_JSP, dims = 1:50)
#ggsave(filename = "./100720/Bleo1_JSP.png", plot = plot_Jackstraw, width = 30, units = "cm")


# UMAP of Possible n(PC)s -> 27, 35, 39, 45

#27
#Bleo1_PC27 <- FindNeighbors(Bleo1, dims = 1:27)
#Bleo1_PC27 <- FindClusters(Bleo1_PC27, resolution = 1.0)
#Bleo1_PC27 <- RunUMAP(Bleo1_PC27, dims = 1:27)

#Plot_PC27 <- DimPlot(Bleo1_PC27, reduction = "umap", group.by = "orig.ident")
#ggsave(filename = "./110720/Bleo1_PC27_UMAP.jpeg", plot = Plot_PC27, width = 30, units = "cm")

#35
#Bleo1_PC35 <- FindNeighbors(Bleo1, dims = 1:35)
#Bleo1_PC35 <- FindClusters(Bleo1_PC35, resolution = 1.0)
#Bleo1_PC35 <- RunUMAP(Bleo1_PC35, dims = 1:35)

# Plot_PC35 <- DimPlot(Bleo1_PC35, reduction = "umap", group.by = "orig.ident")
# ggsave(filename = "./110720/Bleo1_PC35_UMAP.jpeg", plot = Plot_PC35, width = 30, units = "cm")
# 
# #39
# Bleo1_PC39 <- FindNeighbors(Bleo1, dims = 1:39)
# Bleo1_PC39 <- FindClusters(Bleo1_PC39, resolution = 1.0)
# Bleo1_PC39 <- RunUMAP(Bleo1_PC39, dims = 1:39)
# 
# Plot_PC39 <- DimPlot(Bleo1_PC39, reduction = "umap", group.by = "orig.ident")
# ggsave(filename = "./110720/Bleo1_PC39_UMAP.jpeg", plot = Plot_PC39, width = 30, units = "cm")
# 
# #45
Bleo1_PC45 <- FindNeighbors(Bleo1, dims = 1:45)
Bleo1_PC45 <- FindClusters(Bleo1_PC45, resolution = 1.0)
Bleo1_PC45 <- RunUMAP(Bleo1_PC45, dims = 1:45)
# 
# Plot_PC45 <- DimPlot(Bleo1_PC45, reduction = "umap", group.by = "orig.ident")
# ggsave(filename = "./110720/Bleo1_PC45_UMAP.jpeg", plot = Plot_PC45, width = 30, units = "cm")



# Bleo1 <- readRDS("../Bleo1_Graphs/Bleo1.RDS")
saveRDS(Bleo1_PC45, "../Bleo1_Graphs/Bleo1.RDS")

# Bleo1 <- ScaleData(Bleo1)
# Bleo1 <- RunPCA(Bleo1, features = VariableFeatures(object = Bleo1))
# Bleo1_PC45 <- FindNeighbors(Bleo1, dims = 1:45)
# Bleo1_PC45 <- FindClusters(Bleo1_PC45, resolution = 1.0)
# Bleo1_PC45 <- RunUMAP(Bleo1_PC45, dims = 1:45)

Bleo1 <- Bleo1_PC45
# Plot big markers

plot_MonoMacs <- FeaturePlot(Bleo1, features = c("Csf1r", "Lyz2", "Itgam", "Cd68", "S100a8", "S100a4"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_MonoMacs.png", plot = plot_MonoMacs, width = 30, height = 20, units = "cm")

plot_MonoLy6c <- FeaturePlot(Bleo1, features = c("Ly6c2", "Ccr2", "Vcan", "Ace", "TreMK4"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_MonoLy6c.png", plot = plot_MonoLy6c, width = 30, height = 20, units = "cm")

plot_NeuTNKB <- FeaturePlot(Bleo1, features = c("Ly6g", "Mpo", "Cd3d", "Cd3e", "Cd8a", "Cd4", "Klra9", "Klrb1b", "Gzma", "Cd79a", "Cd79b"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_NeuTNKB.png", plot = plot_NeuTNKB, width = 30, height = 40, units = "cm")

plot_bloods_cDC <- FeaturePlot(Bleo1, features = c("Fcer1a", "Il4", "Cpa3", "Kdr", "Pecam1", "Icam2", "Cdh5", "Hba-a1", "Hbb-a2", "Xcr1", "Cd209a", "Siglech"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_bloodsCDC.png", plot = plot_bloods_cDC, width = 30, height = 40, units = "cm")

plot_MonoLy6cHi <- FeaturePlot(Bleo1, features = c("Ly6c2", "Ccr2", "Vcan"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_MonoLy6cHIGH.png", plot = plot_MonoLy6cHi, width = 30, height = 10, units = "cm")

plot_MonoLy6cLo <- FeaturePlot(Bleo1, features = c("Ace", "TreMK4"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_MonoLy6cLOW.png", plot = plot_MonoLy6cLo, width = 20, height = 10, units = "cm")

plot_Neu <- FeaturePlot(Bleo1, features = c("Ly6g", "Mpo"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_Neu.png", plot = plot_Neu, width = 20, height = 10, units = "cm")

plot_T <- FeaturePlot(Bleo1, features = c("Cd3d", "Cd3e", "Cd8a", "Cd4"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_T.png", plot = plot_T, width = 20, height = 20, units = "cm")

plot_NK <- FeaturePlot(Bleo1, features = c( "Klra9", "Klrb1b", "Gzma"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_NK.png", plot = plot_NK, width = 30, height = 10, units = "cm")

plot_B <- FeaturePlot(Bleo1, features = c("Cd79a", "Cd79b"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_B.png", plot = plot_B, width = 20, height = 10, units = "cm")

plot_DC <- FeaturePlot(Bleo1, features = c("Xcr1", "Cd209a", "Siglech"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_DC.png", plot = plot_DC, width = 30, height = 10, units = "cm")

plot_baso <- FeaturePlot(Bleo1, features = c("Fcer1a", "Il4", "Cpa3"), ncol = 3)
ggsave(filename = "./170720_Bleo1_FP_baso.png", plot = plot_baso, width = 30, height = 10, units = "cm")

plot_EC <- FeaturePlot(Bleo1, features = c("Kdr", "Pecam1", "Icam2", "Cdh5"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_EC.png", plot = plot_EC, width = 20, height = 20, units = "cm")

plot_HBC <- FeaturePlot(Bleo1, features = c("Hba-a1", "Hbb-a2"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_HBC.png", plot = plot_HBC, width = 20, height = 10, units = "cm")

plot_Macro <- FeaturePlot(Bleo1, features = c("Mmp12", "Ccl12", "Mgl2", "Pf4"), ncol = 2)
ggsave(filename = "./170720_Bleo1_FP_Macro.png", plot = plot_Macro, width = 20, height = 20, units = "cm")

