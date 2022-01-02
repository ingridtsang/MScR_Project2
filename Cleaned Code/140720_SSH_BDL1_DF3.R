# Based on 110720_SSH_BDL_DF.R -> but fixed saveRDS in UMAP and added indie UMAPS to save just to follow the coding
library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Load data
MliHealthy_A.data <- Read10X(data.dir = "./Filtered Matrices/H-liverA/")
MliHealthy_A <- CreateSeuratObject(counts = MliHealthy_A.data, project = "H-liver_A")

MliHealthy_B.data <- Read10X(data.dir = "./Filtered Matrices/H-liverB/")
MliHealthy_B <- CreateSeuratObject(counts = MliHealthy_B.data, project = "H-liver_B")

MBhealthy_A.data <- Read10X(data.dir = "./Filtered Matrices/H-bloodA/")
MBhealthy_A <- CreateSeuratObject(counts = MBhealthy_A.data, project = "H-blood_A")

MBhealthy_B.data <- Read10X(data.dir = "./Filtered Matrices/H-bloodB/")
MBhealthy_B <- CreateSeuratObject(counts = MBhealthy_B.data, project = "H-blood_B")

MliBDL_A.data <- Read10X(data.dir = "./Filtered Matrices/BDL-liverA/")
MliBDL_A <- CreateSeuratObject(counts = MliBDL_A.data, project = "BDL-liver_A")

MliBDL_B.data <- Read10X(data.dir = "./Filtered Matrices/BDL-liverB/")
MliBDL_B <- CreateSeuratObject(counts = MliBDL_B.data, project = "BDL-liver_B")

MB_BDL_A.data <- Read10X(data.dir = "./Filtered Matrices/BDL-bloodA/")
MB_BDL_A <- CreateSeuratObject(counts = MB_BDL_A.data, project = "BDL-blood_A")

MB_BDL_B.data <- Read10X(data.dir = "./Filtered Matrices/BDL-bloodB/")
MB_BDL_B <- CreateSeuratObject(counts = MB_BDL_B.data, project = "BDL-blood_B")

# Pre-processing
#mt gene 

MliHealthy_A[["percent.mt"]] <- PercentageFeatureSet(MliHealthy_A, pattern = "^mt-")
MliHealthy_B[["percent.mt"]] <- PercentageFeatureSet(MliHealthy_B, pattern = "^mt-")
MBhealthy_A[["percent.mt"]] <- PercentageFeatureSet(MBhealthy_A, pattern = "^mt-")
MBhealthy_B[["percent.mt"]] <- PercentageFeatureSet(MBhealthy_B, pattern = "^mt-")

MliBDL_A[["percent.mt"]] <- PercentageFeatureSet(MliBDL_A, pattern = "^mt-")
MliBDL_B[["percent.mt"]] <- PercentageFeatureSet(MliBDL_B, pattern = "^mt-")
MB_BDL_A[["percent.mt"]] <- PercentageFeatureSet(MB_BDL_A, pattern = "^mt-")
MB_BDL_B[["percent.mt"]] <- PercentageFeatureSet(MB_BDL_B, pattern = "^mt-")


#QC Filter

MliHealthy_A <- subset(MliHealthy_A, subset = nFeature_RNA > 417 & percent.mt < 20)
MliHealthy_B <- subset(MliHealthy_B, subset = nFeature_RNA > 417 & percent.mt < 20)
MBhealthy_A <- subset(MBhealthy_A, subset = nFeature_RNA > 417 & percent.mt < 20)
MBhealthy_B <- subset(MBhealthy_B,subset = nFeature_RNA > 417 & percent.mt < 20)

MliBDL_A <- subset(MliBDL_A, subset = nFeature_RNA > 417 & percent.mt < 20)
MliBDL_B <- subset(MliBDL_B, subset = nFeature_RNA > 417 & percent.mt < 20)
MB_BDL_A <- subset(MB_BDL_A, subset = nFeature_RNA > 417 & percent.mt < 20)
MB_BDL_B <- subset(MB_BDL_B, subset = nFeature_RNA > 417 & percent.mt < 20)

#LogNorm, FS, Scale
MliHealthy_A <- NormalizeData(MliHealthy_A, normalization.method = "LogNormalize", scale.factor = 10000)
MliHealthy_A <- FindVariableFeatures(MliHealthy_A, selection.method = "vst", nfeatures = 2000)
MliHealthy_A <- ScaleData(MliHealthy_A)

MliHealthy_B <- NormalizeData(MliHealthy_B, normalization.method = "LogNormalize", scale.factor = 10000)
MliHealthy_B <- FindVariableFeatures(MliHealthy_B, selection.method = "vst", nfeatures = 2000)
MliHealthy_B <- ScaleData(MliHealthy_B)

MBhealthy_A <- NormalizeData(MBhealthy_A, normalization.method = "LogNormalize", scale.factor = 10000)
MBhealthy_A <- FindVariableFeatures(MBhealthy_A, selection.method = "vst", nfeatures = 2000)
MBhealthy_A <- ScaleData(MBhealthy_A)

MBhealthy_B <- NormalizeData(MBhealthy_B, normalization.method = "LogNormalize", scale.factor = 10000)
MBhealthy_B <- FindVariableFeatures(MBhealthy_B, selection.method = "vst", nfeatures = 2000)
MBhealthy_B <- ScaleData(MBhealthy_B)


MliBDL_A <- NormalizeData(MliBDL_A, normalization.method = "LogNormalize", scale.factor = 10000)
MliBDL_A <- FindVariableFeatures(MliBDL_A, selection.method = "vst", nfeatures = 2000)
MliBDL_A <- ScaleData(MliBDL_A)

MliBDL_B <- NormalizeData(MliBDL_B, normalization.method = "LogNormalize", scale.factor = 10000)
MliBDL_B <- FindVariableFeatures(MliBDL_B, selection.method = "vst", nfeatures = 2000)
MliBDL_B <- ScaleData(MliBDL_B)


MB_BDL_A <- NormalizeData(MB_BDL_A, normalization.method = "LogNormalize", scale.factor = 10000)
MB_BDL_A <- FindVariableFeatures(MB_BDL_A, selection.method = "vst", nfeatures = 2000)
MB_BDL_A <- ScaleData(MB_BDL_A)

MB_BDL_B <- NormalizeData(MB_BDL_B, normalization.method = "LogNormalize", scale.factor = 10000)
MB_BDL_B <- FindVariableFeatures(MB_BDL_B, selection.method = "vst", nfeatures = 2000)
MB_BDL_B <- ScaleData(MB_BDL_B)

#PCA
MliHealthy_A_PCA <- RunPCA(MliHealthy_A, features = VariableFeatures(object = MliHealthy_A))
MliHealthy_B_PCA <- RunPCA(MliHealthy_B, features = VariableFeatures(object = MliHealthy_B))

MBhealthy_A_PCA <- RunPCA(MBhealthy_A, features = VariableFeatures(object = MBhealthy_A))
MBhealthy_B_PCA <- RunPCA(MBhealthy_B, features = VariableFeatures(object = MBhealthy_B))

MliBDL_A_PCA <- RunPCA(MliBDL_A, features = VariableFeatures(object = MliBDL_A))
MliBDL_B_PCA <- RunPCA(MliBDL_B, features = VariableFeatures(object = MliBDL_B))

MB_BDL_A_PCA <- RunPCA(MB_BDL_A, features = VariableFeatures(object = MB_BDL_A))
MB_BDL_B_PCA <- RunPCA(MB_BDL_B, features = VariableFeatures(object = MB_BDL_B))


#UMAP, removed space after [object] in saveRDS; make UMAPs to track progress of code
MliHealthy_A_PCA <- FindNeighbors(MliHealthy_A_PCA, dims = 1:47)
MliHealthy_A_PCA <- FindClusters(MliHealthy_A_PCA, resolution = 1.0)
MliHealthy_A_UMAP <- RunUMAP(MliHealthy_A_PCA, dims = 1:47)
saveRDS(MliHealthy_A_UMAP, file = "./110720/BDL1_MLiH_A_UMAP.rds")

Plot_MliH_A_UMAP <- DimPlot(MliHealthy_A_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiH_A_UMAP.jpeg", plot = Plot_MliH_A_UMAP, width = 30, units = "cm")


MliHealthy_B_PCA <- FindNeighbors(MliHealthy_B_PCA, dims = 1:47)
MliHealthy_B_PCA <- FindClusters(MliHealthy_B_PCA, resolution = 1.0)
MliHealthy_B_UMAP <- RunUMAP(MliHealthy_B_PCA, dims = 1:47)
saveRDS(MliHealthy_B_UMAP, file = "./140720/BDL1_MLiH_B_UMAP.rds")

Plot_MliH_B_UMAP <- DimPlot(MliHealthy_B_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiH_B_UMAP.jpeg", plot = Plot_MliH_B_UMAP, width = 30, units = "cm")

MBhealthy_A_PCA <- FindNeighbors(MBhealthy_A_PCA, dims = 1:47)
MBhealthy_A_PCA <- FindClusters(MBhealthy_A_PCA, resolution = 1.0)
MBhealthy_A_UMAP <- RunUMAP(MBhealthy_A_PCA, dims = 1:47)
saveRDS(MBhealthy_A_UMAP, file = "./140720/BDL1_MBHealthy_A_UMAP.rds")

Plot_MBH_A_UMAP <- DimPlot(MBhealthy_A_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBhealthy_A_UMAP.jpeg", plot = Plot_MBH_A_UMAP, width = 30, units = "cm")

MBhealthy_B_PCA <- FindNeighbors(MBhealthy_B_PCA, dims = 1:47)
MBhealthy_B_PCA <- FindClusters(MBhealthy_B_PCA, resolution = 1.0)
MBhealthy_B_UMAP <- RunUMAP(MBhealthy_B_PCA, dims = 1:47)
saveRDS(MBhealthy_B_UMAP, file = "./140720/BDL1_MBHealthy_B_UMAP.rds")

Plot_MBH_B_UMAP <- DimPlot(MBhealthy_B_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBhealthy_B_UMAP.jpeg", plot = Plot_MBH_B_UMAP, width = 30, units = "cm")

MliBDL_A_PCA <- FindNeighbors(MliBDL_A_PCA, dims = 1:47)
MliBDL_A_PCA <- FindClusters(MliBDL_A_PCA, resolution = 1.0)
MliBDL_A_UMAP <- RunUMAP(MliBDL_A_PCA, dims = 1:47)
saveRDS(MliBDL_A_UMAP, file = "./140720/BDL1_MLiBDL_A_UMAP.rds")

Plot_MliBDL_A_UMAP <- DimPlot(MliBDL_A_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiBDL_A_UMAP.jpeg", plot = Plot_MliBDL_A_UMAP, width = 30, units = "cm")

MliBDL_B_PCA <- FindNeighbors(MliBDL_B_PCA, dims = 1:47)
MliBDL_B_PCA <- FindClusters(MliBDL_B_PCA, resolution = 1.0)
MliBDL_B_UMAP <- RunUMAP(MliBDL_B_PCA, dims = 1:47)
saveRDS(MliBDL_B_UMAP, file = "./140720/BDL1_MLiBDL_B_UMAP.rds")

Plot_MliBDL_B_UMAP <- DimPlot(MliBDL_B_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiBDL_B_UMAP.jpeg", plot = Plot_MliBDL_B_UMAP, width = 30, units = "cm")

MB_BDL_A_PCA <- FindNeighbors(MB_BDL_A_PCA, dims = 1:47)
MB_BDL_A_PCA <- FindClusters(MB_BDL_A_PCA, resolution = 1.0)
MB_BDL_A_UMAP <- RunUMAP(MB_BDL_A_PCA, dims = 1:47)
saveRDS(MB_BDL_A_UMAP, file = "./140720/BDL1_MBBDL_A_UMAP.rds")

Plot_MBBDL_A_UMAP <- DimPlot(MB_BDL_A_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBBDL_A_UMAP.jpeg", plot = Plot_MBBDL_A_UMAP, width = 30, units = "cm")

MB_BDL_B_PCA <- FindNeighbors(MB_BDL_B_PCA, dims = 1:47)
MB_BDL_B_PCA <- FindClusters(MB_BDL_B_PCA, resolution = 1.0)
MB_BDL_B_UMAP <- RunUMAP(MB_BDL_B_PCA, dims = 1:47)
saveRDS(MB_BDL_B_UMAP, file = "./140720/BDL1_MBBDL_B_UMAP.rds")

Plot_MBBDL_A_UMAP <- DimPlot(MB_BDL_A_UMAP, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBBDL_A_UMAP.jpeg", plot = Plot_MBBDL_A_UMAP, width = 30, units = "cm")


# nExp_POI; Assuming 5% doublet formation rate(set your rate according to 10X manual)
MliHealthy_A_nExp_poi <- round(0.05*length(MliHealthy_A_UMAP$orig.ident)) 
MliHealthy_B_nExp_poi <- round(0.05*length(MliHealthy_B_UMAP$orig.ident)) 
MBhealthy_A_nExp_poi <- round(0.05*length(MBhealthy_A_UMAP$orig.ident)) 
MBhealthy_B_nExp_poi <- round(0.05*length(MBhealthy_B_UMAP$orig.ident)) 

MliBDL_A_nExp_poi <- round(0.05*length(MliBDL_A_UMAP$orig.ident)) 
MliBDL_B_nExp_poi <- round(0.05*length(MliBDL_B_UMAP$orig.ident)) 
MB_BDL_A_nExp_poi <- round(0.05*length(MB_BDL_A_UMAP$orig.ident)) 
MB_BDL_B_nExp_poi <- round(0.05*length(MB_BDL_B_UMAP$orig.ident)) 




# pK Identification (no ground-truth)
#select pc according to elbow plot
sweep.res.list_MliHealthy_A_UMAP <- paramSweep_v3(MliHealthy_A_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MliHealthy_A_UMAP <- summarizeSweep(sweep.res.list_MliHealthy_A_UMAP, GT = FALSE)
bcmvn_MliHealthy_A_UMAP <- find.pK(sweep.stats_MliHealthy_A_UMAP)

sweep.res.list_MliHealthy_B_UMAP <- paramSweep_v3(MliHealthy_B_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MliHealthy_B_UMAP <- summarizeSweep(sweep.res.list_MliHealthy_B_UMAP, GT = FALSE)
bcmvn_MliHealthy_B_UMAP <- find.pK(sweep.stats_MliHealthy_B_UMAP)


sweep.res.list_MBhealthy_A_UMAP <- paramSweep_v3(MBhealthy_A_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MBhealthy_A_UMAP <- summarizeSweep(sweep.res.list_MBhealthy_A_UMAP, GT = FALSE)
bcmvn_MBhealthy_A_UMAP <- find.pK(sweep.stats_MBhealthy_A_UMAP)

sweep.res.list_MBhealthy_B_UMAP <- paramSweep_v3(MBhealthy_B_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MBhealthy_B_UMAP <- summarizeSweep(sweep.res.list_MBhealthy_B_UMAP, GT = FALSE)
bcmvn_MBhealthy_B_UMAP <- find.pK(sweep.stats_MBhealthy_B_UMAP)


sweep.res.list_MliBDL_A_UMAP <- paramSweep_v3(MliBDL_A_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MliBDL_A_UMAP <- summarizeSweep(sweep.res.list_MliBDL_A_UMAP, GT = FALSE)
bcmvn_MliBDL_A_UMAP <- find.pK(sweep.stats_MliBDL_A_UMAP)

sweep.res.list_MliBDL_B_UMAP <- paramSweep_v3(MliBDL_B_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MliBDL_B_UMAP <- summarizeSweep(sweep.res.list_MliBDL_B_UMAP, GT = FALSE)
bcmvn_MliBDL_B_UMAP <- find.pK(sweep.stats_MliBDL_B_UMAP)

sweep.res.list_MB_BDL_A_UMAP <- paramSweep_v3(MB_BDL_A_UMAP, PCs = 1:47, sct = FALSE)
sweep.stats_MB_BDL_A_UMAP <- summarizeSweep(sweep.res.list_MB_BDL_A_UMAP, GT = FALSE)
bcmvn_MB_BDL_A_UMAP <- find.pK(sweep.stats_MB_BDL_A_UMAP)

#  rmD + save objects

MliHealthy_A_rmD <- doubletFinder_v3(MliHealthy_A_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MliHealthy_A_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MliHealthy_A_rmD, file = "./140720/BDL1_MLiH_A_rmD.rds")

MliHealthy_B_rmD <- doubletFinder_v3(MliHealthy_B_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MliHealthy_B_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MliHealthy_B_rmD, file = "./140720/BDL1_MLiH_B_rmD.rds")

MBhealthy_A_rmD <- doubletFinder_v3(MBhealthy_A_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MBhealthy_A_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MBhealthy_A_rmD, file = "./140720/BDL1_MBhealthy_A_rmD.rds")

MBhealthy_B_rmD <- doubletFinder_v3(MBhealthy_B_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MBhealthy_B_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MBhealthy_B_rmD, file = "./140720/BDL1_MBhealthy_B_rmD.rds")

MliBDL_A_rmD <- doubletFinder_v3(MliBDL_A_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MliBDL_A_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MliBDL_A_rmD, file = "./140720/BDL1_MLiBDL_A_rmD.rds")

MliBDL_B_rmD <- doubletFinder_v3(MliBDL_B_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MliBDL_B_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MliBDL_B_rmD, file = "./140720/BDL1_MLiBDL_B_rmD.rds")

MB_BDL_A_rmD <- doubletFinder_v3(MB_BDL_A_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MB_BDL_A_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MB_BDL_A_rmD, file = "./140720/BDL1_MBBDL_A_rmD.rds")

MB_BDL_B_rmD <- doubletFinder_v3(MB_BDL_B_UMAP, PCs = 1:47, pN = 0.25, pK = 0.06, nExp = MB_BDL_B_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MB_BDL_B_rmD, file = "./140720/BDL1_MBBDL_B_rmD.rds")


# DF Plots
plot_MliHealthy_A_DF <- DimPlot(MliHealthy_A_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_124")
ggsave(filename = "./140720/BDL1_MLiH_A_DF_UMAP.jpeg", plot = plot_MliHealthy_A_DF, width = 30, units = "cm")

plot_MliHealthy_B_DF <- DimPlot(MliHealthy_B_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_328")
ggsave(filename = "./140720/BDL1_MLiH_B_DF_UMAP.jpeg", plot = plot_MliHealthy_B_DF, width = 30, units = "cm")


plot_MBhealthy_A_DF <- DimPlot(MBhealthy_A_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_147")
ggsave(filename = "./140720/BDL1_MBhealthy_A_DF_UMAP.jpeg", plot = plot_MBhealthy_A_DF, width = 30, units = "cm")

plot_MBhealthy_B_DF <- DimPlot(MBhealthy_B_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_429")
ggsave(filename = "./140720/BDL1_MBhealthy_B_DF_UMAP.jpeg", plot = plot_MBhealthy_B_DF, width = 30, units = "cm")


plot_MliBDL_A_DF <- DimPlot(MliBDL_A_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_121")
ggsave(filename = "./140720/BDL1_MLiBDL_A_DF_UMAP.jpeg", plot = plot_MliBDL_A_DF, width = 30, units = "cm")

plot_MliBDL_B_DF <- DimPlot(MliBDL_B_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_144")
ggsave(filename = "./140720/BDL1_MLiBDL_B_DF_UMAP.jpeg", plot = plot_MliBDL_B_DF, width = 30, units = "cm")


plot_MB_BDL_A_DF <- DimPlot(MB_BDL_A_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_63")
ggsave(filename = "./140720/BDL1_MB_BDL_A_DF_UMAP.jpeg", plot = plot_MB_BDL_A_DF, width = 30, units = "cm")

plot_MB_BDL_B_DF <- DimPlot(MB_BDL_B_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_219")
ggsave(filename = "./140720/BDL1_MB_BDL_B_DF_UMAP.jpeg", plot = plot_MB_BDL_B_DF, width = 30, units = "cm")

# filter doublets
MliHealthy_A_rmD <- subset(MliHealthy_A_rmD, subset = DF.classifications_0.25_0.06_124 =="Singlet")
MliHealthy_B_rmD <- subset(MliHealthy_B_rmD, subset = DF.classifications_0.25_0.06_328 =="Singlet")

MBhealthy_A_rmD <- subset(MBhealthy_A_rmD, subset = DF.classifications_0.25_0.06_147 =="Singlet")
MBhealthy_B_rmD <- subset(MBhealthy_B_rmD, subset = DF.classifications_0.25_0.06_429 =="Singlet")

MliBDL_A_rmD <- subset(MliBDL_A_rmD, subset = DF.classifications_0.25_0.06_121 =="Singlet")
MliBDL_B_rmD <- subset(MliBDL_B_rmD, subset = DF.classifications_0.25_0.06_144 =="Singlet")

MB_BDL_A_rmD <- subset(MB_BDL_A_rmD, subset = DF.classifications_0.25_0.06_63 =="Singlet")
MB_BDL_B_rmD <- subset(MB_BDL_B_rmD, subset = DF.classifications_0.25_0.06_219 =="Singlet")

# rmD Plots

plot_MliHealthy_A_rmD_UMAP <- DimPlot(MliHealthy_A_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiH_A_rmD_UMAP.jpeg", plot = plot_MliHealthy_A_rmD_UMAP, width = 30, units = "cm")

plot_MliHealthy_B_rmD_UMAP <- DimPlot(MliHealthy_B_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiH_B_rmD_UMAP.jpeg", plot = plot_MliHealthy_B_rmD_UMAP, width = 30, units = "cm")


plot_MBhealthy_A_rmD_UMAP <- DimPlot(MBhealthy_A_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBhealthy_A_rmD_UMAP.jpeg", plot = plot_MBhealthy_A_rmD_UMAP, width = 30, units = "cm")

plot_MBhealthy_B_rmD_UMAP <- DimPlot(MBhealthy_B_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBhealthy_B_rmD_UMAP.jpeg", plot = plot_MBhealthy_B_rmD_UMAP, width = 30, units = "cm")


plot_MliBDL_A_rmD_UMAP <- DimPlot(MliBDL_A_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiBDL_A_rmD_UMAP.jpeg", plot = plot_MliBDL_A_rmD_UMAP, width = 30, units = "cm")

plot_MliBDL_B_rmD_UMAP <- DimPlot(MliBDL_B_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MLiBDL_B_rmD_UMAP.jpeg", plot = plot_MliBDL_B_rmD_UMAP, width = 30, units = "cm")

plot_MB_BDL_A_rmD_UMAP <- DimPlot(MB_BDL_A_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBBDL_A_rmD_UMAP.jpeg", plot = plot_MB_BDL_A_rmD_UMAP, width = 30, units = "cm")

plot_MB_BDL_B_rmD_UMAP <- DimPlot(MB_BDL_B_rmD, reduction = "umap")
ggsave(filename = "./140720/BDL1_MBBDL_B_rmD_UMAP.jpeg", plot = plot_MB_BDL_B_rmD_UMAP, width = 30, units = "cm")


# Print nEXP_POI
print("MliHealthyA")
print(MliHealthy_A_nExp_poi)
print("MliHealthyB")
print(MliHealthy_B_nExp_poi)
print(" MBhealthyA")
print(MBhealthy_A_nExp_poi)
print(" MBhealthyB")
print(MBhealthy_B_nExp_poi)
print(" MliBDL_A")
print(MliBDL_A_nExp_poi)
print(" MliBDL_B")
print(MliBDL_B_nExp_poi)
print(" MB_BDL_A")
print(MB_BDL_A_nExp_poi)
print(" MB_BDL_B")
print(MB_BDL_B_nExp_poi)