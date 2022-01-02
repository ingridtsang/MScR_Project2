library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Pre-processing
#mt gene 

MLhealthy[["percent.mt"]] <- PercentageFeatureSet(MLhealthy, pattern = "^mt-")
MBhealthy[["percent.mt"]] <- PercentageFeatureSet(MBhealthy, pattern = "^mt-")
MLbleo[["percent.mt"]] <- PercentageFeatureSet(MLbleo, pattern = "^mt-")
MBbleo[["percent.mt"]] <- PercentageFeatureSet(MBbleo, pattern = "^mt-")

#QC Filter

MLhealthy <- subset(MLhealthy, subset = nFeature_RNA > 550 & nFeature_RNA < 5760 & percent.mt < 21)
MBhealthy <- subset(MBhealthy, subset = nFeature_RNA > 550 & nFeature_RNA < 5760 & percent.mt < 21)
MLbleo <- subset(MLbleo, subset = nFeature_RNA > 550 & nFeature_RNA < 5760 & percent.mt < 21)
MBbleo <- subset(MBbleo, subset = nFeature_RNA > 550 & nFeature_RNA < 5760 & percent.mt < 21)

#LogNorm, FS, Scale
MLhealthy <- NormalizeData(MLhealthy, normalization.method = "LogNormalize", scale.factor = 10000)
MLhealthy <- FindVariableFeatures(MLhealthy, selection.method = "vst", nfeatures = 2000)
MLhealthy <- ScaleData(MLhealthy)

MBhealthy <- NormalizeData(MBhealthy, normalization.method = "LogNormalize", scale.factor = 10000)
MBhealthy <- FindVariableFeatures(MBhealthy, selection.method = "vst", nfeatures = 2000)
MBhealthy <- ScaleData(MBhealthy)

MLbleo <- NormalizeData(MLbleo, normalization.method = "LogNormalize", scale.factor = 10000)
MLbleo <- FindVariableFeatures(MLbleo, selection.method = "vst", nfeatures = 2000)
MLbleo <- ScaleData(MLbleo)

MBbleo <- NormalizeData(MBbleo, normalization.method = "LogNormalize", scale.factor = 10000)
MBbleo <- FindVariableFeatures(MBbleo, selection.method = "vst", nfeatures = 2000)
MBbleo <- ScaleData(MBbleo)


#PCA
MLhealthy_PCA <- RunPCA(MLhealthy, features = VariableFeatures(object = MLhealthy))
MBhealthy_PCA <- RunPCA(MBhealthy, features = VariableFeatures(object = MBhealthy))

MLbleo_PCA <- RunPCA(MLbleo, features = VariableFeatures(object = MLbleo))
MBbleo_PCA <- RunPCA(MBbleo, features = VariableFeatures(object = MBbleo))


#UMAP
MLhealthy_PCA <- FindNeighbors(MLhealthy_PCA, dims = 1:45)
MLhealthy_PCA <- FindClusters(MLhealthy_PCA, resolution = 1.0)
MLhealthy_UMAP <- RunUMAP(MLhealthy_PCA, dims = 1:45)
saveRDS(MLhealthy_UMAP, file = "./110720/MLhealthy_UMAP.RDS")

MLhealthy_nExp_poi <- round(0.05*length(MLhealthy_UMAP$orig.ident)) 
print("MLhealthy")
print(MLhealthy_nExp_poi)

MBhealthy_PCA <- FindNeighbors(MBhealthy_PCA, dims = 1:45)
MBhealthy_PCA <- FindClusters(MBhealthy_PCA, resolution = 1.0)
MBhealthy_UMAP <- RunUMAP(MBhealthy_PCA, dims = 1:45)
saveRDS(MBhealthy_UMAP, file = "./110720/MBhealthy_UMAP.RDS")

MBhealthy_nExp_poi <- round(0.05*length(MBhealthy_UMAP$orig.ident)) 
print("MBhealthy")
print(MBhealthy_nExp_poi)

MLbleo_PCA <- FindNeighbors(MLbleo_PCA, dims = 1:45)
MLbleo_PCA <- FindClusters(MLbleo_PCA, resolution = 1.0)
MLbleo_UMAP <- RunUMAP(MLbleo_PCA, dims = 1:45)
saveRDS(MLbleo_UMAP, file = "./110720/MLbleo_UMAP.RDS")

MLbleo_nExp_poi <- round(0.05*length(MLbleo_UMAP$orig.ident)) 
print("MLbleo")
print(MLbleo_nExp_poi)

MBbleo_PCA <- FindNeighbors(MBbleo_PCA, dims = 1:45)
MBbleo_PCA <- FindClusters(MBbleo_PCA, resolution = 1.0)
MBbleo_UMAP <- RunUMAP(MBbleo_PCA, dims = 1:45)
saveRDS(MBbleo_UMAP, file = "./110720/MBbleo_UMAP.RDS")

MBbleo_nExp_poi <- round(0.05*length(MBbleo_UMAP$orig.ident)) 
print("MBbleo")
print(MBbleo_nExp_poi)

# pK Identification (no ground-truth)
#select pc according to elbow plot
sweep.res.list_MLhealthy <- paramSweep_v3(MLhealthy_UMAP, PCs = 1:45, sct = FALSE)
sweep.stats_MLhealthy <- summarizeSweep(sweep.res.list_MLhealthy, GT = FALSE)
bcmvn_MLhealthy <- find.pK(sweep.stats_MLhealthy)

sweep.res.list_MBhealthy <- paramSweep_v3(MBhealthy_UMAP, PCs = 1:45, sct = FALSE)
sweep.stats_MBhealthy <- summarizeSweep(sweep.res.list_MBhealthy, GT = FALSE)
bcmvn_MBhealthy <- find.pK(sweep.stats_MBhealthy)

sweep.res.list_MLbleo <- paramSweep_v3(MLbleo_UMAP, PCs = 1:45, sct = FALSE)
sweep.stats_MLbleo <- summarizeSweep(sweep.res.list_MLbleo, GT = FALSE)
bcmvn_MLbleo <- find.pK(sweep.stats_MLbleo)

sweep.res.list_MBbleo <- paramSweep_v3(MBbleo_UMAP, PCs = 1:45, sct = FALSE)
sweep.stats_MBbleo <- summarizeSweep(sweep.res.list_MBbleo, GT = FALSE)
bcmvn_MBbleo <- find.pK(sweep.stats_MBbleo)

# set pk according to the pk plot
MLhealthy_nExp_poi <- round(0.05*length(MLhealthy_UMAP$orig.ident))   ## Assuming 5% doublet formation rate(set your rate according to 10X manual)
MBhealthy_nExp_poi <- round(0.05*length(MBhealthy_UMAP$orig.ident))
MLbleo_nExp_poi <- round(0.05*length(MLbleo_UMAP$orig.ident))
MBbleo_nExp_poi <- round(0.05*length(MBbleo_UMAP$orig.ident))

MLhealthy_rmD <- doubletFinder_v3(MLhealthy_UMAP, PCs = 1:45, pN = 0.25, pK = 0.06, nExp = MLhealthy_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MLhealthy_rmD, file = "./110720/Bleo1_MLhealthy_rmD.RDS")

MBhealthy_rmD <- doubletFinder_v3(MBhealthy_UMAP, PCs = 1:45, pN = 0.25, pK = 0.06, nExp = MBhealthy_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MBhealthy_rmD, file = "./110720/Bleo1_MBhealthy_rmD.RDS")

MLbleo_rmD <- doubletFinder_v3(MLbleo_UMAP, PCs = 1:45, pN = 0.25, pK = 0.06, nExp = MLbleo_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MLbleo_rmD, file = "./110720/Bleo1_MLbleo_rmD.RDS")

MBbleo_rmD <- doubletFinder_v3(MBbleo_UMAP, PCs = 1:45, pN = 0.25, pK = 0.06, nExp = MBbleo_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
saveRDS(MBbleo_rmD, file = "./110720/Bleo1_MBbleo_rmD.RDS")

