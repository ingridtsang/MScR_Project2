library(DoubletFinder)
library(umap)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Doublet removal on single datasets

MKhealthy.data <- Read10X(data.dir = "./Filtered Matrices/H-kidney/")
MKhealthy <- CreateSeuratObject(counts = MKhealthy.data, project = "H-kidney")

MBhealthy.data <- Read10X(data.dir = "./Filtered Matrices/H-blood/")
MBhealthy <- CreateSeuratObject(counts = MBhealthy.data, project = "H-blood")

MKinjured.data <- Read10X(data.dir = "./Filtered Matrices/Inj-kidney/")
MKinjured <- CreateSeuratObject(counts = MKinjured.data, project = "Inj-kidney")

MBinjured.data <- Read10X(data.dir = "./Filtered Matrices/Inj-blood/")
MBinjured <- CreateSeuratObject(counts = MBinjured.data, project = "Inj-blood")


# Pre-processing
# mt gene 
MKhealthy[["percent.mt"]] <- PercentageFeatureSet(MKhealthy, pattern = "^mt-")
MBhealthy[["percent.mt"]] <- PercentageFeatureSet(MBhealthy, pattern = "^mt-")

MKinjured[["percent.mt"]] <- PercentageFeatureSet(MKinjured, pattern = "^mt-")
MBinjured[["percent.mt"]] <- PercentageFeatureSet(MBinjured, pattern = "^mt-")

# QC

MKhealthy <- subset(MKhealthy, subset = nFeature_RNA > 450 & nFeature_RNA < 6000 & percent.mt < 15)
MBhealthy <- subset(MBhealthy, subset = nFeature_RNA > 450 & nFeature_RNA < 6000 & percent.mt < 15)

MKinjured <- subset(MKinjured, subset = nFeature_RNA > 450 & nFeature_RNA < 6000 & percent.mt < 15)
MBinjured <- subset(MBinjured, subset = nFeature_RNA > 450 & nFeature_RNA < 6000 & percent.mt < 15)

# LogNorm, FS, Scale
MKhealthy <- NormalizeData(MKhealthy, normalization.method = "LogNormalize", scale.factor = 10000)
MKhealthy <- FindVariableFeatures(MKhealthy, selection.method = "vst", nfeatures = 2000)
MKhealthy <- ScaleData(MKhealthy)

MBhealthy <- NormalizeData(MBhealthy, normalization.method = "LogNormalize", scale.factor = 10000)
MBhealthy <- FindVariableFeatures(MBhealthy, selection.method = "vst", nfeatures = 2000)
MBhealthy <- ScaleData(MBhealthy)

MKinjured <- NormalizeData(MKinjured, normalization.method = "LogNormalize", scale.factor = 10000)
MKinjured <- FindVariableFeatures(MKinjured, selection.method = "vst", nfeatures = 2000)
MKinjured <- ScaleData(MKinjured)

MBinjured <- NormalizeData(MBinjured, normalization.method = "LogNormalize", scale.factor = 10000)
MBinjured <- FindVariableFeatures(MBinjured, selection.method = "vst", nfeatures = 2000)
MBinjured <- ScaleData(MBinjured)

#PCA
MKhealthy_PCA <- RunPCA(MKhealthy, features = VariableFeatures(object = MKhealthy))
MBhealthy_PCA <- RunPCA(MBhealthy, features = VariableFeatures(object = MBhealthy))

MKinjured_PCA <- RunPCA(MKinjured, features = VariableFeatures(object = MKinjured))
MBinjured_PCA <- RunPCA(MBinjured, features = VariableFeatures(object = MBinjured))


#UMAP
MKhealthy_PCA <- FindNeighbors(MKhealthy_PCA, dims = 1:43)
MKhealthy_PCA <- FindClusters(MKhealthy_PCA, resolution = 1.0)
MKhealthy_UMAP <- RunUMAP(MKhealthy_PCA, dims = 1:43)
saveRDS(MKhealthy_UMAP, file = "./030720/MKhealthy_UMAP.RDS")

MKhealthy_nExp_poi <- round(0.05*length(MKhealthy_UMAP$orig.ident)) 
print("MKhealthy")
print(MKhealthy_nExp_poi)

MBhealthy_PCA <- FindNeighbors(MBhealthy_PCA, dims = 1:43)
MBhealthy_PCA <- FindClusters(MBhealthy_PCA, resolution = 1.0)
MBhealthy_UMAP <- RunUMAP(MBhealthy_PCA, dims = 1:43)
saveRDS(MBhealthy_UMAP, file = "./030720/MBhealthy_UMAP.RDS")

MBhealthy_nExp_poi <- round(0.05*length(MBhealthy_UMAP$orig.ident)) 
print("MBhealthy")
print(MBhealthy_nExp_poi)

MKinjured_PCA <- FindNeighbors(MKinjured_PCA, dims = 1:43)
MKinjured_PCA <- FindClusters(MKinjured_PCA, resolution = 1.0)
MKinjured_UMAP <- RunUMAP(MKinjured_PCA, dims = 1:43)
saveRDS(MKinjured_UMAP, file = "./030720/MKinjured_UMAP.RDS")

MKinjured_nExp_poi <- round(0.05*length(MKinjured_UMAP$orig.ident)) 
print("MKinjured")
print(MKinjured_nExp_poi)

MBinjured_PCA <- FindNeighbors(MBinjured_PCA, dims = 1:43)
MBinjured_PCA <- FindClusters(MBinjured_PCA, resolution = 1.0)
MBinjured_UMAP <- RunUMAP(MBinjured_PCA, dims = 1:43)
saveRDS(MBinjured_UMAP, file = "./030720/MBinjured_UMAP.RDS")

MBinjured_nExp_poi <- round(0.05*length(MBinjured_UMAP$orig.ident)) 
print("MBinjured")
print(MBinjured_nExp_poi)


#pK Identification (no ground-truth)
#select pc according to elbow plot
sweep.res.list_MKhealthy <- paramSweep_v3(MKhealthy_UMAP, PCs = 1:43, sct = FALSE)
sweep.stats_MKhealthy <- summarizeSweep(sweep.res.list_MKhealthy, GT = FALSE)
bcmvn_MKhealthy <- find.pK(sweep.stats_MKhealthy)

sweep.res.list_MBhealthy <- paramSweep_v3(MBhealthy_UMAP, PCs = 1:43, sct = FALSE)
sweep.stats_MBhealthy <- summarizeSweep(sweep.res.list_MBhealthy, GT = FALSE)
bcmvn_MBhealthy <- find.pK(sweep.stats_MBhealthy)

sweep.res.list_MKinjured <- paramSweep_v3(MKinjured_UMAP, PCs = 1:43, sct = FALSE)
sweep.stats_MKinjured <- summarizeSweep(sweep.res.list_MKinjured, GT = FALSE)
bcmvn_MKinjured <- find.pK(sweep.stats_MKinjured)

sweep.res.list_MBinjured <- paramSweep_v3(MBinjured_UMAP, PCs = 1:43, sct = FALSE)
sweep.stats_MBinjured <- summarizeSweep(sweep.res.list_MBinjured, GT = FALSE)
bcmvn_MBinjured <- find.pK(sweep.stats_MBinjured)

# set pk according to the pk plot
MKhealthy_nExp_poi <- round(0.05*length(MKhealthy_UMAP$orig.ident))   ## Assuming 5% doublet formation rate(set your rate according to 10X manual)
MBhealthy_nExp_poi <- round(0.05*length(MBhealthy_UMAP$orig.ident))
MKinjured_nExp_poi <- round(0.05*length(MKinjured_UMAP$orig.ident))
MBinjured_nExp_poi <- round(0.05*length(MBinjured_UMAP$orig.ident))

MKhealthy_rmD <- doubletFinder_v3(MKhealthy_UMAP, PCs = 1:43, pN = 0.25, pK = 0.06, nExp = MKhealthy_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
MBhealthy_rmD <- doubletFinder_v3(MBhealthy_UMAP, PCs = 1:43, pN = 0.25, pK = 0.06, nExp = MBhealthy_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
MKinjured_rmD <- doubletFinder_v3(MKinjured_UMAP, PCs = 1:43, pN = 0.25, pK = 0.06, nExp = MKinjured_nExp_poi, reuse.pANN = FALSE, sct = FALSE)
MBinjured_rmD <- doubletFinder_v3(MBinjured_UMAP, PCs = 1:43, pN = 0.25, pK = 0.06, nExp = MBinjured_nExp_poi, reuse.pANN = FALSE, sct = FALSE)

plot_MKhealthy_rmD <- DimPlot(MKhealthy_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_307")
ggsave(filename = "./030720/MKhealthy_DF_UMAP.jpeg", plot = plot_MKhealthy_rmD, width = 25, units = "cm")

plot_MBhealthy_rmD <- DimPlot(MBhealthy_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_160")
ggsave(filename = "./030720/MBhealthy_DF_UMAP.jpeg", plot = plot_MBhealthy_rmD, width = 25, units = "cm")

plot_MKinjured_rmD <- DimPlot(MKinjured_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_242")
ggsave(filename = "./030720/MKinjured_DF_UMAP.jpeg", plot = plot_MKinjured_rmD, width = 25, units = "cm")

plot_MBinjured_rmD <- DimPlot(MBinjured_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_287")
ggsave(filename = "./030720/MBinjured_DF_UMAP.jpeg", plot = plot_MBinjured_rmD, width = 25, units = "cm")

# filter doublets
MKhealthy_rmD <- subset(MKhealthy_rmD, subset = DF.classifications_0.25_0.06_307 =="Singlet")
MBhealthy_rmD <- subset(MBhealthy_rmD, subset = DF.classifications_0.25_0.06_160 =="Singlet")
MKinjured_rmD <- subset(MKinjured_rmD, subset = DF.classifications_0.25_0.06_242 =="Singlet")
MBinjured_rmD <- subset(MBinjured_rmD, subset = DF.classifications_0.25_0.06_287 =="Singlet")

plot_MKhealthy_rmD_UMAP <- DimPlot(MKhealthy_rmD, reduction = "umap")
ggsave(filename = "./030720/MKhealthy_rmD_UMAP.jpeg", plot = plot_MKhealthy_rmD_UMAP, width = 25, units = "cm")

plot_MBhealthy_rmD_UMAP <- DimPlot(MBhealthy_rmD, reduction = "umap")
ggsave(filename = "./030720/MBhealthy_rmD_UMAP.jpeg", plot = plot_MBhealthy_rmD_UMAP, width = 25, units = "cm")

plot_MKinjured_rmD_UMAP <- DimPlot(MKinjured_rmD, reduction = "umap")
ggsave(filename = "./030720/MKinjured_rmD_UMAP.jpeg", plot = plot_MKinjured_rmD_UMAP, width = 25, units = "cm")

plot_MBinjured_rmD_UMAP <- DimPlot(MBinjured_rmD, reduction = "umap")
ggsave(filename = "./030720/MBinjured_rmD_UMAP.jpeg", plot = plot_MBinjured_rmD_UMAP, width = 25, units = "cm")
