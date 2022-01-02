library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)



path = "[dataset_path]"

dataset.data <- Read10X(data.dir = path)
dataset <- CreateSeuratObject(counts = dataset.data, project = "[dataset_name]")


# Set %mt

dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^mt-")

# LogNorm, FS, Scale
dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
dataset <- ScaleData(dataset)

# PCA
dataset_PCA <- RunPCA(dataset, features = VariableFeatures(object = dataset))

# UMAP -> Use dims = 50 for now
dataset_PCA <- FindNeighbors(dataset_PCA, dims = 1:50)
dataset_PCA <- FindClusters(dataset_PCA, resolution = 1.0)
dataset_UMAP <- RunUMAP(dataset_PCA, dims = 1:50)

# set pk according to the pk plot -> Need to figure out nExp_poi so you can 
dataset_nExp_poi <- round(0.05*length(dataset_UMAP$orig.ident))   ## Assuming 5% doublet formation rate(set your rate according to 10X manual)

# pK Identification (no ground-truth)
#select pc according to elbow plot
sweep.res.list_dataset_UMAP <- paramSweep_v3(dataset_UMAP, PCs = 1:50, sct = FALSE)
sweep.stats_dataset_UMAP <- summarizeSweep(sweep.res.list_dataset_UMAP, GT = FALSE)
bcmvn_dataset_UMAP <- find.pK(sweep.stats_dataset_UMAP)

#-----------save dataset_UMAP--------------------

# Doublet Finder
dataset_rmD <- doubletFinder_v3(dataset_UMAP, PCs = 1:50, pN = 0.25, pK = 0.06, nExp = dataset_nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Make graphs - DoubletFinder
plot_dataset_DF <- DimPlot(dataset_rmD, reduction = "umap", group.by = paste0("DF.classifications_0.25_0.06_", dataset_nExp_poi))
ggsave(filename = "[filename].png", width = 30, units = "cm")

# Remove doublets input dataset_nExp_poi by hand
code = paste0("DF.classifications_0.25_0.06_", dataset_nExp_poi) #test to see if this works
dataset_rmD <- subset(dataset_rmD, subset = code =="Singlet")
saveRDS(dataset_rmD, file = "[filename].RDS")

# Make graphs - OG
plot_dataset_OG <- DimPlot(dataset_UMAP, reduction = "umap")
ggsave(filename = "[filename].jpeg", plot = plot_dataset_OG, width = 30, units = "cm")

# Make graphs - rmD
plot_dataset_rmD_UMAP <- DimPlot(dataset_rmD, reduction = "umap")
ggsave(filename = "[filename].jpeg", plot = plot_dataset_rmD_UMAP, width = 30, units = "cm")

#-----------save dataset_rmD--------------------