library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

path <-"[datasetDF_path]"

dataset_DF_n <- readRDS(path, "[RDS Object name.RDS]")


# Merge datasets into experiment
Experiment <- merge(dataset_DF_1, y = c(dataset_DF_2, dataset_DF_3, dataset_DF_4), add.cell.ids = c("dataset_DF_1", "dataset_DF_2", "dataset_DF_3", "dataset_DF_4"), project = "Experiment") 
unique(sapply(X = strsplit(colnames(Experiment), split = "_"), FUN = "[", 1))

# set %mt
Experiment[["percent.mt"]] <- PercentageFeatureSet(Experiment, pattern = "^MT-")

# QC
plot_QC_Vlnplot <- VlnPlot(Experiment, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave(filename = "[desired path and name].png", width = 25, units = "cm")

plot_QC_Scatter1 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot_QC_Scatter2 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2
ggsave(filename = "[desired path]QCscatter.png", width = 25, units = "cm")

#Filter cells 
Experiment <- subset(Experiment, subset = nFeature_RNA > [desired cut-off] & percent.mt < [desired cut-off] & nCount_RNA < [desired cut-off])

#Normalise
Experiment <-  NormalizeData(Experiment, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
Experiment <- FindVariableFeatures(Experiment, selection.method = "vst", nfeatures = 2000)
Experiment_top15 <- head(VariableFeatures(Experiment), 15)
Plot_FS_15 <- LabelPoints(plot = VariableFeaturePlot(Experiment), points = Experiment_top15, repel = TRUE)
ggsave(filename = "[desired path and name].png", width = 30, units = "cm")

#Scale 
Experiment <- ScaleData(Experiment)

```

```{r}

# Plot post-QC graphs

plot_postQC_Vlnplot <- VlnPlot(Experiment, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave("[desired path and name].png", width = 24, units ="cm")

plot_postQC_Scatter1 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot_postQC_Scatter2 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot_postQC_Scatter <- plot_postQC_Scatter1 + plot_postQC_Scatter2
ggsave("[desired path and name].png",width = 24, units ="cm")

#PCA
Experiment <- RunPCA(Experiment, features = VariableFeatures(object = Experiment))
capture.output(print(Experiment[["pca"]], dims = 1:10, nfeatures = 5), file = "[desired path and name].txt")

saveRDS(Experiment, file = "[desired path and name].RDS")

#n(PC) Plots
png(file="[desired path and name].jpeg", width = 30, height = 60, units = "cm", res = 800)
DimHeatmap(Experiment, dims = [lower_dim:higher_dim], cells = 500, balanced = TRUE, ncol = 3)
dev.off()

plot_Elbow <- ElbowPlot(Experiment, ndims = 50)
ggsave(filename = "[desired path and name].png", plot = plot_Elbow, width = 30, units = "cm")

# UMAP
Experiment <- FindNeighbors(Experiment, dims = 1:[selected_dim])
Experiment <- FindClusters(Experiment, resolution = 1.00)
Experiment <- RunUMAP(Experiment, dims = 1:[selected_dim])

saveRDS(Experiment, file = "[desired path and name].RDS")

plot_UMAP <- DimPlot(Experiment, reduction = "umap", group.by = "orig.ident")
ggsave("[desired path and name].png")

plot_UMAP_clusters <- DimPlot(Experiment, reduction = "umap", label = T)
ggsave("[desired path and name].png")


# Marker Genes

Experiment.markers <- FindAllMarkers(Experiment, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Experiment.markersOUTPUT <- Experiment.markers %>% group_by(cluster) %>% arrange(pct.2, .by_group = T)
write.csv(Experiment.markersOUTPUT, file = "[desired path and name].csv")


Plot_FP_QCmetrics <- FeaturePlot(Experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("[desired path and name].png", width = 35, height = 15, units ="cm")


# Plot big markers

plot_celltype <- FeaturePlot(Experiment, features = c("gene1", "gene2", "gene3"), ncol = 3, pt.size = 0.4)
ggsave(filename = "[desired path and name].png", plot = plot_MonoMacs, width = 30, height = 22, units = "cm")
