### Before starting, Command-F "dataset" and "Experiment" and replace with actual name of your dataset/experiments. Command-F "[[desired path and name]"  and change filename of outputs manually. Do the same with anything in square brackets [] under "Set variables". ### 

#Load packages
library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

#Set Variables
path <-"[datasetDF_path]" #Set path of datasets which have had doublets removed in DoubletFinder (DF)
dataset_DF_n <- readRDS(path, "[desired path and name.RDS]") # read in the RDS files of each DF-cleaned dataset 


# Merge datasets into one experiment. Four datasets have been specified as an example but add/delete as necessary in "y = c(...)" and "add.cekk.ids = c(...)"
Experiment <- merge(dataset_DF_1, y = c(dataset_DF_2, dataset_DF_3, dataset_DF_4), add.cell.ids = c("dataset_DF_1", "dataset_DF_2", "dataset_DF_3", "dataset_DF_4"), project = "Experiment") 
unique(sapply(X = strsplit(colnames(Experiment), split = "_"), FUN = "[", 1))

# set %mitochondrial genes. "^MT-" for animal datasets and "^mt-" for human (I think) 
Experiment[["percent.mt"]] <- PercentageFeatureSet(Experiment, pattern = "^MT-")

# Quality Control Plots
plot_QC_Vlnplot <- VlnPlot(Experiment, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident") #Create Violin Plots of 3 QC metrics in one image
ggsave(filename = "[desired path and name].png", width = 25, units = "cm") #Save as png

plot_QC_Scatter1 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") #Create scatterplot of nCountRNA vs nFeature_RNA
plot_QC_Scatter2 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident") #Create scatterplot of nCountRNA vs percent.mt (i.e. percentage of mitochondrial genes)
plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2 #Combine the 2 scatter plots into one object/image
ggsave(filename = "[desired path and name].png", width = 25, units = "cm") #Save as png
```
# Using the QC plots, set desired min and max cut-offs for each QC metric and input into next line of code. 

```{r}
#Filter cells 
Experiment <- subset(Experiment, subset = nFeature_RNA > [desired cut-off] & percent.mt < [desired cut-off] & nCount_RNA < [desired cut-off]) #input "Desired cut-off"

#Normalise
Experiment <-  NormalizeData(Experiment, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
Experiment <- FindVariableFeatures(Experiment, selection.method = "vst", nfeatures = 2000) 
Experiment_top15 <- head(VariableFeatures(Experiment), 15) #Set the top15 most variable genes
Plot_FS_15 <- LabelPoints(plot = VariableFeaturePlot(Experiment), points = Experiment_top15, repel = TRUE) #Plot a scatter of variability in gene expression, with the top15 most variable genes labelled
ggsave(filename = "[desired path and name].png", width = 30, units = "cm")

#Scale 
Experiment <- ScaleData(Experiment)

```

```{r}

# Plot graphs to check how the dataset has been affected by QC

plot_postQC_Vlnplot <- VlnPlot(Experiment, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave("[desired path and name].png", width = 24, units ="cm")

plot_postQC_Scatter1 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot_postQC_Scatter2 <- FeatureScatter(Experiment, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot_postQC_Scatter <- plot_postQC_Scatter1 + plot_postQC_Scatter2
ggsave("[desired path and name].png",width = 24, units ="cm")

#Run PCA analysis and save the experiment as an RDS
Experiment <- RunPCA(Experiment, features = VariableFeatures(object = Experiment))
capture.output(print(Experiment[["pca"]], dims = 1:10, nfeatures = 5), file = "[desired path and name].txt")

saveRDS(Experiment, file = "[desired path and name].RDS")

#n(PC) Plots to figure out which dimention should be included in future analysis. 
png(file="[desired path and name].jpeg", width = 30, height = 60, units = "cm", res = 800)
DimHeatmap(Experiment, dims = [lower_dim:higher_dim], cells = 500, balanced = TRUE, ncol = 3) #Set a range between "lower_dim" and "higher_dim" (inclusive) to check which dimension(s) are(n't) affecting your data. Range should ideally be no bigger than 20. Repeat up to dim 50 max
dev.off()

plot_Elbow <- ElbowPlot(Experiment, ndims = 50) #Elbow plot is a quicker and less computationally-intensive way of figuring out which dimension(s) to include in UMAP analysis
ggsave(filename = "[desired path and name].png", plot = plot_Elbow, width = 30, units = "cm")
```

```{r}

# Make UMAP plots, repeat as necessary with different "[selected-dim]" to find the right clustering
Experiment <- FindNeighbors(Experiment, dims = 1:[selected_dim])
Experiment <- FindClusters(Experiment, resolution = 1.00)
Experiment <- RunUMAP(Experiment, dims = 1:[selected_dim])

saveRDS(Experiment, file = "[desired path and name].RDS")

plot_UMAP <- DimPlot(Experiment, reduction = "umap", group.by = "orig.ident") #UMAP grouped by sample/datasets
ggsave("[desired path and name].png")

plot_UMAP_clusters <- DimPlot(Experiment, reduction = "umap", label = T) #UMAP grouped by cluster
ggsave("[desired path and name].png")


# Find Marker Genes
Experiment.markers <- FindAllMarkers(Experiment, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Experiment.markersOUTPUT <- Experiment.markers %>% group_by(cluster) %>% arrange(pct.2, .by_group = T)
write.csv(Experiment.markersOUTPUT, file = "[desired path and name].csv")

#Plot QC metrics onto clusters
Plot_FP_QCmetrics <- FeaturePlot(Experiment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("[desired path and name].png", width = 35, height = 15, units ="cm")


# Plot big markers
plot_celltype <- FeaturePlot(Experiment, features = c("gene1", "gene2", "gene3"), ncol = 3, pt.size = 0.4)
ggsave(filename = "[desired path and name].png", plot = plot_MonoMacs, width = 30, height = 22, units = "cm")
