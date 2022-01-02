library(DoubletFinder)
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Load indie datasets
MLhealthy_rmD <- readRDS(file = "./110720/Bleo1_MLhealthy_rmD.RDS")
MBhealthy_rmD <- readRDS(file = "./110720/Bleo1_MBhealthy_rmD.RDS")

MLbleo_rmD <- readRDS(file = "./110720/Bleo1_MLbleo_rmD.RDS")
MBbleo_rmD <- readRDS(file = "./110720/Bleo1_MBbleo_rmD.RDS")

# Can't read from connection MBbleo_rmD therefore
#MBbleo_UMAP <- readRDS(file = "./110720/MBbleo_UMAP.RDS")
#MBbleo_nExp_poi <- round(0.05*length(MBbleo_UMAP$orig.ident))
#MBbleo_rmD <- doubletFinder_v3(MBbleo_UMAP, PCs = 1:45, pN = 0.25, pK = 0.06, nExp = MBbleo_nExp_poi, reuse.pANN = FALSE, sct = FALSE)


# DF Plots
plot_MLhealthy_DF <- DimPlot(MLhealthy_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_182")
ggsave(filename = "./110720/Bleo1_MLhealthy_DF_UMAP.jpeg", plot = plot_MLhealthy_DF, width = 30, units = "cm")

plot_MBhealthy_DF <- DimPlot(MBhealthy_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_159")
ggsave(filename = "./110720/Bleo1_MBhealthy_DF_UMAP.jpeg", plot = plot_MBhealthy_DF, width = 30, units = "cm")

plot_MLbleo_DF <- DimPlot(MLbleo_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_142")
ggsave(filename = "./110720/Bleo1_MLbleo_DF_UMAP.jpeg", plot = plot_MLbleo_DF, width = 30, units = "cm")

plot_MBbleo_DF <- DimPlot(MBbleo_rmD, reduction = "umap",group.by = "DF.classifications_0.25_0.06_178")
ggsave(filename = "./110720/Bleo1_MBbleo_DF_UMAP.jpeg", plot = plot_MBbleo_DF, width = 30, units = "cm")

# filter doublets
MLhealthy_rmD <- subset(MLhealthy_rmD, subset = DF.classifications_0.25_0.06_182 =="Singlet")
MBhealthy_rmD <- subset(MBhealthy_rmD, subset = DF.classifications_0.25_0.06_159 =="Singlet")
MLbleo_rmD <- subset(MLbleo_rmD, subset = DF.classifications_0.25_0.06_142 =="Singlet")
MBbleo_rmD <- subset(MBbleo_rmD, subset = DF.classifications_0.25_0.06_178 =="Singlet")

plot_MLhealthy_rmD_UMAP <- DimPlot(MLhealthy_rmD, reduction = "umap")
ggsave(filename = "./110720/Bleo1_MLhealthy_rmD_UMAP.jpeg", plot = plot_MLhealthy_rmD_UMAP, width = 30, units = "cm")

plot_MBhealthy_rmD_UMAP <- DimPlot(MBhealthy_rmD, reduction = "umap")
ggsave(filename = "./110720/Bleo1_MBhealthy_rmD_UMAP.jpeg", plot = plot_MBhealthy_rmD_UMAP, width = 30, units = "cm")

plot_MLbleo_rmD_UMAP <- DimPlot(MLbleo_rmD, reduction = "umap")
ggsave(filename = "./110720/Bleo1_MLbleo_rmD_UMAP.jpeg", plot = plot_MLbleo_rmD_UMAP, width = 30, units = "cm")

plot_MBbleo_rmD_UMAP <- DimPlot(MBbleo_rmD, reduction = "umap")
ggsave(filename = "./110720/Bleo1_MBbleo_rmD_UMAP.jpeg", plot = plot_MBbleo_rmD_UMAP, width = 30, units = "cm")

#Over-write 110720_SSH_Bleo_DF files with subsetted datasets
saveRDS(MLhealthy_rmD, file = "./110720/Bleo1_MLhealthy_rmD.RDS")
saveRDS(MBhealthy_rmD, file = "./110720/Bleo1_MBhealthy_rmD.RDS")
saveRDS(MLbleo_rmD, file = "./110720/Bleo1_MLbleo_rmD.RDS")
saveRDS(MBbleo_rmD, file = "./110720/Bleo1_MBbleo_rmD.RDS")
