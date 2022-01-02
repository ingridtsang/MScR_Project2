library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

# Load Data @PCA
mA5 <- readRDS("./Kidney_mA5-GH/Test1_Outputs/170720_mA5_PCA.RDS")

# UMAP
mA5t2 <- FindNeighbors(mA5, dims = 1:45)
mA5t2 <- FindClusters(mA5t2, resolution = 1.00)
mA5t2 <- RunUMAP(mA5t2, dims = 1:45)

saveRDS(mA5t2, file = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5Test2_UMAP.RDS")

plot_UMAP <- DimPlot(mA5t2, reduction = "umap", group.by = "orig.ident")
ggsave("./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_UMAPogid.png", plot = plot_UMAP)

plot_UMAP_clusters <- DimPlot(mA5t2, reduction = "umap", label = T)
ggsave("./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_UMAPclusters.png", plot = plot_UMAP_clusters)

# Plot post-QC graphs

plot_QC_Vlnplot <- VlnPlot(mA5t2, features= c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
ggsave("./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_QCVln.png", plot = plot_QC_Vlnplot, width = 24, units ="cm")

plot_QC_Scatter1 <- FeatureScatter(mA5t2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot_QC_Scatter2 <- FeatureScatter(mA5t2, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot_QC_Scatter <- plot_QC_Scatter1 + plot_QC_Scatter2
ggsave("./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_QCScatter.png", plot = plot_QC_Scatter, width = 24, units ="cm")


Plot_FP_QCmetrics <- FeaturePlot(mA5t2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggsave("./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_QCFP.png", plot = Plot_FP_QCmetrics, width = 35, height = 15, units ="cm")


# Marker Genes

mA5t2.markers <- FindAllMarkers(mA5t2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mA5t2.markersOUTPUT <- mA5t2.markers %>% group_by(cluster) %>% arrange(pct.2, .by_group = T)
write.csv(mA5t2.markersOUTPUT, file = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5Test2_biomarkers.csv")


# Plot big markers

plot_MonoMacs <- FeaturePlot(mA5t2, features = c("Csf1r", "Lyz2", "Itgam", "Cd68", "S100a8", "S100a4"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_MonoMacs.png", plot = plot_MonoMacs, width = 30, height = 22, units = "cm")

plot_MonoLy6c <- FeaturePlot(mA5t2, features = c("Ly6c2", "Ccr2", "Vcan", "Ace", "Treml4"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_MonoLy6c.png", plot = plot_MonoLy6c, width = 30, height = 22, units = "cm")

plot_NeuTNKB <- FeaturePlot(mA5t2, features = c("Ly6g", "Mpo", "Cd3d", "Cd3e", "Cd8a", "Cd4", "Klra9", "Klrb1b", "Gzma", "Cd79a", "Cd79b"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_NeuTNKB.png", plot = plot_NeuTNKB, width = 30, height = 44, units = "cm")

plot_bloods_cDC <- FeaturePlot(mA5t2, features = c("Fcer1a", "Il4", "Cpa3", "Kdr", "Pecam1", "Icam2", "Cdh5", "Hba-a1", "Hbb-a2", "Xcr1", "Cd209a", "Siglech"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_bloodsCDC.png", plot = plot_bloods_cDC, width = 30, height = 44, units = "cm")

plot_MonoLy6cHi <- FeaturePlot(mA5t2, features = c("Ly6c2", "Ccr2", "Vcan"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_MonoLy6cHIGH.png", plot = plot_MonoLy6cHi, width = 30, height = 11, units = "cm")

plot_MonoLy6cLo <- FeaturePlot(mA5t2, features = c("Ace", "Treml4"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_MonoLy6cLOW.png", plot = plot_MonoLy6cLo, width = 20, height = 11, units = "cm")

plot_Neu <- FeaturePlot(mA5t2, features = c("Ly6g", "Mpo"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_Neu.png", plot = plot_Neu, width = 20, height = 11, units = "cm")

plot_T <- FeaturePlot(mA5t2, features = c("Cd3d", "Cd3e", "Cd8a", "Cd4"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_T.png", plot = plot_T, width = 20, height = 22, units = "cm")

plot_NK <- FeaturePlot(mA5t2, features = c( "Klra9", "Klrb1b", "Gzma"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_NK.png", plot = plot_NK, width = 30, height = 11, units = "cm")

plot_B <- FeaturePlot(mA5t2, features = c("Cd79a", "Cd79b"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_B.png", plot = plot_B, width = 30, height = 11, units = "cm")

plot_DC <- FeaturePlot(mA5t2, features = c("Xcr1", "Cd209a", "Siglech"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_DC.png", plot = plot_DC, width = 30, height = 11, units = "cm")

plot_baso <- FeaturePlot(mA5t2, features = c("Fcer1a", "Il4", "Cpa3"), ncol = 3, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_baso.png", plot = plot_baso, width = 30, height = 11, units = "cm")

plot_EC <- FeaturePlot(mA5t2, features = c("Kdr", "Pecam1", "Icam2", "Cdh5"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_EC.png", plot = plot_EC, width = 20, height = 22, units = "cm")

plot_HBC <- FeaturePlot(mA5t2, features = c("Hba-a1", "Hbb-a2"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_HBC.png", plot = plot_HBC, width = 20, height = 11, units = "cm")

plot_Macro <- FeaturePlot(mA5t2, features = c("Mmp12", "Ccl12" , "Mgl2", "Pf4"), ncol = 2, pt.size = 0.4)
ggsave(filename = "./Kidney_mA5-GH/Test2_Outputs/230720_mA5t2_FP_Macro.png", plot = plot_Macro, width = 20, height = 22, units = "cm")


