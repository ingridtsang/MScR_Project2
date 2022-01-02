# Load required Libraries

library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)

#Load directory
mA5t3 <- readRDS("./230720_mA5Test2_UMAP.RDS")
### cell type annotation ###
cluster.annotation <- c("Blood B cell 1", "Kidney macrophage 1", "Kidney macrophage 2", "UUO blood Ly6Chi monocyte", 
                        "Blood Ly6Clo monocyte", "CD8+ T cell 1", "Kidney macrophage 3", "Kidney macrophage 4", 
                        "Cd4+ T cell 1", "NK cell 1", "Kidney monocyte", "Healthy blood Ly6Chi monoyte", "Kidney cDC2 1", 
                        "Cd4+ T cell 2", "Kidney cDC2 2", "Kidney B cell", "Kidney cDC1", 
                        "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", 
                        "Kidney IFN-primed macrophage", "Kidney endothelial cells", "Blood neutrophils", "NK cell 2", 
                        "pDC","Neutrophil/B cell doublet","Basophil","Proliferating","Neutrophil/T cell doublet", "ILC2", 
                        "CCR7+ cDC", "Cd8+ T cell 2", "Monocyte/T cell doublet", "Platelet")

names(cluster.annotation) <- levels(mA5t3)
mA5t3 <- RenameIdents(mA5t3, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(mA5t3),
                        celltype = Idents(mA5t3))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
mA5t3 <- AddMetaData(mA5t3, cell.data, col.name = "celltype")

# save annotated UMAP
annotated.umap.plot <- DimPlot(mA5t3, reduction = "umap", split.by = "orig.ident", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
ggsave(filename = "../mA5t3ANNOTATIONTRIAL_UMAP_OGID.png", plot = annotated.umap.plot, width = 30, units = "cm")

### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(mA5t3),
                        celltype = Idents(mA5t3))

lineage.annotation <- c("B cell",	"MP",	"MP",	"MP",	"MP",	"T cell/ILC",	"MP",	"MP",	"T cell/ILC",	"T cell/ILC",	
                          "MP",	"MP",	"MP",	"T cell/ILC",	"MP",	"B cell",	"MP",	"MP",	"MP",	"MP",	"Endothelia",	
                          "Neutrophil",	"T cell/ILC",	"MP",	"Doublet",	"Basophil",	"Proliferating",	"Doublet",	
                          "T cell/ILC",	"MP",	"T cell/ILC",	"Doublet",	"Platelet")


lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
mA5t3 <- AddMetaData(mA5t3, meta.data, col.name = "lineage")
# save annotated UMAP
lineage.umap.plot <- DimPlot(mA5t3, reduction = "umap", group.by = "lineage", split.by = "orig.ident", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = "../mA5t3ANNOTATIONTRIAL_UMAPlineageOGID.png", plot = lineage.umap.plot, width = 30, units = "cm")
