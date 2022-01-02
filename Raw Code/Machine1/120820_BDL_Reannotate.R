```{r}
# Load required Libraries

library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)
library(data.table)
library(plyr)

#Set paths
BDL2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
BDL2Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/Outputs/"

#Set date

Today.date <- "070820" 

# Set colour palette

BDL2.colpal <- c("#8EB26E",	"#88BCB2",	"#33A02C",	"#00CDB3",	"#C8828D",	"#F4A1AF",	"#B2DF8A",	"#AFF1E1")

# Load Data
BDL2 <- readRDS(paste0(BDL2.path, Today.date, "_BDL2_UMAP.RDS"))

BDL2.markers <- read.csv(paste0(BDL2Output.path, Today.date, "_BDL2_ALLmarkers.csv"))
```


```{r}
### cell type annotation ###
BDL2_cluster.annotation <- c("Liver Kupffer Cell 1",	"Blood Ly6Clo monocyte 1",	"BDL Blood Ly6Chi monocyte 1",	"Healthy Blood Ly6Chi monocyte 1",	"Liver macrophage 2",	"B cell 1",	"Liver macrophage 4",	"Liver macrophage 1",	"B cell 2",	"NK cell 1",	"BDL Blood Ly6Clo monocyte 1",	"CD8+ T cell 1",	"T cell 2",	"Cd4+ T cell 1",	"Proliferating",	"cDC2",	"NK cell 2",	"Neutrophil 1",	"Healthy Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Chi monocyte 2",	"Liver Endothelial Cells",	"Liver macrophage 3",	"Platelet 1",	"pDC",	"Blood Ly6Clo monocyte 2",	"B cell 3",	"Blood IFN-primed Ly6Chi monocyte",	"Liver Kupffer Cell 2",	"cDC1",	"Neutrophil 2",	"Liver macrophage 5",	"Monocyte/B cell doublet",	"Basophil",	"Monocyte/B cell doublet 1",	"Monocyte/T cell doublet 2",	"Platelet 2")

names(BDL2_cluster.annotation) <- levels(BDL2)
BDL2 <- RenameIdents(BDL2, BDL2_cluster.annotation)

BDL2_cell.data <- data.table(barcode = colnames(BDL2),
                        Cell.type = Idents(BDL2))

BDL2_cell.data <- data.frame(BDL2_cell.data, row.names = BDL2_cell.data$barcode)
BDL2_cell.data$barcode <- NULL
BDL2 <- AddMetaData(BDL2, BDL2_cell.data, col.name = "Cell.type")

### save annotated UMAP ###
annotated.umap.plot.split <- DimPlot(BDL2, reduction = "umap", split.by = "Names", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_cluster_OGIDsplit.png"), plot = annotated.umap.plot.split, width = 45, height = 10, units = "cm")

annotated.umap.plot.merge <- DimPlot(BDL2, reduction = "umap", label = T, label.size = 2, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_cluster.png"), plot = annotated.umap.plot.merge, width = 30, units = "cm")


### cell lineage annotation ###
BDL2_cell.data <- data.table(barcode = colnames(BDL2),
                        Cell.type = Idents(BDL2))

BDL2_lineage.annotation <- c("MP",	"MP",	"MP",	"MP",	"MP",	"B cell",	"MP",	"MP",	"B cell",	"T cell/ILC",	"MP",	"T cell/ILC",	"T cell/ILC",	"T cell/ILC",	"Prolferating",	"MP",	"MP",	"MP",	"MP",	"MP",	"Endothelia",	"MP",	"Platelet",	"MP",	"MP",	"B cell",	"MP",	"MP",	"MP",	"Neutrophil",	"MP",	"Doublet",	"Basophil",	"Doublet",	"Doublet",	"Platelet")


BDL2_lineage.data <- data.table(Cell.type = BDL2_cluster.annotation, lineage = BDL2_lineage.annotation)

meta.data <- merge(BDL2_cell.data, BDL2_lineage.data, by = "Cell.type")

meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$Cell.type <- NULL

BDL2 <- AddMetaData(BDL2, meta.data, col.name = "lineage")

# save annotated UMAP
lineage.umap.plot <- DimPlot(BDL2, reduction = "umap", group.by = "lineage", split.by = "Names", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_lineage_OGIDsplit.png"), plot = lineage.umap.plot, width = 45, height = 10, units = "cm")

lineage.umap.plot <- DimPlot(BDL2, reduction = "umap", group.by = "lineage", label = T, label.size = 3, pt.size = 0.01) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, Today.date, "_BDL2_aUMAP_lineage.png"), plot = lineage.umap.plot, width = 30, units = "cm")

saveRDS(BDL2, file = paste0(BDL2.path, Today.date, "_BDL2_annotated.RDS"))


```


```{r}
#annotate

annotate <- read.table(file = paste0(BDL2.path,"BDL2_Annotations",".csv"), header = T, sep = ",", row.names = 1)

map_Celltype <- function(row){
  cluster <- row[7]
  cluster <- as.numeric(cluster) + 1
  print(cluster)
  Cell.type <- annotate$Cell.type[cluster]
  print(Cell.type)
  return(Cell.type)
}

BDL2.markers <- data.frame(BDL2.markers)

BDL2.markers$Cell.type <- apply(X = BDL2.markers,1, FUN = map_Celltype)


map_Lineage <- function(row){
  cluster <- row[7]
  cluster <- as.numeric(cluster) + 1
  lineage <- annotate$lineage[cluster]
  return(lineage)
}

BDL2.markers <- data.frame(BDL2.markers)

BDL2.markers$Lineage <- apply(X = BDL2.markers,1, FUN = map_Lineage)


```


