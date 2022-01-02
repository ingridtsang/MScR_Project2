# Load library
library(uwot)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(cowplot)
library(data.table)

#Set paths

mA5origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
BDL2origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
CCl4origin.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"

mA5Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Kidney_UUO/"
Bleo2Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Lung_Bleo/"
BDL2Output.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_BDL/"
CCl4Output.path  <- "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Liver_CCl4/"

#Set date

Load.date <- "070820"
Save.date <- "160820"
date.1 <- "200820"
date.2 <- "210820"

#Colour palettes
#mA5_lin.colpal <- c("#ffe5cc",	"#ffcb99",	"#ffb266",	"#ff9832",	"#ff7f00",	"#cc6500",	"#994c00",	"#663200",	"#331900")
#Bleo2_lin.colpal <-	c("#a58ac2",	"#9062c1",	"#7744ac",	"#5d3688",	"#442763",	"#3f245c")
#BDL2_lin.colpal <- c("#d6ecd4",	"#add9aa",	"#84c680",	"#5bb356",	"#33a02c",	"#2d9027",	"#23701e",	"#195016",	"#0f300d")
#CCl4_lin.colpal <- c("#d0e5f7",	"#a1cbf0",	"#73b1e9",	"#4497e2",	"#167edb",	"#1371c5",	"#0f5899",	"#0b3f6d",	"#062541")

mA5_lin.colpal <- c("#ffe5cc",	"#ffcb99",	"#ffb266",	"#ff9832",	"#ff7f00",		"#994c00",	"#663200",	"#331900")
Bleo2_lin.colpal <-	c("#a58ac2",	"#9062c1",	"#7744ac",	"#5d3688",	"#3f245c")
BDL2_lin.colpal <- c("#d6ecd4",	"#add9aa",	"#84c680",	"#5bb356",	"#33a02c",		"#23701e",	"#195016",	"#0f300d")
CCl4_lin.colpal <- c("#d0e5f7",	"#a1cbf0",	"#73b1e9",	"#167edb",	"#1371c5",	"#0f5899",	"#0b3f6d",	"#062541")
#
#
##Load annotated datasets
#amA5 <- readRDS(paste0(mA5origin.path, Load.date, "_mA5_annotated.RDS"))
#aBleo2 <- readRDS(paste0(Bleo2origin.path, Load.date, "_Bleo2_annotated.RDS"))
#aBDL2 <- readRDS(paste0(BDL2origin.path, date.2, "_BDL2_annotated.RDS"))
#aCCl4 <- readRDS(paste0(CCl4origin.path, Load.date, "_CCl4_annotated.RDS"))
#
#
##Clean datasets + save
#
#amA5_cleaned <- subset(amA5 , subset = lineage == "Doublet", invert = TRUE)
#saveRDS(amA5_cleaned, file = paste0(mA5Output.path, date.1, "_amA5_cleaned.RDS"))
#
#aBleo2_cleaned <- subset(aBleo2, subset = lineage != "Doublet")
#saveRDS(aBleo2_cleaned, file = paste0(Bleo2Output.path, date.1, "_aBleo2_cleaned.RDS"))
#
#aBDL2_cleaned <- subset(aBDL2, subset = lineage != "Doublet")
#saveRDS(aBDL2_cleaned, file = paste0(BDL2Output.path, date.1, "_aBDL2_cleaned.RDS"))
#
#aCCl4_cleaned <- subset(aCCl4, subset = lineage != "Doublet")
#saveRDS(aCCl4_cleaned, file = paste0(CCl4Output.path, date.1, "_aCCl4_cleaned.RDS"))
#
## Find and save marker gene lists
#mA5CLEANED.markers <- FindAllMarkers(amA5_cleaned, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(mA5CLEANED.markers, file = paste0(mA5Output.path, date.1, "_amA5CLEANED_ALLmarkers.csv"))
#
#Bleo2CLEANED.markers <- FindAllMarkers(aBleo2_cleaned, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(Bleo2CLEANED.markers, file = paste0(Bleo2Output.path, date.1, "_aBleo2CLEANED_ALLmarkers.csv"))
#
#BDL2CLEANED.markers <- FindAllMarkers(aBDL2_cleaned, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(BDL2CLEANED.markers, file = paste0(BDL2Output.path, date.1, "_aBDL2CLEANED_ALLmarkers.csv"))
#
#CCl4CLEANED.markers <- FindAllMarkers(aCCl4_cleaned, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(CCl4CLEANED.markers, file = paste0(CCl4Output.path, date.1, "_aCCl4CLEANED_ALLmarkers.csv"))
#
## Make function to annotate marker gene list
mapCLEANED_Lineage <- function(row){
  cluster <- as.character(row[7])
  print(cluster)
  for(i in 1:length(annotate$Cell.type)){
    print(length(annotate$Cell.type))
    annotate_cell_type <- annotate$Cell.type[i]
    print(annotate_cell_type)
    if(cluster==annotate_cell_type){
      return(as.character(annotate$lineage[i]))
      print(as.character(annotate$lineage[i]))
    }
  }
}
#
## Annotate cleaned marker gene list
#annotate <- read.csv(paste0(mA5Output.path,"amA5CLEANED_Annotations",".csv"), header = T, sep = ",", row.names = 1)
#mA5CLEANED.markers <- read.csv(paste0(mA5Output.path, date.1, "_amA5CLEANED_ALLmarkers.csv"))
#mA5CLEANED.markers$lineage <- apply(X = mA5CLEANED.markers, 1, FUN = mapCLEANED_Lineage)
#write.csv(mA5CLEANED.markers, file = paste0(mA5Output.path, date.1, "_mA5CLEANEDAnno_ALLmarkers.csv"))
#
#annotate <- read.csv(paste0(Bleo2Output.path,"aBleo2CLEANED_Annotations",".csv"), header = T, sep = ",", row.names = 1)
#Bleo2CLEANED.markers <- read.csv(paste0(Bleo2Output.path, date.1, "_aBleo2CLEANED_ALLmarkers.csv"))
#Bleo2CLEANED.markers$lineage <- apply(X = Bleo2CLEANED.markers, 1, FUN = mapCLEANED_Lineage)
#Bleo2CLEANED.markers$lineage <- sapply(Bleo2CLEANED.markers$lineage, toString)
#write.csv(Bleo2CLEANED.markers, file = paste0(Bleo2Output.path, date.1, "_aBleo2CLEANEDAnno_ALLmarkers.csv"))
#
#annotate <- read.csv(paste0(BDL2Output.path,"aBDL2CLEANED_Annotations.csv"), header = T, sep = ",", row.names = 1)
#BDL2CLEANED.markers <- read.csv(paste0(BDL2Output.path, date.1, "_aBDL2CLEANED_ALLmarkers.csv"))
#BDL2CLEANED.markers$lineage <- apply(X = BDL2CLEANED.markers, 1, FUN = mapCLEANED_Lineage)
#BDL2CLEANED.markers$lineage <- sapply(BDL2CLEANED.markers$lineage, toString)
#write.csv(BDL2CLEANED.markers, file = paste0(BDL2Output.path, date.1, "_aBDL2CLEANEDAnno_ALLmarkers.csv"))
#
#annotate <- read.csv(paste0(CCl4Output.path, "aCCl4CLEANED_Annotations.csv"), header = T, sep = ",", row.names = 1)
#CCl4CLEANED.markers <- read.csv(paste0(CCl4Output.path, date.1, "_aCCl4CLEANED_ALLmarkers.csv"))
#CCl4CLEANED.markers$lineage <- apply(X = CCl4CLEANED.markers, 1, FUN = mapCLEANED_Lineage)
#CCl4CLEANED.markers$lineage <- sapply(CCl4CLEANED.markers$lineage, toString)
#write.csv(CCl4CLEANED.markers, file = paste0(CCl4Output.path, date.1, "_aCCl4CLEANEDAnno_ALLmarkers.csv"))
#
## Load Data
amA5_cleaned <- readRDS(paste0(mA5Output.path, date.1, "_amA5_cleaned.RDS"))

aBleo2_cleaned <- readRDS(paste0(Bleo2Output.path, date.1, "_aBleo2_cleaned.RDS"))

aBDL2_cleaned <- readRDS(paste0(BDL2Output.path, date.1, "_aBDL2_cleaned.RDS"))

aCCl4_cleaned <- readRDS(paste0(CCl4Output.path, date.1, "_aCCl4_cleaned.RDS"))


##Make ALL Heatmaps (Supp) Cell-type

mA5CLEANED.markers <- read.csv(paste0(mA5Output.path, date.1, "_amA5CLEANEDAnno_ALLmarkers.csv"))
Bleo2CLEANED.markers <- read.csv(paste0(Bleo2Output.path, date.1, "_aBleo2CLEANEDAnno_ALLmarkers.csv"))
BDL2CLEANED.markers <- read.csv( paste0(BDL2Output.path, date.1, "_aBDL2CLEANEDAnno_ALLmarkers.csv"))
CCl4CLEANED.markers <- read.csv(paste0(CCl4Output.path, date.1, "_aCCl4CLEANEDAnno_ALLmarkers.csv"))

#Re_order CLEANED
#amA5_cleaned$Cell.type <- factor(amA5_cleaned$Cell.type, levels = c("Kidney macrophage 1", "Kidney macrophage 2", "Kidney macrophage 3", "Kidney macrophage 4", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney monocyte", "Blood Ly6Clo monocyte", "Healthy blood Ly6Chi monoyte", "UUO blood Ly6Chi monocyte","Basophil", "Blood B cell 1", "Kidney B cell", "Blood neutrophils", "Cd4+ T cell 1", "Cd4+ T cell 2", "CD8+ T cell 1", "Cd8+ T cell 2", "NK cell 1", "NK cell 2", "ILC2", "pDC", "CCR7+ cDC", "Kidney cDC1", "Kidney cDC2 1", "Kidney cDC2 2", "Platelet", "Proliferating", "Kidney endothelial cells"))
#
#aBleo2_cleaned$Cell.type <- factor(aBleo2_cleaned$Cell.type, levels = c( "Lung alveolar macrophage 1", "Lung alveolar macrophage 2", "Proliferating Lung Alveolar macrophage", "Lung interstitial macrophage 1", "Lung interstitial macrophage 2", "Lung macrophage 3", "Lung macrophage 4", "Lung Monocyte","Blood IFN-primed Ly6Chi Monocyte", "Blood Ly6Chi Monocyte", "Blood Ly6Clo Monocyte 1", "Blood Ly6Clo Monocyte 2", "Lcn2+ Monocyte",  "Basophil", "Blood B cell 1", "Blood B cell 2", "Lung B Cell 1", "Cd4+ T cell 1", "Blood CD4+ T cell 2", "Lung CD4+ T cell 2",  "CD8+ T cell 2", "Blood CD8+ T cell", "Lung CD8+ T cell", "Proliferating T cell", "ILC2", "NK cell 1", "CCR7+ cDC", "cDC1", "cDC2", "pDC", "Proliferating cDC", "Lung endothelial cells"))

#aBDL2_cleaned$Cell.type <- factor(aBDL2_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1",	"Liver Kupffer Cell 2",	"Liver macrophage 1",	"Liver macrophage 2",	"Liver macrophage 3",	"Liver macrophage 4",	"Liver macrophage 5",	"Blood Ly6Clo monocyte 1",	"Blood Ly6Clo monocyte 2",	"BDL Blood Ly6Clo monocyte 1",	"Healthy Blood Ly6Chi monocyte 1",	"Healthy Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Chi monocyte 1",	"BDL Blood Ly6Chi monocyte 2",	"Blood IFN-primed Ly6Chi monocyte",	"Proliferating",	"pDC",	"cDC1",	"cDC2",	"B cell 1",	"B cell 2",	"B cell 3",	"Cd4+ T cell 1",	"CD8+ T cell 1",	"T cell 2",	"NK cell 1",	"NK cell 2",	"Neutrophil 1",	"Neutrophil 2",	"Basophil",	"Platelet 1",	"Platelet 2",	"Liver Endothelial Cells"))

#aCCl4_cleaned$Cell.type <- factor(aCCl4_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Blood IFN-primed Ly6Chi monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","Healthy Blood Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Clo monocyte","Liver monocyte","B cell 1","B cell 2","B cell 3","Cd4+ T cell 1","CD8+ T cell 1","CD8+ T cell 2","Basophil","Eosinophil","Neutrophil","NK cell 1","NK cell 2","pDC","cDC1","cDC2","Liver CCR7+ cDC", "Proliferating","Liver Endothelial Cells"))

#Gene Heatmaps - Cell.type
#mA5CLEANED_top10 <- mA5CLEANED.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#mA5CLEANED_top3 <- mA5CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
#print(mA5CLEANED_top3$gene)
#mA5_cellHM <- DoHeatmap(amA5_cleaned, group.by = "Cell.type", features = mA5CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
#ggsave(filename = paste0(mA5Output.path, date.1, "_amA5CLEANED_Supp_celltypeHM.png"), plot = mA5_cellHM, width = 50, height = 32, units = "cm")
##-
#
#Bleo2CLEANED_top10 <- Bleo2CLEANED.markers %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_logFC)
#Bleo2CLEANED_top3 <- Bleo2CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
#print(Bleo2CLEANED_top3$gene)
#Bleo2_cellHM <- DoHeatmap(aBleo2_cleaned, group.by = "Cell.type", features = Bleo2CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
#ggsave(filename = paste0(Bleo2Output.path, date.1, "_Bleo2CLEANED_Supp_celltypeHM.png"), plot = Bleo2_cellHM,  width = 45, height = 32, units = "cm")
#
#-
#
#BDL2CLEANED_top10 <- BDL2CLEANED.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#BDL2CLEANED_top3 <- BDL2CLEANED_top10  %>% top_n(n = 3, wt = pct.2)
#print(BDL2CLEANED_top3$gene)
#BDL2_cellHM <- DoHeatmap(aBDL2_cleaned, features = BDL2CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
#ggsave(filename = paste0(BDL2Output.path, date.1, "_BDL2CLEANED_Supp_celltypeHM.png"), plot = BDL2_cellHM,  width = 45, height = 32, units = "cm")
#
###-
#CCl4CLEANED_top10 <- CCl4CLEANED.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#CCl4CLEANED_top3 <- CCl4CLEANED_top10%>% top_n(n = 3, wt = pct.2)
#print(CCl4CLEANED_top3$gene)
#CCl4_cellHM <- DoHeatmap(aCCl4_cleaned, group.by = "Cell.type", features = CCl4CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 3, angle = 90) + NoLegend()
#ggsave(filename = paste0(CCl4Output.path, date.1, "_CCl4CLEANED_Supp_celltypeHM.png"), plot = CCl4_cellHM,  width = 45, height = 32, units = "cm")


#Gene Heatmap -  lineage
#mA5CLEANED_top10 <- mA5CLEANED.markers %>% group_by(lineage) %>% top_n(n = 10, wt = avg_logFC)
#mA5CLEANED_top3 <- mA5CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
#print(mA5CLEANED_top3$gene)
#mA5_cellHM <- DoHeatmap(amA5_cleaned, group.by = "lineage", features = mA5CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 4, angle = 90, group.colors = mA5_lin.colpal) + NoLegend()
#ggsave(filename = paste0(mA5Output.path, date.1, "_amA5CLEANED_Supp_lineageHM.png"), plot = mA5_cellHM, width = 45, height = 18, units = "cm")
##-
#
#Bleo2CLEANED_top10 <- Bleo2CLEANED.markers %>% group_by(lineage)  %>% top_n(n = 10, wt = avg_logFC)
#Bleo2CLEANED_top3 <- Bleo2CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
#print(Bleo2CLEANED_top3$gene)
#Bleo2_cellHM <- DoHeatmap(aBleo2_cleaned, group.by = "lineage", features = Bleo2CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 4, angle = 90, group.colors = Bleo2_lin.colpal) + NoLegend()
#ggsave(filename = paste0(Bleo2Output.path, date.1, "_Bleo2CLEANED_Supp_lineageHM.png"), plot = Bleo2_cellHM, width = 45, height = 18, units = "cm")

#-

#BDL2CLEANED_top10 <- BDL2CLEANED.markers %>% group_by(lineage) %>% top_n(n = 10, wt = avg_logFC)
#BDL2CLEANED_top3 <- BDL2CLEANED_top10  %>% top_n(n = 3, wt = pct.2)
#print(BDL2CLEANED_top3$gene)
#BDL2_cellHM <- DoHeatmap(aBDL2_cleaned, group.by = "lineage", features = BDL2CLEANED_top3$gene, assay = "RNA", slot = "data", label = T, size = 4, angle = 90, group.colors = BDL2_lin.colpal) + NoLegend()
#ggsave(filename = paste0(BDL2Output.path, date.1, "_BDL2CLEANED_Supp_lineageHM.png"), plot = BDL2_cellHM, width = 45, height = 18, units = "cm")

#-
#CCl4CLEANED_top10 <- CCl4CLEANED.markers %>% group_by(lineage) %>% top_n(n = 10, wt = avg_logFC)
#CCl4CLEANED_top3 <- CCl4CLEANED_top10%>% top_n(n = 3, wt = pct.2)
#print(CCl4CLEANED_top3$gene)
#CCl4_cellHM <- DoHeatmap(aCCl4_cleaned, group.by = "lineage", features = CCl4CLEANED_top3$gene, assay = "RNA", slot = "data",  size = 4, angle = 90, group.colors = CCl4_lin.colpal) + NoLegend()
#ggsave(filename = paste0(CCl4Output.path, date.1, "_CCl4CLEANED_lineageHM.png"), plot = CCl4_cellHM, width = 45, height = 18, units = "cm")
#
#
##Re_order CLEANED
#amA5_cleaned$Cell.type <- factor(amA5_cleaned$Cell.type, levels = c("Kidney macrophage 1", "Kidney macrophage 2", "Kidney macrophage 3", "Kidney macrophage 4", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney monocyte", "Blood Ly6Clo monocyte", "Healthy blood Ly6Chi monoyte", "UUO blood Ly6Chi monocyte","Basophil", "Blood B cell 1", "Kidney B cell", "Blood neutrophils", "Cd4+ T cell 1", "Cd4+ T cell 2", "CD8+ T cell 1", "Cd8+ T cell 2", "NK cell 1", "NK cell 2", "ILC2", "pDC", "CCR7+ cDC", "Kidney cDC1", "Kidney cDC2 1", "Kidney cDC2 2", "Platelet", "Proliferating", "Kidney endothelial cells"))
#
#aBleo2_cleaned$Cell.type <- factor(aBleo2_cleaned$Cell.type, levels = c( "Lung alveolar macrophage 1", "Lung alveolar macrophage 2", "Proliferating Lung Alveolar macrophage", "Lung interstitial macrophage 1", "Lung interstitial macrophage 2", "Lung macrophage 3", "Lung macrophage 4", "Lung Monocyte","Blood IFN-primed Ly6Chi Monocyte", "Blood Ly6Chi Monocyte", "Blood Ly6Clo Monocyte 1", "Blood Ly6Clo Monocyte 2", "Lcn2+ Monocyte",  "Basophil", "Blood B cell 1", "Blood B cell 2", "Lung B Cell 1", "Cd4+ T cell 1", "Blood CD4+ T cell 2", "Lung CD4+ T cell 2",  "CD8+ T cell 2", "Blood CD8+ T cell", "Lung CD8+ T cell", "Proliferating T cell", "ILC2", "NK cell 1", "CCR7+ cDC", "cDC1", "cDC2", "pDC", "Proliferating cDC", "Lung endothelial cells"))

#aBDL2_cleaned$Cell.type <- factor(aBDL2_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1", "Liver Kupffer Cell 2", "Liver macrophage 1", "Liver macrophage 2", "Liver macrophage 3", "Liver macrophage 4", "Liver macrophage 5","BDL Blood Ly6Chi monocyte 1",	"BDL Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Clo monocyte 1", "Blood IFN-primed Ly6Chi monocyte", "Blood Ly6Clo monocyte 1", "Blood Ly6Clo monocyte 2", "Healthy Blood Ly6Chi monocyte 1", "Healthy Blood Ly6Chi monocyte 2", "B cell 1", "B cell 2", "B cell 3", "Basophil", "Cd4+ T cell 1", "CD8+ T cell 1", "T cell 2", "Neutrophil 1", "Neutrophil 2", "NK cell 1", "NK cell 2", "pDC", "cDC1", "cDC2", "Proliferating", "Platelet 1", "Platelet 2", "Liver Endothelial Cells"))

#aCCl4_cleaned$Cell.type <- factor(aCCl4_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Blood IFN-primed Ly6Chi monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","Healthy Blood Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Clo monocyte","Liver monocyte","B cell 1","B cell 2","B cell 3","Cd4+ T cell 1","CD8+ T cell 1","CD8+ T cell 2","Basophil","Eosinophil","Neutrophil","NK cell 1","NK cell 2","pDC","cDC1","cDC2","Liver CCR7+ cDC", "Proliferating","Liver Endothelial Cells"))
#
#
##Gene Heatmaps - MP
mA5CLEANED_top10 <- subset(mA5CLEANED.markers, subset = lineage == "MP") %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
mA5CLEANED_top3 <- mA5CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
print(mA5CLEANED_top3$gene)
mA5_MP <- subset(amA5_cleaned, subset = lineage =="MP")
mA5_MP$Cell.type <- factor(amA5_cleaned$Cell.type, levels = c("Kidney macrophage 1", "Kidney macrophage 2", "Kidney macrophage 3", "Kidney macrophage 4", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney monocyte", "Blood Ly6Clo monocyte", "Healthy blood Ly6Chi monoyte", "UUO blood Ly6Chi monocyte", "pDC", "CCR7+ cDC", "Kidney cDC1", "Kidney cDC2 1", "Kidney cDC2 2"))
print(unique(mA5_MP$Cell.type))
mA5_cellHM <- DoHeatmap(mA5_MP, group.by = "Cell.type", features = mA5CLEANED_top3$gene, size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(mA5Output.path, date.1, "_amA5CLEANED_Fig_MPHM.png"), plot = mA5_cellHM, width = 45, height = 30, units = "cm")
##-
#
Bleo2CLEANED_top10 <- subset(Bleo2CLEANED.markers, subset = lineage == "MP") %>% group_by(cluster)  %>% top_n(n = 10, wt = avg_logFC)
Bleo2CLEANED_top3 <- Bleo2CLEANED_top10 %>% top_n(n = 3, wt = pct.2)
print(Bleo2CLEANED_top3$gene)
Bleo2_MP <- subset(aBleo2_cleaned, subset = lineage =="MP")
aBleo2_cleaned$Cell.type <- factor(aBleo2_cleaned$Cell.type, levels = c( "Lung alveolar macrophage 1", "Lung alveolar macrophage 2", "Proliferating Lung Alveolar macrophage", "Lung interstitial macrophage 1", "Lung interstitial macrophage 2", "Lung macrophage 3", "Lung macrophage 4", "Lung Monocyte","Blood IFN-primed Ly6Chi Monocyte", "Blood Ly6Chi Monocyte", "Blood Ly6Clo Monocyte 1", "Blood Ly6Clo Monocyte 2", "Lcn2+ Monocyte",  "CCR7+ cDC", "cDC1", "cDC2", "pDC", "Proliferating cDC"))
print(unique(Bleo2_MP$Cell.type))
Bleo2_cellHM <- DoHeatmap(Bleo2_MP,  group.by = "Cell.type", features = Bleo2CLEANED_top3$gene, size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(Bleo2Output.path, date.1, "_Bleo2CLEANED_Fig_MPHM.png"), plot = Bleo2_cellHM, width = 45, units = "cm")
#
##-

BDL2CLEANED_top10 <- subset(BDL2CLEANED.markers, subset = lineage == "MP") %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
BDL2CLEANED_top3 <- BDL2CLEANED_top10  %>% top_n(n = 3, wt = pct.2)
print(BDL2CLEANED_top3$gene)
BDL2_MP <- subset(aBDL2_cleaned, subset = lineage == "MP")
BDL2_MP$Cell.type <- factor(BDL2_MP$Cell.type, levels = c("Liver Kupffer Cell 1", "Liver Kupffer Cell 2", "Liver macrophage 1", "Liver macrophage 2", "Liver macrophage 3", "Liver macrophage 4", "Liver macrophage 5","BDL Blood Ly6Chi monocyte 1",	"BDL Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Clo monocyte 1", "Blood IFN-primed Ly6Chi monocyte", "Blood Ly6Clo monocyte 1", "Blood Ly6Clo monocyte 2", "Healthy Blood Ly6Chi monocyte 1", "Healthy Blood Ly6Chi monocyte 2","pDC", "cDC1", "cDC2"))
print(unique(BDL2_MP$Cell.type))
Idents(BDL2_MP) <- "Cell.type"
BDL2_cellHM <- DoHeatmap(BDL2_MP,  group.by = "Cell.type", features = BDL2CLEANED_top3$gene, slot = "data", assay = "RNA", size = 3, angle = 90, label = T) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, date.1, "_BDL2CLEANED_Fig_MPHMtop3.png"), plot = BDL2_cellHM, width = 45, units = "cm")


BDL2CLEANED_top10 <- subset(BDL2CLEANED.markers, subset = lineage == "MP") %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
BDL2CLEANED_top5 <- BDL2CLEANED_top10  %>% top_n(n = 5, wt = pct.2)
print(BDL2CLEANED_top5$gene)
BDL2_MP <- subset(aBDL2_cleaned, subset = lineage == "MP")
BDL2_MP$Cell.type <- factor(BDL2_MP$Cell.type, levels = c("Liver Kupffer Cell 1", "Liver Kupffer Cell 2", "Liver macrophage 1", "Liver macrophage 2", "Liver macrophage 3", "Liver macrophage 4", "Liver macrophage 5","BDL Blood Ly6Chi monocyte 1",	"BDL Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Clo monocyte 1", "Blood IFN-primed Ly6Chi monocyte", "Blood Ly6Clo monocyte 1", "Blood Ly6Clo monocyte 2", "Healthy Blood Ly6Chi monocyte 1", "Healthy Blood Ly6Chi monocyte 2","pDC", "cDC1", "cDC2"))
print(unique(BDL2_MP$Cell.type))
Idents(BDL2_MP) <- "Cell.type"
BDL2_cellHM <- DoHeatmap(BDL2_MP,  group.by = "Cell.type", features = BDL2CLEANED_top5$gene, slot = "data", assay = "RNA", size = 3, angle = 90, label = T) + NoLegend()
ggsave(filename = paste0(BDL2Output.path, date.1, "_BDL2CLEANED_Fig_MPHMtop5.png"), plot = BDL2_cellHM, width = 50, height = 25, units = "cm")




#-
CCl4CLEANED_top10 <- subset(CCl4CLEANED.markers, subset = lineage == "MP") %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
CCl4CLEANED_top3 <- CCl4CLEANED_top10%>% top_n(n = 3, wt = pct.2)
print(CCl4CLEANED_top3$gene)
aCCl4_cleaned$Cell.type <- factor(aCCl4_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Blood IFN-primed Ly6Chi monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","Healthy Blood Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Clo monocyte","Liver monocyte","B cell 1","B cell 2","B cell 3","Cd4+ T cell 1","CD8+ T cell 1","CD8+ T cell 2","Basophil","Eosinophil","Neutrophil","NK cell 1","NK cell 2","pDC","cDC1","cDC2","Liver CCR7+ cDC", "Proliferating","Liver Endothelial Cells"))

CCl4_MP <- subset(aCCl4_cleaned, subset = lineage =="MP")
CCl4_MP$Cell.type <- factor(aCCl4_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Blood IFN-primed Ly6Chi monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","Healthy Blood Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Clo monocyte","Liver monocyte","pDC","cDC1","cDC2","Liver CCR7+ cDC"))
print(unique(CCl4_MP$Cell.type))
CCl4_cellHM <- DoHeatmap(CCl4_MP,  group.by = "Cell.type", features = CCl4CLEANED_top3$gene, slot = "data", assay = "RNA", size = 3, angle = 90) + NoLegend()
ggsave(filename = paste0(CCl4Output.path, date.1, "_CCl4CLEANED_Fig_MPHM.png"), plot = CCl4_cellHM, width = 45, units = "cm")
