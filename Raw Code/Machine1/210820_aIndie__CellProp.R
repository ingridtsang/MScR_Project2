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

#Set colour
mA5.colpal <- c("#CA9858",	"#FF7F00",	"#C8828D",	"#FFCA62")

Bleo2.colpal <- c("#9683B0",	"#6A3D9A",	"#C8828D",	"#DA98D6")

BDL2.colpal <- c("#8EB26E",	"#88BCB2",	"#33A02C",	"#00CDB3",	"#C8828D",	"#F4A1AF",	"#B2DF8A",	"#AFF1E1")

CCl4.colpal <- c("#74909E",	"#167EDB",	"#C8828D",	"#A6CEE3")

#Load data
amA5_cleaned <- readRDS(paste0(mA5Output.path, date.1, "_amA5_cleaned.RDS"))

aBleo2_cleaned <- readRDS(paste0(Bleo2Output.path, date.1, "_aBleo2_cleaned.RDS"))

aBDL2_cleaned <- readRDS(paste0(BDL2Output.path, date.1, "_aBDL2_cleaned.RDS"))

aCCl4_cleaned <- readRDS(paste0(CCl4Output.path, date.1, "_aCCl4_cleaned.RDS"))

# Individual Visualisations CELL -> Proportion
### Dataset membership vary by cell type ###

mA5 <- subset(amA5_cleaned, subset = lineage == "MP")
mA5$Cell.type <- factor(mA5$Cell.type, levels = c("Kidney macrophage 1", "Kidney macrophage 2", "Kidney macrophage 3", "Kidney macrophage 4", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney monocyte", "Blood Ly6Clo monocyte", "Healthy blood Ly6Chi monoyte", "UUO blood Ly6Chi monocyte", "pDC", "CCR7+ cDC", "Kidney cDC1", "Kidney cDC2 1", "Kidney cDC2 2"))
mA5$Names <- factor(mA5$Names, levels = c("UUO-Kidney", "H-Kidney", "UUO-Blood", "H-Blood"))

Idents(mA5) <- "Cell.type"
table(Idents(mA5))
mA5_CellNameProp <- prop.table(table(Idents(mA5), mA5$Names))
mA5_CellNameProp
mA5_DataProp <- as.data.frame(mA5_CellNameProp)

names(mA5_DataProp)[names(mA5_DataProp) == "Var1"] <- "Cell_type"
names(mA5_DataProp)[names(mA5_DataProp) == "Var2"] <- "Dataset"
names(mA5_DataProp)[names(mA5_DataProp) == "Freq"] <- "Proportion"

plot_mA5_DataPropBAR <- ggplot(mA5_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = mA5.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(mA5Output.path, date.2, "_amA5CLEANED_F_MP_CellPropBAR_Cell.png"), plot = plot_mA5_DataPropBAR, width = 30, units = "cm")

plot_mA5_DataPropBAR_Dataset <- ggplot(mA5_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(mA5Output.path, date.2, "_mA5CLEANED_F_MP_CellPropBAR_Data.png"), plot = plot_mA5_DataPropBAR_Dataset, width = 30, units = "cm")

#-
Bleo2 <- subset(aBleo2_cleaned, subset = lineage == "MP")
Bleo2$Cell.type <- factor(Bleo2$Cell.type, levels = c( "Lung alveolar macrophage 1", "Lung alveolar macrophage 2", "Proliferating Lung Alveolar macrophage", "Lung interstitial macrophage 1", "Lung interstitial macrophage 2", "Lung macrophage 3", "Lung macrophage 4", "Lung Monocyte","Blood IFN-primed Ly6Chi Monocyte", "Blood Ly6Chi Monocyte", "Blood Ly6Clo Monocyte 1", "Blood Ly6Clo Monocyte 2", "Lcn2+ Monocyte",  "CCR7+ cDC", "cDC1", "cDC2", "pDC", "Proliferating cDC"))
Bleo2$Names <- factor(Bleo2$Names, levels = c("Bleo-Lung", "H-Lung", "Bleo-Blood", "H-Blood"))

Idents(Bleo2) <- "Cell.type"
table(Idents(Bleo2))
Bleo2_CellNameProp <- prop.table(table(Idents(Bleo2), Bleo2$Names))
Bleo2_CellNameProp
Bleo2_DataProp <- as.data.frame(Bleo2_CellNameProp)

names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Var1"] <- "Cell_type"
names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Var2"] <- "Dataset"
names(Bleo2_DataProp)[names(Bleo2_DataProp) == "Freq"] <- "Proportion"


plot_Bleo2_DataPropBAR <- ggplot(Bleo2_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = Bleo2.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(Bleo2Output.path, date.2, "_Bleo2CLEANED_DataPropBAR_Celltype.png"), plot = plot_Bleo2_DataPropBAR, width = 30, units = "cm")

plot_Bleo2_DataPropBAR_Dataset <- ggplot(Bleo2_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(Bleo2Output.path, date.2, "_Bleo2CLEANED_DataPropBAR_Dataset.png"), plot = plot_Bleo2_DataPropBAR_Dataset, width = 30, units = "cm")

#-
BDL2 <- subset(aBDL2_cleaned, subset = lineage == "MP")
BDL2$Cell.type <- factor(BDL2$Cell.type, levels = c("Liver Kupffer Cell 1", "Liver Kupffer Cell 2", "Liver macrophage 1", "Liver macrophage 2", "Liver macrophage 3", "Liver macrophage 4", "Liver macrophage 5","BDL Blood Ly6Chi monocyte 1",	"BDL Blood Ly6Chi monocyte 2",	"BDL Blood Ly6Clo monocyte 1", "Blood IFN-primed Ly6Chi monocyte", "Blood Ly6Clo monocyte 1", "Blood Ly6Clo monocyte 2", "Healthy Blood Ly6Chi monocyte 1", "Healthy Blood Ly6Chi monocyte 2","pDC", "cDC1", "cDC2"))
BDL2$Names <- factor(BDL2$Names, levels = c("BDL-Liver-A", "BDL-Liver-B", "H-Liver-A", "H-Liver-B", "BDL-Blood-A", "BDL-Blood-B", "H-Blood-A", "H-Blood-B"))


Idents(BDL2) <- "Cell.type"
table(Idents(BDL2))
BDL2_CellNameProp <- prop.table(table(Idents(BDL2), BDL2$Names))
BDL2_CellNameProp
BDL2_DataProp <- as.data.frame(BDL2_CellNameProp)

names(BDL2_DataProp)[names(BDL2_DataProp) == "Var1"] <- "Cell_type"
names(BDL2_DataProp)[names(BDL2_DataProp) == "Var2"] <- "Dataset"
names(BDL2_DataProp)[names(BDL2_DataProp) == "Freq"] <- "Proportion"

plot_BDL2_DataPropBAR <- ggplot(BDL2_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = BDL2.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(BDL2Output.path, date.2, "_BDL2CLEANED_DataPropBAR_Celltype.png"), plot = plot_BDL2_DataPropBAR, width = 30, units = "cm")

plot_BDL2_DataPropBAR_Dataset <- ggplot(BDL2_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(BDL2Output.path, date.2, "_BDL2CLEANED_DataPropBAR_Dataset.png"), plot = plot_BDL2_DataPropBAR_Dataset, width = 30, units = "cm")


#-
CCl4 <- subset(aCCl4_cleaned, subset = lineage == "MP")
CCl4$Cell.type <- factor(aCCl4_cleaned$Cell.type, levels = c("Liver Kupffer Cell 1","CCl4 Liver Kupffer Cell 2","Liver macrophage 1","Blood IFN-primed Ly6Chi monocyte","CCl4 Blood Ly6Chi monocyte 1","CCl4 Blood Ly6Chi monocyte 2","CCl4 Blood Ly6Chi monocyte 3","Healthy Blood Ly6Chi monocyte","Blood Ly6Clo monocyte","CCl4 Blood Ly6Clo monocyte","Liver monocyte","pDC","cDC1","cDC2","Liver CCR7+ cDC"))
CCl4$Names <- factor(CCl4$Names, levels = c("CCl4-Liver", "H-Liver", "CCl4-Blood", "H-Blood"))

Idents(CCl4) <- "Cell.type"
table(Idents(CCl4))
CCl4_CellNameProp <- prop.table(table(Idents(CCl4), CCl4$Names))
CCl4_CellNameProp
CCl4_DataProp <- as.data.frame(CCl4_CellNameProp)


names(CCl4_DataProp)[names(CCl4_DataProp) == "Var1"] <- "Cell_type"
names(CCl4_DataProp)[names(CCl4_DataProp) == "Var2"] <- "Dataset"
names(CCl4_DataProp)[names(CCl4_DataProp) == "Freq"] <- "Proportion"

plot_CCl4_DataPropBAR <- ggplot(CCl4_DataProp, aes(fill=Dataset, y=Proportion, x=Cell_type)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = CCl4.colpal) + theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(CCl4Output.path, date.2, "_CCl4CLEANED_DataPropBAR_Celltype.png"), plot = plot_CCl4_DataPropBAR, width = 30, units = "cm")

plot_CCl4_DataPropBAR_Dataset <- ggplot(CCl4_DataProp, aes(fill=Cell_type, y=Proportion, x=Dataset)) + geom_bar(position="fill", stat="identity")
ggsave(filename = paste0(CCl4Output.path, date.2, "_CCl4CLEANED_DataPropBAR_Dataset.png"), plot = plot_CCl4_DataPropBAR_Dataset, width = 30, units = "cm")
