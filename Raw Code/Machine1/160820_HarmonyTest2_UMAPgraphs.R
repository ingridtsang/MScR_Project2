# Title     : TODO
# Objective : TODO
# Created by: mjdioli
# Created on: 16/08/2020


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

mA5.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Kidney-UUO/"
Bleo2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Lung-Bleo/"
BDL2.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-BDL/"
CCl4.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Liver-CCl4/"
harmony.path <- "/home/mjdioli/Documents/Ingrid/Harmony_Test2_SSH/Merge/"

#Load data
harTest2_UMAP <- readRDS(paste0(harmony.path, "120820_HarmonyTest2_UMAP.RDS"))

harTest2_UMAP$Condition <-  factor(harTest2_UMAP$Condition, levels = c("Healthy", "Fibrotic"))


#plots Merged, Condition

har_UMAP_CondName <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "Condition", split.by = "Experiment")
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_CondExp.png"),  width = 50, height = 15, units = "cm")

har_UMAP_NameCond <- DimPlot(harTest2_UMAP, reduction = "umap", group.by = "Experiment", split.by = "Condition")
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_ExpCond.png"),  width = 45, height = 15, units = "cm")

har_UMAP_Cond <- DimPlot(harTest2_UMAP, reduction = "umap", split.by = "Condition", label = T)
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_Cond.png"), width = 40, units = "cm")

#plots Indie, Condition
UUO <- subset(x = harTest2_UMAP, subset = Experiment == "Kidney.UUO")
har_UMAP_UUO <- DimPlot(UUO, reduction = "umap", split.by = "Condition", label = T)
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_CondUUO.png"), width = 40, height = 15, units = "cm")

Bleo <- subset(x = harTest2_UMAP, subset = Experiment == "Lung.Bleo")
har_UMAP_Bleo <- DimPlot(Bleo,  reduction = "umap", split.by = "Condition", label = T)
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_CondBleo.png"), width = 40, height = 15, units = "cm")

BDL <- subset(x = harTest2_UMAP, subset = Experiment == "Liver.BDL")
har_UMAP_BDL <- DimPlot(BDL,  reduction = "umap", split.by = "Condition", label = T)
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_CondBDL.png"), width = 40, height = 15, units = "cm")

CCl4 <- subset(x = harTest2_UMAP, subset = Experiment == "Liver.CCl4")
har_UMAP_CCl4 <- DimPlot(CCl4,  reduction = "umap", split.by = "Condition", label = T)
ggsave(filename = paste0(harmony.path, "160820_Test2_UMAP_CondCCl4.png"), width = 40, height = 15, units = "cm")