# Load library
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

merge.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/"
Experiment.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Merge/RegressExp/"

IndieOutput.path <-  "/home/mjdioli/Documents/Ingrid/Harmony_Test3/Graphs/R_Indie/"

#Set date

Load.date <- "070820"
Save.date <- "160820"
date.1 <- "200820"
date.2 <- "210820"
date.3 <- "220820"
date.4 <- "240820"
date.5 <- "280820"

#Set colour
mA5Names.colpal <- c("#CA9858",	"#FF7F00",	"#C8828D",	"#FFCA62")

Bleo2Names.colpal <- c("#9683B0",	"#6A3D9A",	"#C8828D",	"#DA98D6")

BDL2Names.colpal <- c("#8EB26E",	"#88BCB2",	"#33A02C",	"#00CDB3",	"#C8828D",	"#F4A1AF",	"#B2DF8A",	"#AFF1E1")

CCl4Names.colpal <- c("#74909E",	"#167EDB",	"#C8828D",	"#A6CEE3")


mA5HM.colpal <- c("#4c1d00",	"#ff6100",	"#ffcfb2")

Bleo2HM.colpal <- c("#9683B0",	"#6A3D9A",	"#C8828D",	"#DA98D6")

BDL2HM.colpal <- c("#8EB26E",	"#88BCB2",	"#33A02C",	"#00CDB3",	"#C8828D",	"#F4A1AF",	"#B2DF8A",	"#AFF1E1")

CCl4HM.colpal <- c("#74909E",	"#167EDB",	"#C8828D",	"#A6CEE3")


d <- c("mA5") #, "Bleo2", "BDL2", "CCl4")



mA5CellSeq <- c("Kidney macrophage 1", "Kidney macrophage 2", "Kidney macrophage 3", "Kidney macrophage 4", "Kidney proliferating macrophage 1", "Kidney proliferating macrophage 2", "Kidney IFN-primed macrophage", "Kidney monocyte", "Blood Ly6Clo monocyte", "Healthy blood Ly6Chi monoyte", "UUO blood Ly6Chi monocyte", "pDC", "CCR7+ cDC", "Kidney cDC1", "Kidney cDC2 1", "Kidney cDC2 2")



dataset <- d
for (d in dataset) {

  path <- paste0("/home/mjdioli/Documents/Ingrid/Harmony_Test3/", d, "/")
  HM <- readRDS(paste0(path, date.1, "_a", d, "_cleaned.RDS"))
  Markers <- read.csv(file = paste0(path, date.1, "_a", d, "CLEANEDAnno_ALLmarkers.csv"))

  top10 <- subset(Markers, subset = lineage == "MP") %>% group_by(cluster) %>% top_n(n = -10, wt = p_val_adj)
  top3 <- top10 %>% top_n(n = 3, wt = pct.2)
  print(top3$gene)
  MP <- subset(HM, subset = lineage =="MP")

  if_else(d =="mA5", l == mA5CellSeq,
          if_else(d == "BDL2", l == BDL2CellSeq,
                  if_else(d == "Bleo2", l == Bleo2CellSeq,
                          if_else(d == "CCl4", l == CCl4CellSeq, print("bug")))))


  HM$Cell.type <- factor(HM$Cell.type, levels = l)
  print(unique(HM$Cell.type))

  plot_top10 <- DoHeatmap(HM, group.by = "Cell.type", features = top10$gene, size = 3, angle = 90)
  plot_top3 <- DoHeatmap(HM, group.by = "Cell.type", features = top3$gene, size = 3, angle = 90)

  if_else(d =="mA5", c == mA5HM.colpal,
          if_else(d == "BDL2", c == BDL2HM.colpal,
                  if_else(d == "Bleo2", c == Bleo2HM.colpal,
                          if_else(d == "CCl4", c == CCl4HM.colpal, print("bug")))))

  plot_top10 <- plot_top10 + scale_fill_gradient2(c)
  plot_top3 <- plot_top3 + scale_fill_gradient2(c)
  ggsave(filename = paste0(IndieOutput.path, date.5,"_", d, "HM_MP10.svg"), plot = plot_top10, width = 45, height = 30, units = "cm")

  ggsave(filename = paste0(IndieOutput.path, date.5,"_", d, "HM_MP3.svg"), plot = plot_top3, width = 45, height = 30, units = "cm")
}
