#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Clusering, Annotation and Intergrate CellenOne Data             ##
## Data: 2025.03.20                                                                ##
## Author: Zhang Haoran                                                            ##
#####################################################################################

##----0.load package--------
library(tidyverse)
library(ggplot2)
library(iq)
library(circlize)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GseaVis)
library(ggbeeswarm)
library(cowplot)
library(protti)
library(Seurat)
library(scCustomize)
library(patchwork)
library(MSnID)
library(kBET)
library(scRNAtoolVis)
library(bbknnR)
library(SCP)
library(reticulate)
use_python("D:/Software/Miniconda3/envs/SCP")
library(sceasy)
library(UCell)
library(Nebulosa)
source("scp_utils.R")
#-------------------------------------------------------------------------------
##-----1. cell clustering--------
###----load data--------
tmp.path <- "res/r/SCP_MG_3186_seurat_zhr_step1_harmony_20250107_label_transferred_v3.rds"
seu <- readRDS(file = tmp.path)
Reductions(seu)

seu <- FindNeighbors(seu,
                     reduction = "harmony",
                     dims = 1:5)

seq <- seq(0.1,1,by = 0.1)
tmp.seu <- seu
for(res in seq){
  tmp.seu <- FindClusters(tmp.seu, resolution = res)
}
library(gglogger)
library(scplotter)
tmp.seu[[]]
p <- ClustreePlot(tmp.seu, prefix = "SCT_snn_res.")
p
myggsavePro(p = p,
            prefix = paste0("res/fig/SCP_MG_","SCP_",
                            "after_QC_",
                            "check_batch",
                            "_scvi_",
                            "_clustree"),
            width = 8,
            height = 10)

##-----2. cell annotation--------
DefaultAssay(tmp.seu) <- "SCT"
tmp.seu <- FindNeighbors(tmp.seu,
                     reduction = "harmony",
                     dims = 1:6)
tmp.seu <- FindClusters(tmp.seu,
                        resolution = 0.5,
                        cluster.name = "harmony_cluster",
                        print.output = TRUE)
tmp.seu@meta.data$harmony_cluster <- as.numeric(tmp.seu@meta.data$harmony_cluster)

p <- DimPlot(tmp.seu,label = T,
             repel = T,
             label.size = 9,
             group.by = "harmony_cluster",
             reduction = "harmonyUMAP",
             pt.size=1.5)+
  labs(title="0.5 res harmony Cluster",x="harmony UMAP 1",y="harmony UMAP 2")+
  scale_color_manual(values = c("#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                                "#EE934E","#D1352B","#9B5B33","#F5CFE4","#B383B9",
                                "#8FA4AE","#FCED82","#F5D2A8","#BBDD78","#f66c79"))+
  theme_cowplot(font_size = 25)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))+
  theme(legend.position = "none")

p

tmp.seu@meta.data$celltype <- as.numeric(as.character(tmp.seu@meta.data$harmony_cluster))
metadata <- tmp.seu[[]]
metadata[metadata$harmony_cluster == 1,][["celltype"]] = 'C1'
metadata[metadata$harmony_cluster == 2,][["celltype"]] = 'C2'
metadata[metadata$harmony_cluster == 3,][["celltype"]] = 'C3'
metadata[metadata$harmony_cluster == 4,][["celltype"]] = 'C1'
metadata[metadata$harmony_cluster == 5,][["celltype"]] = 'C5'
metadata[metadata$harmony_cluster == 6,][["celltype"]] = 'C4'
metadata[metadata$harmony_cluster == 7,][["celltype"]] = 'C6'
metadata[metadata$harmony_cluster == 8,][["celltype"]] = 'C2'
metadata[metadata$harmony_cluster == 9,][["celltype"]] = 'C1'
metadata[metadata$harmony_cluster == 10,][["celltype"]] = 'C1'
metadata[metadata$harmony_cluster == 11,][["celltype"]] = 'C5'
metadata[metadata$harmony_cluster == 12,][["celltype"]] = 'C4'
metadata[metadata$harmony_cluster == 13,][["celltype"]] = 'C7'
tmp.seu@meta.data <- metadata
Idents(tmp.seu) <- tmp.seu[[]]$celltype

p <- DimPlot(tmp.seu,label = T,
             repel = T,
             label.size = 9,
             group.by = "celltype",
             reduction = "harmonyUMAP",
             pt.size=1.5)+
  labs(title="0.5 res harmony Cluster",x="harmony UMAP 1",y="harmony UMAP 2")+
  scale_color_manual(values = c("#FCED82","#3C77AF","#7DBFA7","#AECDE1",
                                "#EE934E","#D1352B","#9B5B33","#F5CFE4","#B383B9",
                                "#8FA4AE","#D2EBC8","#F5D2A8","#BBDD78","#f66c79"))+
  theme_cowplot(font_size = 25)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))+
  theme(legend.position = "none")

p

##-----3.Integrate cell morphylogical data--------
excel_mor <- read.csv('cellenOne_MG_morphology.csv',
                      header = TRUE)
colnames(excel_mor)[1] <- "id"

tmp.seu <- tmp.seu[,!tmp.seu$loading_date == "20241014"]
tmp.seu <- tmp.seu[,!tmp.seu$nCount_SCP > 2000000]

metadata <- tmp.seu[[]] %>% 
  mutate(id=gsub("^.*_10min_","",rownames(tmp.seu[[]])))
setdiff(metadata$id,excel_mor$id)

metadata <- left_join(metadata,excel_mor,by="id")
head(metadata)
head(tmp.seu[[]])
tail(metadata)
tail(tmp.seu[[]])

tmp.seu[[]] <- metadata
SCP::FeatureDimPlot(tmp.seu,c("Diameter","Elongation","Circularity",
                              "Intensity","nFeature_SCP","nCount_SCP"),
                    pt.size=1,ncol=2,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="harmony UMAP 1",
                    ylab="harmony UMAP 2")

FeatureStatPlot(tmp.seu,stat.by = c("Diameter","Elongation","Circularity",
                                    "Intensity","nFeature_SCP","nCount_SCP"), 
                palcolor = c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                             "#EE934E","#a5917f","#F5CFE4","#634d7b",
                             "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                group.by = "celltype", plot_type = "box",
                comparisons = list(c("Homeostatic", "Aged_1"), 
                                   c("Homeostatic","Aged_2"),
                                   c("Homeostatic","Phagocytic")))

