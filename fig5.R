#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Data Visualization                                              ##
## Data: 2025.03.31                                                                ##
## Author: Haoran Zhang                                                            ##
#####################################################################################

##----0.load package--------
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(Seurat)
library(SCP)
library(DEP)
use_python("D:/Software/Miniconda3/envs/SCP")
library(UCell)
source("scp_utils.R")

##----fig 5A, D
scplotter::CellDimPlot(tmp.seu,
            group_by = "celltype", facet_by = "celltype", reduction = "harmonyUMAP",
            highlight = 'celltype == "3"', theme = "theme_blank", 
            legend.position = "none",highlight_stroke=0,highlight_size=0.5,
            palcolor=c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                       "#EE934E","#a5917f","#F5CFE4","#634d7b",
                       "#8FA4AE","#F5D2A8","#FCED82","#BBDD78")
)
##----fig 5B
SCP::FeatureDimPlot(tmp.seu,c("Mtnr1a","Eif5a"),pt.size=0.5,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="UMAP 1",
                    ylab="UMAP 2")
##----fig 5C
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
DefaultAssay(ttt) <- "SCP"
p <- DotPlot(ttt,
        features = c("Creb1","Mtnr1a","Gsn","Jakmip2","Ccdc88b",
                     "Dsg1c","Dsg1b","Pgm2l1","Nckap1","RO60","Aloxe3","Alox12b",
                     "Cdsn","Rps27a","Cox16","Lyz1","Hdgfl3","Ddx3y","Agl",
                     "Ppp2r1a","Eif5a","Gsdma","Prss3b","Dpp7","Dst","Plbd1"),
        group.by = "celltype",
        scale = T,
        cols =my_col)+
  RotatedAxis()+
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral")[11:1])(100))+
  theme(legend.position = "top")
p

genes <- p$data %>% subset(id=="3") %>% arrange(desc(pct.exp))
genes <- genes$features.plot
DotPlot(ttt,
        features = genes,
        group.by = "celltype",
        scale = T,
        cols =my_col)+
  RotatedAxis()+
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral")[11:1])(100))+
  theme(legend.position = "top")

##----fig 5E
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
ht1 <- SCP::GroupHeatmap(ttt,
                    features = c(
                      "Arg1","Lamp1","Wdr77","Anxa2","Anxa7","Anxa11","Lgalsl",
                      "Lgals3","Vim","Pld2","Gba1",
                      "Aloxe3","Alox12b","Tgm1","Tgm3",
                      "Prdx2","Prdx4","Flna","Flnb","Got1"),border=F,
                    group_palcolor=c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                                      "#EE934E","#a5917f","#F5CFE4","#634d7b",
                                      "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                    group.by = c("celltype")
)
ht1$plot
##----fig 5F
SCP::FeatureDimPlot(tmp.seu,c("Wdr77","Lamp1"),pt.size=0.5,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="UMAP 1",
                    ylab="UMAP 2")
##----fig 5G
library(AUCell)
library(msigdbr)
msigdbr(species = "Mus musculus") %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) %>% 
  print(n=100)

GOMF_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:MF") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
str(GOMF_list)

peptidase <- GOMF_list$GOMF_PEPTIDASE_ACTIVITY
tmp.seu <- AddModuleScore(tmp.seu,
                          features = peptidase,
                          ctrl = 100,
                          name = "Peptidase_Score")
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
SCP::FeatureDimPlot(tmp.seu, 
                    features = c("Peptidase_Score1"),
                    reduction = "harmonyUMAP",pt.size=0.5,ncol = 4,
                    palcolor = sc.hic.orange(100))

##----fig 5H
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
scplotter::FeatureStatPlot(ttt, features = c("Peptidase_Score1"), 
                ident = "celltype", plot_type = "box",
                palcolor = c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                             "#EE934E","#a5917f","#F5CFE4","#634d7b",
                             "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                comparisons = list(c("1", "2"),c("6","5") ,
                                   c("4","5"),
                                   c("3","6"),
                                   c("2","6")))
##----fig 5I
p <- DotPlot(ttt,
        features = c("Arg1","Lgals3","Apoe","Hexb","P2ry12","Fscn1",
                     "Abi3","Bin1","Hpgds","Basp1","Ckb","Rap1gds1",
                     "Crybb1","Golm1","Sall1","F11r","Itga6","Csf1r",
                     "Pde3b","Abcc3","Cst3","Cx3cr1","Tgfbr1","Tmem119"
                     ),
        group.by = "celltype",
        scale = T)+
   RotatedAxis()+
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral")[11:1])(100))+
  theme(legend.position = "top")
p

