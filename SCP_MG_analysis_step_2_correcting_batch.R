#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Correct batch effects                                            ##
## Data: 2025.03.20                                                                ##
## Author: Wang Longteng, Zhang Haoran                                             ##
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
seu <- readRDS("res/r/SCP_MG_3186_seurat_post_QC_428cutoff_zhr_20241228.rds")
tmp.seu <- seu
scp_count<-GetAssayData(tmp.seu, slot ="counts")
dim(scp_count)
non_empty_genes <- rowSums(scp_count)>0
tmp.seu<-tmp.seu[non_empty_genes,]
nrow(tmp.seu)
dim(GetAssayData(tmp.seu,slot ="counts"))
tmp.seu <- SCTransform(seu,assay = "SCP",
                       vars.to.regress = "nFeature_SCP")
VariableFeaturePlot(object = tmp.seu)

tmp.seu <- RunPCA(tmp.seu,npcs = 20)
SCP::CellDimPlot(tmp.seu,pt.size = 2,dims = c(1,2),
                 palcolor = c("#b2df8a","#a6cee3","#1f78b4"),
                 group.by = "chromatogram_infor",reduction = "pca")
VizDimLoadings(object = tmp.seu, dims = 9, reduction = "pca",nfeatures = 70)
DimHeatmap(object = tmp.seu, dims = 9, cells = 3000, 
           balanced = TRUE,nfeatures = 100,ncol=2)

pca_embeddings <- Embeddings(tmp.seu, reduction = "pca")
# remove PC1 and PC2
pca_embeddings <- pca_embeddings[, c(3:20)]

tmp.seu[["pca"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PC_", assay = DefaultAssay(tmp.seu))
library(harmony)
tmp.seu <- RunHarmony(tmp.seu, group.by.vars = "prepare_date",
                      reduction.use = "pca", dims.use = 1:18)

tmp.seu <- RunUMAP(tmp.seu,dims = c(1:6),assay = "SCT",
                   reduction = "harmony",
                   reduction.name = "harmonyUMAP")

tmp.seu[[]]$chromatogram_infor <- factor(tmp.seu[[]]$chromatogram_infor,
                                         levels = c("Io01","Io02","Io03"))
SCP::CellDimPlot(tmp.seu,pt.size = 0.5,palcolor = c("#b2df8a","#a6cee3","#1f78b4"),
                 group.by = "chromatogram_infor",reduction = "harmonyUMAP")
SCP::CellDimPlot(tmp.seu,pt.size = 0.5,palcolor = col.spectral(6)[4:6],
                 group.by = "age",reduction = "harmonyUMAP")

tmp.homeostatic <- c("P2ry12","Tgfbr1","Tmem119","Cx3cr1")
SCP::FeatureDimPlot(tmp.seu,tmp.homeostatic,pt.size=1,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="harmony UMAP 1",
                    ylab="harmony UMAP 2")


tmp.BAM <- c("Mrc1","Stab1",
             "Cbr2","Ptprc",
             "Csf1r",
             "Snx6",
             "Snx2",
             "Dab2","Ap1b1")
SCP::FeatureDimPlot(tmp.seu,tmp.BAM,pt.size=1,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="harmony UMAP 1",
                    ylab="harmony UMAP 2")


tmp.presyn <- c("Syp","Snap25","Syt1","Gphn","Syn1","Vamp2")
SCP::FeatureDimPlot(tmp.seu,tmp.presyn,pt.size=1,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="harmony UMAP 1",
                    ylab="harmony UMAP 2")

#-------------------check batch------------------------------------------------
tmp.cells <- WhichCells(tmp.seu,expression = age== "24M")
ttt <- tmp.seu[,tmp.cells]
table(ttt$brain_region,ttt$prepare_date)

tmp.data.plot <- Embeddings(ttt,reduction = "harmonyUMAP")
tmp.data.plot <- cbind(tmp.data.plot,ttt[[]])

p <- ggplot(tmp.data.plot,aes(harmonyUMAP_1,harmonyUMAP_2,
                              color=brain_region))+
  geom_point(size=2,stroke=0.5)+
  scale_color_manual(values =  c("#00758b","#9cca62"),
                     name="Brain Region")+
  labs(title="harmonyUMAP\n(2M male)",x="harmony UMAP 1",y="harmony UMAP 2")+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_text(size = 18))
p

SCP::CellDimPlot(ttt, group.by = "brain_region", 
            reduction = "harmonyUMAP",palcolor = c("#00758b","#9cca62"),pt.size=2)
 # save as 6*7, landscape
p <- ggplot(tmp.data.plot,aes(harmonyUMAP_1,harmonyUMAP_2,
                              fill=prepare_date,
                              shape=brain_region))+
  geom_point(size=3,stroke=0.5,alpha=0.5)+
  scale_fill_manual(values =  c("#5bb1cb","#f7962e"),
                    guide = guide_legend(override.aes = list(shape = 21)),
                    label=c("Batch 1","Batch 2"))+
  scale_shape_manual(values = c(21,24),name="Brain Region")+
  labs(title="harmonyUMAP\n(2M male)",x="harmony UMAP 1",y="harmony UMAP 2")+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_text(size = 18))
p
SCP::CellDimPlot(ttt, group.by = "prepare_date", 
                 reduction = "harmonyUMAP",palcolor = c("#5bb1cb","#f7962e"),pt.size=2)
# save as 6*7, landscape

#-----------------------label transfering---------------------------------------
library(qs)
rna <- qread("hammond_seurat/Hammond_et-al-2019_Seurat_Converted_v4.qs")
rna <- subset(rna, subset = Sex == "Male")
rna <- SCTransform(rna)
DefaultAssay(tmp.seu) <- "SCT"
scp_seu <- tmp.seu

markers <- read.csv("markers_rna.csv", header=FALSE)
markers <- as.vector(markers)
markers <- markers$V1
markers <- markers[markers %in% rownames(rna)]
markers <- markers[markers %in% rownames(scp_seu)]

rna_2M <- subset(rna, Age %in% c("P30", "P100"))
rna_14M <- subset(rna, Age %in% c("P100", "P540"))
rna_24M <- subset(rna, Age == "P540")

#    predict ID----------------------------------------------------------------------

scp_2M_0815 <- subset(scp_seu, prepare_date == "0815" & age == "2M")
scp_2M_1015 <- subset(scp_seu, prepare_date == "1015" & age == "2M")
scp_14M_0911 <- subset(scp_seu, prepare_date == "0911" & age == "14M")
scp_14M_0913 <- subset(scp_seu, prepare_date == "0913" & age == "14M")
scp_24M_0822 <- subset(scp_seu, prepare_date == "0822" & age == "24M")
scp_24M_1105 <- subset(scp_seu, prepare_date == "1105" & age == "24M")

anchors_2M_0815 <- FindTransferAnchors(reference = rna_2M,  
                                       query = scp_2M_0815, 
                                       features = markers,
                                       normalization.method = "SCT",
                                       recompute.residuals = FALSE)

anchors_2M_1015 <- FindTransferAnchors(reference = rna_2M,  
                                       query = scp_2M_1015, 
                                       features = markers,
                                       normalization.method = "SCT",
                                       recompute.residuals = FALSE)

anchors_14M_0911 <- FindTransferAnchors(reference = rna_14M,  
                                        query = scp_14M_0911, 
                                        features = markers,
                                        normalization.method = "SCT",
                                        recompute.residuals = FALSE)

anchors_14M_0913 <- FindTransferAnchors(reference = rna_14M,  
                                        query = scp_14M_0913, 
                                        features = markers,
                                        normalization.method = "SCT",
                                        recompute.residuals = FALSE)

anchors_24M_0822 <- FindTransferAnchors(reference = rna_24M,  
                                        query = scp_24M_0822, 
                                        features = markers,
                                        normalization.method = "SCT",
                                        recompute.residuals = FALSE)

anchors_24M_1105 <- FindTransferAnchors(reference = rna_24M,  
                                        query = scp_24M_1105, 
                                        features = markers,
                                        normalization.method = "SCT",
                                        recompute.residuals = FALSE)




prediction_scp_2M_0815 <- TransferData(anchorset = anchors_2M_0815,  
                                       refdata = rna_2M$Paper_Cluster)

scp_2M_0815 <- AddMetaData(scp_2M_0815, metadata = prediction_scp_2M_0815)

prediction_scp_2M_1015 <- TransferData(anchorset = anchors_2M_1015,  
                                       refdata = rna_2M$Paper_Cluster)

scp_2M_1015 <- AddMetaData(scp_2M_1015, metadata = prediction_scp_2M_1015)


prediction_scp_14M_0911 <- TransferData(anchorset = anchors_14M_0911,  
                                        refdata = rna_14M$Paper_Cluster)

scp_14M_0911 <- AddMetaData(scp_14M_0911, metadata = prediction_scp_14M_0911)

prediction_scp_14M_0913 <- TransferData(anchorset = anchors_14M_0913,  
                                        refdata = rna_14M$Paper_Cluster)

scp_14M_0913 <- AddMetaData(scp_14M_0913, metadata = prediction_scp_14M_0913)


prediction_scp_24M_0822 <- TransferData(anchorset = anchors_24M_0822,  
                                        refdata = rna_24M$Paper_Cluster)

scp_24M_0822 <- AddMetaData(scp_24M_0822, metadata = prediction_scp_24M_0822)

prediction_scp_24M_1105 <- TransferData(anchorset = anchors_24M_1105,  
                                        refdata = rna_24M$Paper_Cluster)

scp_24M_1105 <- AddMetaData(scp_24M_1105, metadata = prediction_scp_24M_1105)


scp_seu$predicted.id <- NA
scp_seu$predicted.id[colnames(scp_2M_0815)] <- scp_2M_0815$predicted.id
scp_seu$predicted.id[colnames(scp_2M_1015)] <- scp_2M_1015$predicted.id
scp_seu$predicted.id[colnames(scp_14M_0911)] <- scp_14M_0911$predicted.id
scp_seu$predicted.id[colnames(scp_14M_0913)] <- scp_14M_0913$predicted.id
scp_seu$predicted.id[colnames(scp_24M_0822)] <- scp_24M_0822$predicted.id
scp_seu$predicted.id[colnames(scp_24M_1105)] <- scp_24M_1105$predicted.id


p <- DimPlot(scp_seu,
             reduction = "harmonyUMAP",
             group.by = "predicted.id",
             label = F,
             label.size = 6,
             repel = TRUE,
             cols = c("#EE934E","#D1352B","#9B5B33","#D2EBC8",
                      "#3C77AF","#7DBFA7","#AECDE1",
                      "#F5CFE4","#FCED82","#BBDD78"),
             pt.size=2) +
  ggtitle("Query transferred labels (SCP)") +
  labs(x="harmony UMAP 1",y="harmony UMAP 2")
p

scp_seu[[]]$predicted.id <- factor(scp_seu[[]]$predicted.id,
                                   levels = c("2a","2b","5","6","7a","7b","7c","8","9","Mono/Mac"))
SCP::CellDimPlot(scp_seu, group.by = "predicted.id", reduction = "harmonyUMAP",
                 theme_use = "theme_blank",pt.size = 1,
                 palcolor=c("#EE934E","#D1352B","#9B5B33","#D2EBC8",
                            "#3C77AF","#7DBFA7","#AECDE1",
                            "#F5CFE4","#FCED82","#BBDD78"))

scplotter::CellStatPlot(scp_seu, plot_type = "sankey", 
                        ident = "celltype",
                        group_by = c("celltype", "predicted.id"),
                        
                        links_alpha = .6)


scp_seu[[]]
scplotter::CellStatPlot(scp_seu, group_by = "age", 
                        ident = "predicted.id",frac = "group",
                        ylab="Proportions (%)",
                        alpha=0.8,
                        palcolor=c("#EE934E","#D1352B","#9B5B33","#D2EBC8",
                                   "#3C77AF","#7DBFA7","#AECDE1",
                                   "#F5CFE4","#FCED82","#BBDD78"),
                        swap = TRUE, position = "stack")
#save as 4*6 inches, portait

#    predict Age----------------------------------------------------------------------

rna[[]]

prediction_scp_2M_0815 <- TransferData(anchorset = anchors_2M_0815,  
                                       refdata = rna_2M$Age)

scp_2M_0815 <- AddMetaData(scp_2M_0815, metadata = prediction_scp_2M_0815)

prediction_scp_2M_1015 <- TransferData(anchorset = anchors_2M_1015,  
                                       refdata = rna_2M$Age)

scp_2M_1015 <- AddMetaData(scp_2M_1015, metadata = prediction_scp_2M_1015)


prediction_scp_14M_0911 <- TransferData(anchorset = anchors_14M_0911,  
                                        refdata = rna_14M$Age)

scp_14M_0911 <- AddMetaData(scp_14M_0911, metadata = prediction_scp_14M_0911)

prediction_scp_14M_0913 <- TransferData(anchorset = anchors_14M_0913,  
                                        refdata = rna_14M$Age)

scp_14M_0913 <- AddMetaData(scp_14M_0913, metadata = prediction_scp_14M_0913)


prediction_scp_24M_0822 <- TransferData(anchorset = anchors_24M_0822,  
                                        refdata = rna_24M$Age)

scp_24M_0822 <- AddMetaData(scp_24M_0822, metadata = prediction_scp_24M_0822)

prediction_scp_24M_1105 <- TransferData(anchorset = anchors_24M_1105,  
                                        refdata = rna_24M$Age)

scp_24M_1105 <- AddMetaData(scp_24M_1105, metadata = prediction_scp_24M_1105)



scp_seu$predicted.id <- NA
scp_seu$predicted.id[colnames(scp_2M_0815)] <- scp_2M_0815$predicted.id
scp_seu$predicted.id[colnames(scp_2M_1015)] <- scp_2M_1015$predicted.id
scp_seu$predicted.id[colnames(scp_14M_0911)] <- scp_14M_0911$predicted.id
scp_seu$predicted.id[colnames(scp_14M_0913)] <- scp_14M_0913$predicted.id
scp_seu$predicted.id[colnames(scp_24M_0822)] <- scp_24M_0822$predicted.id
scp_seu$predicted.id[colnames(scp_24M_1105)] <- scp_24M_1105$predicted.id


#save as 4*6 inches, portait

scp_seu[[]]$predicted.id <- factor(scp_seu[[]]$predicted.id,
                                   levels = c("P30","P100","P540"))
scplotter::CellStatPlot(scp_seu, group_by = "age", 
                        ident = "predicted.id",frac = "group",
                        palcolor=col.spectral(6)[4:6],
                        alpha = 0.8,ylab="Proportions (%)",
                        swap = TRUE, position = "stack")
#save as 4*6 inches, portait

saveRDS(scp_seu,
        "res/r/SCP_MG_3184_seurat_zhr_step1_harmony_20241230_label_transferred.rds")
