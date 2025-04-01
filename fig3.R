#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Data Visualization                                              ##
## Data: 2025.03.31                                                                ##
## Author: Zhang Haoran                                                            ##
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

##----fig 3A
tmp.seu[[]]$celltype <- factor(tmp.seu[[]]$celltype,
                                 levels = c("1",
                                            "2","3","4","5","6","7"))

SCP::CellDimPlot(tmp.seu, group.by = "celltype", reduction = "harmonyUMAP",
                 theme_use = "theme_blank",pt.size = 1,
                 palcolor=c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                            "#EE934E","#a5917f","#F5CFE4","#634d7b",
                            "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"))
##----fig 3B
tmp.cells <- WhichCells(tmp.seu,expression = age== "2M")
ttt <- tmp.seu[,tmp.cells]
table(ttt$brain_region,ttt$prepare_date)
SCP::CellDimPlot(ttt, group.by = "brain_region", 
            reduction = "harmonyUMAP",palcolor = c("#00758b","#9cca62"),pt.size=2)
 # save as 6*7, landscape

tmp.cells <- WhichCells(tmp.seu,expression = age== "14M")
ttt <- tmp.seu[,tmp.cells]
table(ttt$brain_region,ttt$prepare_date)
SCP::CellDimPlot(ttt, group.by = "brain_region", 
            reduction = "harmonyUMAP",palcolor = c("#00758b","#9cca62"),pt.size=2)
 # save as 6*7, landscape

tmp.cells <- WhichCells(tmp.seu,expression = age== "24M")
ttt <- tmp.seu[,tmp.cells]
table(ttt$brain_region,ttt$prepare_date)
SCP::CellDimPlot(ttt, group.by = "brain_region", 
            reduction = "harmonyUMAP",palcolor = c("#00758b","#9cca62"),pt.size=2)
 # save as 6*7, landscape

##----fig 3C
tmp.data.plot <- tmp.seu[[]] %>%
  group_by(celltype,brain_region) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(brain_region) %>%
  mutate(brain_region_sum = sum(n)) %>%
  ungroup() %>%
  mutate(brain_region_pct = n/brain_region_sum) %>%
  group_by(celltype) %>%
  mutate(pct = brain_region_pct/sum(brain_region_pct)) %>%
  ungroup()
  # group_by(celltype,brain_region) %>%
  # summarise(n = n()) %>%
  # mutate(pct = n/sum(n)) %>%
  # ungroup()


p <- ggplot(tmp.data.plot,
            aes(celltype,pct,fill=brain_region))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Cell cluster")+
  #scale_x_discrete(label = c("Hippocampus","Prefrontal cortex"))+
  scale_fill_manual(values = c("#00758b","#9cca62"),
                    name="Brain Region")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(#angle = 45,
                                 hjust = 1,
                                 vjust = 1))
p# save as 4*6, landscape


##----fig 3D
###Cluster Propotion
tmp.data.plot <- tmp.seu[[]] %>%
  group_by(celltype,age) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(age) %>% 
  mutate(age_sum = sum(n)) %>% 
  ungroup() %>% 
  mutate(age_pct = n/age_sum) %>% 
  group_by(celltype) %>%
  mutate(pct = age_pct/sum(age_pct)) %>% 
  ungroup()
  # mutate(pct = n/sum(n)) %>%
  # ungroup()

tmp.data.plot$celltype <- factor(tmp.data.plot$celltype,
                                 levels = c("1",
                                   "2","3","4","5","6","7"))

p <- ggplot(tmp.data.plot,
            aes(celltype,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Cell cluster")+
  #scale_x_discrete(label = c("Hippocampus","Prefrontal cortex"))+
  scale_fill_manual(values = col.spectral(6)[4:6],
                    name="Brain Region")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(#angle = 45,
                                   hjust = 1,
                                   vjust = 1))
  
p # save as 4*5inches, landscape


##----fig 3E
###Age Propotion
tmp.data.plot <- tmp.seu[[]] %>%
  group_by(age,celltype) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

p <- ggplot(tmp.data.plot,
            aes(age,pct,fill=celltype))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Cell cluster")+
  scale_fill_manual(values = c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                               "#EE934E","#a5917f","#F5CFE4","#634d7b",
                               "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                    name="Cell type")+
  # scale_fill_manual(values = c("#ee8172","#83d0e2","#4dbdab","#7788ac","#f7b9a6",
  #                            "#a9b2cb","#b2dfd4","#a5917f","#e5b04a","#634d7b",
  #                            "#318FB5","#bdbdbd","#0773b0","#2f973d","#ec7c1d",
  #                            "#005086","#c21071"),
  # name="CCA Cluster")+
  #scale_x_discrete(label = c("Hippocampus","Prefrontal cortex"))+
  #scale_fill_manual(values = col.spectral(6)[4:6])+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 30)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))
# theme(plot.title = element_text(hjust = 0.5),
#       axis.text.x = element_text(angle = 45,
#                                  hjust = 1,
#                                  vjust = 1))
p #save as 8*8 inches,landscape


##----fig 3F
library(vegan)
age_abundance <- tmp.seu[[]] %>%
  group_by(celltype,age) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  pivot_wider(id_cols = age,
              names_from = celltype,
              values_from = n) %>% 
  column_to_rownames("age")


# Shannon-Wiener index
shannon_index <- diversity(age_abundance, index = "shannon")
print(shannon_index)
# Simpson index
simpson_index <- diversity(age_abundance, index = "simpson")

# Beta diversity
bray_curtis <- vegdist(age_abundance, method = "bray")
print(bray_curtis)
# 
tmp.data.plot <- cbind(shannon_index,simpson_index) %>% 
  as.data.frame() %>% 
  rownames_to_column("age") %>% 
  pivot_longer(cols = -age, 
               names_to = "diversity_index", 
               values_to = "value")
tmp.data.plot$diversity_index <- plyr::mapvalues(tmp.data.plot$diversity_index,
                                                 from = c("shannon_index","simpson_index"),
                                                 to = c("Shannon-Wiener Index","Inverse Simpson Index"))


tmp.data.plot$age <- factor(tmp.data.plot$age,
                            levels=c("2M","14M","24M"))

library(DEP)
p <- ggplot(tmp.data.plot,
            aes(x=age,y=value,
                color=diversity_index,
                group=diversity_index))+
  geom_point(size=3)+geom_line(size=2)+
  scale_color_manual(values = col.spectral(6)[4:6],
                     name="Diversity Index")+
  labs(title = "Diversity analysis",x="Age",y="Index Value")+
  theme_DEP2()
p # save as 4*6 inches, landscape

##----fig 3G
DefaultAssay(tmp.seu) <- "SCP"
Idents(tmp.seu) <- "celltype"
DEG <- FindAllMarkers(tmp.seu, logfc.threshold = 0.5, only.pos = T,min.pct = 0.3,
                      test.use = 'wilcox') %>% 
  mutate(pct.diff=pct.1-pct.2)

palcolor <- c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
              "#EE934E","#a5917f","#F5CFE4","#634d7b",
              "#8FA4AE","#F5D2A8","#FCED82","#BBDD78")
DEGs <- DEG[with(DEG, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
ht <- FeatureHeatmap(
  srt = tmp.seu, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$cluster,
  height = 5, width = 4,feature_split_palcolor=palcolor,group_palcolor=palcolor
)
print(ht$plot)# save as A4, landscape

#GO enrichment analysis
tmp.seu[[]]$celltype <- factor(tmp.seu[[]]$celltype,levels = c("1","2","3","4",
                                                       "5","6","7"))
ego_df.list <- lapply(unique(tmp.seu[[]]$celltype),function(xx){
  print(xx)
  sig_deg_up <- subset(DEG,avg_log2FC > 0.5&cluster== xx)
  ego <- enrichGO(gene          = row.names(sig_deg_up),
                  #universe     = row.names(dge.celltype),
                  OrgDb         = 'org.Mm.eg.db',
                  keyType       = 'SYMBOL',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
  ego_df <- data.frame(ego) %>% 
    mutate(cluster=xx)
  return(ego_df)
})

ego_all <- Reduce(rbind,ego_df.list)
