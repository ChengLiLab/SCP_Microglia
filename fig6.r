#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Main Figure                                                     ##
## Data: 2025.03.21                                                                ##
## Author: Guangxin Zhang, Haoran Zhang                                            ##
#####################################################################################


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
library(ggside)
library(tidyverse)
library(DEP)
setwd('/Users/guangxinzhang/Documents/new_analyze')
source("code/scp_utils.R")

seurat_obj <- readRDS("seurat_object/SCP_MG_3184_seurat_zhr_step1_harmony_20250217_label_transferred_morpho_score_v5.rds")

####################### Fig 6B #######################

df <- data.frame(
  Diameter = seurat_obj@meta.data$Diameter,
  nFeature_SCP = seurat_obj@meta.data$nFeature_SCP
)

cor_test <- cor.test(df$Diameter, df$nFeature_SCP)
cor_value <- cor_test$estimate  
p_value <- cor_test$p.value     

age_info <- as.factor(seurat_obj$age)
color_age <- col.spectral(6)[4:6]

p <- ggplot(df, mapping = aes(Diameter, nFeature_SCP)) +
  geom_point(size = 0.9, color = "#7BA4CC", alpha = 0.4) +
  geom_point(mapping = aes(color = age_info), size = 1) +
  scale_color_manual(values = color_age) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 4))) +
  labs(x = 'Diameter (μm) ', y = 'Number of Protein Groups') +
  theme_classic(base_size = 20) + 
  theme(
    legend.position = c(0.65, 0.7),  
    legend.background = element_blank(),
    legend.text = element_text(size = 20),
    legend.box.margin = margin(t = -20),
    legend.key.spacing = unit(0, 'cm'),
    legend.title = element_blank()
  )

p_6b <- p + 
  geom_xsidehistogram(
    aes(y = after_stat(density)),
    bins = 20, 
    fill = "#7BA4CA"
  ) +
  geom_xsidedensity(
    aes(y = stat(density)),
    color = "#CD6453", size = 2
  ) +
  scale_xsidey_continuous(labels = NULL) +
  geom_ysidehistogram(
    aes(x = after_stat(density)),
    bins = 20, 
    fill = "#7BA4CA"
  ) +
  geom_ysidedensity(
    aes(x = stat(density)),
    color = "#CD6453", size = 2
  ) +
  scale_ysidex_continuous(labels = NULL) +  
  theme(
    ggside.panel.scale = 0.2,
    ggside.axis.line = element_blank(),
    ggside.axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = element_line(size = 2),
    axis.ticks = element_line(size = 2)
  ) +
  annotate(
    "text", 
    x = min(df$Diameter) + 0.35, 
    y = max(df$nFeature_SCP) , 
    label = paste("R =", round(cor_value, 2), "\nP =", signif(p_value, 3)),
    hjust = 0, 
    vjust = 1, 
    size = 7, 
    color = "black"
  )

print(p_6b)




####################### Fig 6C #######################

p_6c <- ggplot(df, aes(x = correlation, y = -log10(p_value))) +
  geom_point(aes(color = significant)) +
  scale_color_manual(values=c("#00aaff","#d2dae2","#ff8802"))+
  labs(title = "Correlation with Diameter (Pearson)", 
       x = "Correlation Coefficient", 
       y = "-log10(p value)") +
  ggrepel::geom_text_repel(
    data = filter(df, significant != "Not sig"), 
    aes(label = toupper(gene)),  
    size = 3, 
    box.padding = unit(0.1, "lines"), 
    point.padding = unit(0.1, "lines"), 
    segment.size = 0.5
  ) +
  ggrepel::geom_text_repel(
    data = filter(df, gene %in% c("Plek","Anxa1","Anxa2","P2ry12","Cltb")), 
    aes(label = toupper(gene)),  
    size = 3, 
    color = "red",
    box.padding = unit(0.1, "lines"), 
    point.padding = unit(0.1, "lines"), 
    segment.size = 0.5
  ) +
  geom_point(
    data = filter(df, gene %in% c("Plek","Anxa1","Anxa2","P2ry12","Cltb")),
    aes(label = toupper(gene)),  
    color = "red"
  ) +
  theme_DEP2() +
  theme(legend.position = "none")

p_6c

####################### Fig 6E #######################

df <- cor_res_dia %>% subset(significant %in% c("Neg","Pos")) %>% 
  group_by(significant) %>% arrange(p_value) %>% 
   ungroup()

ego_df.list <- lapply(unique(df$significant),function(xx){
  print(xx)
  pos_cor <- subset(df,significant== xx)
  ego <- enrichGO(gene          = pos_cor$gene,
                  OrgDb         = 'org.Mm.eg.db',
                  keyType       = 'SYMBOL',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
  ego_df <- data.frame(ego) %>% 
    mutate(significant=xx)
  return(ego_df)
})

ego_all_dia <- Reduce(rbind,ego_df.list)

ego_top <- ego_all_dia %>%
  group_by(significant) %>%
  arrange(p.adjust) %>% 
  dplyr::filter(ONTOLOGY == "BP") %>%
  slice_head(n = 15) %>%
  ungroup()

neglog10_trans <- scales::trans_new(
  "neglog10",
  transform = function(x) -log10(x),
  inverse   = function(x) 10^(-x)
)

ego_BP_show <- ego_all_dia %>% 
  subset(ID %in% c("GO:0050764","GO:1902774","GO:1902774","GO:0045333","GO:0008380","GO:0031507",
                 "GO:0016072","GO:0005681","GO:0002699"))

ego_BP_show$negLogP <- -log10(ego_BP_show$p.adjust + 1e-300)

p_6e <- scplotter::EnrichmentPlot(
  ego_BP_show,
  plot_type = "comparison", 
  group_by = "significant",
  color_by = "negLogP"            
) +
  labs(title = "Diameter \n GO:BP enrichment analysis") +
  scale_fill_gradientn(
    colours = paletteer::paletteer_d("RColorBrewer::YlOrRd"),
    name = "-log10(p.adj)"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p_6e



####################### Fig 6F #######################

diameter <- scplotter::FeatureStatPlot(seurat_obj, features = c("Diameter"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
   labs(y = "Diameter (\u03bc m)", x = "Cell Clusters")

ident_pro <- scplotter::FeatureStatPlot(seurat_obj, features = c("nFeature_SCP"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
 labs(y = "# of Identified Protein Groups",x = "Cell Clusters")

 intensity <- scplotter::FeatureStatPlot(seurat_obj, features = c("nCount_SCP"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
  labs(y = "Total Intensity",x = "Cell Clusters")

 clta <- scplotter::FeatureStatPlot(seurat_obj, features = c("Clta"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
   labs(y = "Expression Level",x = "Cell Clusters")

 cltc <- scplotter::FeatureStatPlot(seurat_obj, features = c("Cltc"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
   labs(y = "Expression Level", x = "Cell Clusters")

 hp1bp3 <- scplotter::FeatureStatPlot(seurat_obj, features = c("Hp1bp3"), ident = "celltype",
   add_bg = TRUE, add_box = TRUE,
   comparisons = list(c("2", "5"),c("4", "5"),c("5", "7")),
   palcolor = groups_color) +
   theme(legend.position = "none") +
   labs(y = "Expression Level", x = "Cell Clusters")

p_6f <- cowplot::plot_grid(diameter, ident_pro, intensity, clta, cltc, hp1bp3, ncol = 3)

p_6f
