#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Data Visualization                                              ##
## Data: 2025.03.31                                                                ##
## Author: Bijia Chen, Haoran Zhang                                                ##
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



scp_seu$group <- ifelse(scp_seu$celltype == 7, "BAM", "Microglia")

##----fig4A,G
scplotter::CellDimPlot(tmp.seu,
            group_by = "celltype", facet_by = "celltype", reduction = "harmonyUMAP",
            highlight = 'celltype == "3"', theme = "theme_blank", 
            legend.position = "none",highlight_stroke=0,highlight_size=0.5,
            palcolor=c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                       "#EE934E","#a5917f","#F5CFE4","#634d7b",
                       "#8FA4AE","#F5D2A8","#FCED82","#BBDD78")
)
##----fig4B
Idents(scp_seu) <- "celltype"

avg_expr <- AggregateExpression(scp_seu, group.by = "celltype")$SCP
for(i in 1:ncol(avg_expr)){colnames(avg_expr)[i] = gsub("g","C", colnames(avg_expr)[i])}

cor_matrix <- cor(as.matrix(avg_expr), method = "spearman")

p <- pheatmap(cor_matrix, 
              cluster_rows = TRUE, 
              cluster_cols = TRUE, 
              color = YlGnBu_gradient(100),
              main = "BAM vs Microglia Subgroups Similarity")
##----fig4C
Idents(scp_seu) <- "group"
diff_genes <- FindMarkers(scp_seu, ident.1 = "BAM", ident.2 = "Microglia")
diff_genes$gene <- rownames(diff_genes)

fc = 1

diff_genes$threshold <- 'Other'
diff_genes[diff_genes$p_val_adj < 0.01 & diff_genes$avg_log2FC > fc,'threshold'] <- 'Up'
diff_genes[diff_genes$p_val_adj < 0.01 & diff_genes$avg_log2FC < -fc,'threshold'] <- 'Down'
table(diff_genes$threshold)

diff_genes$threshold <- factor(diff_genes$threshold,levels = c('Up','Down','Other'))

p_volcano <- ggplot(diff_genes,aes(x=avg_log2FC,y=-log10(p_val_adj))) +
  geom_point(aes(color=threshold), size = 1, alpha=0.5) +
  scale_color_manual(values = c("Up" = "#d35230",
                                "Down" = "#2b7cd3",
                                "Other" = "gray")) +
  ggrepel::geom_text_repel(data = diff_genes[(diff_genes$p_val_adj<0.01 & abs(diff_genes$avg_log2FC) > fc) & 
                                               (-log10(diff_genes$p_val_adj) > 50 | abs(diff_genes$avg_log2FC) > 5) | 
                                               diff_genes$gene %in% c("P2ry12","Slc2a5","Tmem119","Cx3cr1",
                                                                      "Mrc1","F13a1","Stab1","Dab2") ,],
                           aes(label = gene),
                           size = 4, 
                           box.padding = 0.5, 
                           point.padding = 0.8, 
                           min.segment.length = 0.5, 
                           segment.color = "black",
                           show.legend = F) +
  geom_vline(xintercept=c(-1,1),col="gray",linetype = 'dotted', linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01),col="gray",linetype = 'dotted', linewidth = 0.5) +
  labs(x="log2 (FoldChange)",
       y="-log10 (p-adj)",
       title = "BAM vs Microglia - Differential Proteins") +
  theme_few() +
  theme(legend.position = 'none')
##----fig4D
sig_genes_up <- diff_genes %>% filter(p_val_adj < 0.01 & avg_log2FC > 1) %>% rownames()
sig_genes_down <- diff_genes %>% filter(p_val_adj < 0.01 & avg_log2FC < -1) %>% rownames()

go_enrich_up <- enrichGO(gene = sig_genes_up,
                         OrgDb = org.Mm.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",  
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05,
                         readable = TRUE)
go_enrich_up_select <- go_enrich_up %>% 
  simplify(.,cutoff=0.7,by="p.adjust",select_fun=min) %>%
  as.data.frame() %>%
  mutate(logp.adjust = -log10(p.adjust)) %>%
  slice_head(n=6) %>%
  mutate(class="BAM")

go_enrich_down <- enrichGO(gene = sig_genes_down,
                           OrgDb = org.Mm.eg.db,
                           keyType = "SYMBOL",
                           ont = "BP",  
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
go_enrich_down_select <- go_enrich_down %>% 
  as.data.frame() %>%
  subset(ID %in% c("GO:0008380","GO:0042063","GO:0010001","GO:0150063","GO:0048880")) %>% 
  mutate(logp.adjust = -log10(p.adjust)) %>%
  mutate(class="MG")

go_merge <- rbind(go_enrich_up_select,go_enrich_down_select) 
scplotter::EnrichmentPlot(go_merge, 
                          plot_type = "comparison", 
                          group_by = "class")+
  scale_fill_gradientn(colours=paletteer::paletteer_d("RColorBrewer::YlOrRd"),
                       name="-log10 (p.adj)")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))

# barplot
go_merge$GeneRatio_numeric <- sapply(strsplit(go_merge$GeneRatio, "/"),   
                                     function(x) as.numeric(x[1]) / as.numeric(x[2]))  
go_merge$Description <- factor(go_merge$Description, levels = rev(go_merge$Description))
p_go <- ggplot(go_merge,aes(x = GeneRatio_numeric,y = Description,fill=logp.adjust))  +  
  geom_bar(stat="identity",width=0.8 ) + 
  theme_few() +
  scale_fill_gradient2(low="#FFFFCCFF",mid="#FD8D3CFF",high ="#b04130",midpoint = 12.5) +  
  labs(x="GeneRatio",y=" ",title=" ",fill="-log10 (p.adj)") + 
  theme(axis.text=element_text(size=10,color="black"),
        axis.title = element_text(size=16),title = element_text(size=13)) +
  facet_wrap(~ class, nrow = 2, scales = "free_y")
##----fig4E,F
tmp.merkers <- c("P2ry12","Slc2a5","Tmem119","Cx3cr1",
                 "Mrc1","F13a1","Stab1","Dab2")
SCP::FeatureDimPlot(tmp.seu,tmp.merkers,pt.size=0.5,
                    palcolor = sc.hic.orange(100),
                    reduction="harmonyUMAP",
                    xlab="UMAP 1",
                    ylab="UMAP 2",
                    ncol = 4)
# Save as A4, landscape

##----fig4H
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
scplotter::FeatureStatPlot(ttt, features = c("Tmem119","Csf1r","P2ry12","P2yr13"), 
                           group_by = "celltype", plot_type = "violin",
                           ident = "celltype",add_point = TRUE,
                           palcolor = c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                                        "#EE934E","#a5917f","#F5CFE4","#634d7b",
                                        "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                           # comparisons = list(c("6", "4"), 
                           #                    c("6","2"),
                           #                    c("1","6")),
                           # ylab = "Expression level",
                           x_text_angle=0)
# save as 6*8, landscape
##----fig4I
ttt <- tmp.seu[,!tmp.seu$celltype == "7"]
ht1 <- SCP::GroupHeatmap(ttt,
                         features = c(
                           "Cx3cr1","C1qa","Itgam","Itgb1","Sirpa"),border=F,
                         group.by = c("celltype")
)
ht1$plot
##----fig4J
DotPlot(ttt,
        features = c("Slc17a7","Snap25","Snap23","Snap29",
                     "Syn1","Syp","Syt1","Stx4","Vamp2",
                     "Vamp4","Gap43","Stxbp1","Stx12",
                     "Stx7","Stxbp5","Gphn","Camk2a",
                     "Septin7"),
        group.by = "celltype",
        scale = T,
        cols =my_col)+
  RotatedAxis()+
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral")[11:1])(100))+
  theme(legend.position = "top")
##----fig4K
SCP::FeatureDimPlot(tmp.seu, 
                    features = c("Pre_syn_Score1","Post_syn_Score1",
                                 "Proteasome_Score1","Peptidase_Score1",
                                 "Spliceosome_Score1","MGnD_Score1",
                                 "LDAM_Score1","Splicing_Score1"),
                    reduction = "harmonyUMAP",pt.size=0.5,ncol = 4,
                    palcolor = sc.hic.orange(100))

##----fig4L
library(AUCell)
library(msigdbr) 
msigdbr(species = "Mus musculus") %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) %>% 
  print(n=100)

GOBP_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
str(GOBP_list)
aerobic_res <- list(GOBP_list$GOBP_AEROBIC_RESPIRATION)

FeatureStatPlot(ttt,stat.by = c("aerobic_Score1"), 
                palcolor = c("#ee8172","#D2EBC8","#3C77AF","#7DBFA7","#AECDE1",
                             "#EE934E","#a5917f","#F5CFE4","#634d7b",
                             "#8FA4AE","#F5D2A8","#FCED82","#BBDD78"),
                group.by = "celltype", plot_type = "box",
                comparisons = list(c("1","2"),
                                  c("4","3"),
                                  c("5", "4"), 
                                  c("5","2")))