##---- 0. load packages --------
library(tidyverse)
library(ggplot2)
library(Seurat)
library(patchwork)
library(reticulate)
library(msigdbr)
library(UCell)
library(ggstatsplot)
library(RColorBrewer)
source("scp_utils.R")

tmp.path <- "res/SCP_MG_3085_seurat_clustering.rds"
tmp.seu <- readRDS(file = tmp.path)

# Figure 4A
scplotter::CellDimPlot(tmp.seu,
                       group_by = "plot_cluster", facet_by = "plot_cluster", reduction = "bbknnUMAP",
                       highlight = 'plot_cluster == "C4"', theme = "theme_blank", 
                       legend.position = "none", highlight_stroke = 0, highlight_size = 0.5,
                       palcolor = c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F"))


# Figure 4C
tmp.seu[["SCP"]] <- as(tmp.seu[["SCP"]], "Assay")
ht1 <- SCP::GroupHeatmap(tmp.seu,
                         features = c(
                           "Cx3cr1", "C1qa", "Itgam", "Itgb1", "Sirpa"), border = F,
                         group.by = c("plot_cluster"),
                         split.by = "plot_cluster",
                         cell_split_palcolor = c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F")
)

ht1$plot


# Figure 4D
msigdbr(species = "Mus musculus") %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) %>% 
  print(n = 100)

GOBP_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
str(GOBP_list)
aerobic_res <- list(GOBP_list$GOBP_AEROBIC_RESPIRATION)

gene_sets <- list(
  aerobic_res = aerobic_res
)
tmp.seu <- AddModuleScore_UCell(
  tmp.seu,
  features = gene_sets,
  name = "_UCell",
  assay = "SCP",
  slot = "counts"
)


SCP::FeatureDimPlot(tmp.seu, 
                    features = c("aerobic_res_UCell"),
                    reduction = "bbknnUMAP", pt.size = 0.5, ncol = 4,
                    palcolor = sc.hic.orange(100))

# Figure 4E
SCP::FeatureStatPlot(tmp.seu, stat.by = c("aerobic_res_UCell"), 
                     palcolor = c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F"),
                     group.by = "plot_cluster", plot_type = "violin", add_box = T,
                     comparisons = list(c("C1", "C2"),
                                        c("C4", "C3"),
                                        c("C1", "C4")
                     ))

kruskal_test_result <- kruskal.test(aerobic_res_UCell ~ plot_cluster, data = tmp.seu[[]])

# Print test results
cat("Kruskal-Wallis test results:\n")
print(kruskal_test_result)

if (kruskal_test_result$p.value < 0.05) {
  cat("\nKruskal-Wallis test is significant (p < 0.05). Performing post-hoc pairwise comparisons (Wilcoxon rank-sum test with Holm correction):\n")
  
  pairwise_results <- pairwise.wilcox.test(tmp.seu[[]]$aerobic_res_UCell, 
                                           tmp.seu[[]]$plot_cluster, 
                                           p.adjust.method = "holm") # Apply Holm correction [11](@ref)
  
  print(pairwise_results)
  
  cat("\nResult interpretation: The table shows p-values after Holm correction.\n")
  cat("Significance is typically based on an adjusted p-value < 0.05.\n")
  
} else {
  cat("\nKruskal-Wallis test result is not significant (p =", round(kruskal_test_result$p.value, 4), ").\n")
  cat("Based on standard statistical practice, post-hoc pairwise comparisons are not typically performed if the overall test is not significant.\n")
}

# Figure 4F, 4G
msigdbr(species = "Mus musculus") %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) %>% 
  print(n = 100)

GOBP_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
str(GOBP_list)
aerobic_res <- list(GOBP_list$GOBP_AEROBIC_RESPIRATION)
synapse_pruning <- list(GOBP_list$GOBP_SYNAPSE_PRUNING)

gene_sets <- list(
  aerobic_res = aerobic_res,
  synapse_pruning = synapse_pruning
)
tmp.seu <- AddModuleScore_UCell(
  tmp.seu,
  features = gene_sets,
  name = "_UCell",
  assay = "SCP",
  slot = "counts"
)

ttt <- tmp.seu[, !tmp.seu$plot_cluster == "C5"]
ttt[["SCP"]] <- as(ttt[["SCP"]], "Assay")
homeo <- FetchData(ttt, slot = "data", vars = c("P2ry12", "Tmem119", "P2ry13", "Csf1r"))
df <- cbind(homeo, ttt[[]])

df2 <- df %>% filter(P2ry12!=0)
p1 <- ggscatterstats(df2,x=P2ry12, y=aerobic_res_UCell,
                     centrality.para="mean",
                     xsidehistogram.args = list(fill = "#febf6f", 
                                                color = "black", na.rm = TRUE), 
                     ysidehistogram.args = list(fill = "#fb9a99",
                                                color = "black", na.rm = TRUE),
                     margins="both",
                     marginal = T,
                     marginal.type = "densigram",
                     ylab = "UCell Score for Aerobic respiration",
                     xlab = "P2RY12 intensity",
                     title = "Aerobic respiration activity vs P2RY12 intensity")

df2 <- df %>% filter(P2ry12!=0 & synapse_pruning_UCell!=0)
p2 <- ggscatterstats(df2,x=P2ry12, y=synapse_pruning_UCell,
                     centrality.para="mean",
                     xsidehistogram.args = list(fill = "#febf6f", 
                                                color = "black", na.rm = TRUE), 
                     ysidehistogram.args = list(fill = "#b2df8a",
                                                color = "black", na.rm = TRUE),
                     margins="both",
                     marginal = T,
                     marginal.type = "densigram",
                     ylab = "UCell Score for Synapse pruning",
                     xlab = "P2RY12 intensity",
                     title = "Synapse pruning activity vs P2RY12 intensity")


p1|p2

# Figure 4H
seurat_obj[[]] <- seurat_obj[[]] %>% 
  mutate(cluster = paste0("C", plot_cluster)) %>% 
  mutate(cluster = factor(plot_cluster, levels = c("C1", "C2", "C3", "C4", "C5")))
cells_to_keep <- rownames(seurat_obj@meta.data)[
  complete.cases(seurat_obj@meta.data[, c("sling_pseudo_1")])
]
seurat_obj_l1 <- subset(seurat_obj, cells = cells_to_keep)
pseudo <- colorRampPalette(c("#AECDE1", "#ee8172", "#EE934E")) 
SCP::FeatureDimPlot(seurat_obj_l1, features = "sling_pseudo_1", 
                    reduction = "bbknnUMAP", pt.size = 1,
                    ylab = "UMAP 2", xlab = "UMAP 1", theme_use = "theme_blank",
                    palcolor = pseudo(50)) +
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
            color = "#36aa87", size = 1.2, inherit.aes = FALSE)

# Figure 4I
DynamicPlot(
  srt = seurat_obj_l1, lineages = c("sling_pseudo_1"), group.by = "cluster",
  features = c("P2ry12", "Tmem119"),
  point_palcolor = Cluster,
  line_palcolor = "#36aa87",
  compare_lineages = TRUE, compare_features = FALSE
)

# Figure 4J
tmp.seu[["SCP"]] <- as(tmp.seu[["SCP"]], "Assay")
SCP::FeatureDimPlot(tmp.seu, c("Wdr77", "Lamp1"), pt.size = 0.5,
                    palcolor = sc.hic.orange(100),
                    reduction = "bbknnUMAP",
                    xlab = "UMAP 1",
                    ylab = "UMAP 2")

# Figure 4K
tmp.seu[["SCP"]] <- as(tmp.seu[["SCP"]], "Assay")
ht1 <- SCP::GroupHeatmap(tmp.seu, slot = "data",
                         features = c(
                           "Arg1", "Lamp1", "Wdr77", "Anxa2", "Anxa7", "Anxa11",
                           "Lgalsl", "Lgals3", "Pld2", "Gba1",
                           "Aloxe3", "Alox12b", "Tgm1", "Tgm3",
                           "Prdx2", "Prdx4", "Flna", "Flnb", "Got1"), border = F,
                         group.by = c("plot_cluster")
)
ht1$plot

# Figure 4L
p <- DotPlot(tmp.seu,
             features = c("Arg1", "Lgals3", "Apoe", "Hexb", "P2ry12", "Fscn1",
                          "Abi3", "Bin1", "Hpgds", "Basp1", "Ckb", "Rap1gds1",
                          "Crybb1", "Golm1", "Sall1", "F11r", "Itga6", "Csf1r",
                          "Pde3b", "Abcc3", "Cst3", "Cx3cr1", "Tgfbr1", "Tmem119"
             ),
             group.by = "plot_cluster",
             scale = T) +
  RotatedAxis() +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(11, "Spectral")[11:1])(100)) +
  theme(legend.position = "top") +
  scale_x_discrete(labels = function(x) toupper(x))
p

# Figure 4M
msigdbr(species = "Mus musculus") %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat) %>% 
  print(n = 100)

GOMF_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:MF") %>% 
  split(x = .$gene_symbol, f = .$gs_name)
str(GOMF_list)
peptidase <- GOMF_list$GOMF_PEPTIDASE_ACTIVITY

gene_sets <- list(
  peptidase = peptidase
)

## AddModuleScore_UCell
tmp.seu <- AddModuleScore_UCell(
  tmp.seu,
  features = gene_sets,
  name = "_UCell",
  assay = "SCP",
  slot = "counts"
)

SCP::FeatureDimPlot(tmp.seu, 
                    features = c("peptidase_UCell"),
                    reduction = "bbknnUMAP", pt.size = 0.5, ncol = 4,
                    palcolor = sc.hic.orange(100))

# Figure 4N
scplotter::FeatureStatPlot(tmp.seu, features = c("peptidase_UCell"), 
                           ident = "plot_cluster", plot_type = "violin", add_box = TRUE,
                           palcolor = c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F"))

combined_data <- FetchData(tmp.seu, vars = c("peptidase_UCell", "plot_cluster"))
shapiro_result <- shapiro.test(combined_data[, 1])
print(shapiro_result)

kruskal_test_result <- kruskal.test(peptidase_UCell ~ plot_cluster, data = combined_data)

# Print test results
cat("Kruskal-Wallis test results:\n")
print(kruskal_test_result)

if (kruskal_test_result$p.value < 0.05) {
  cat("\nKruskal-Wallis test is significant (p < 0.05). Performing post-hoc pairwise comparisons (Wilcoxon rank-sum test with Holm correction):\n")
  
  pairwise_results <- pairwise.wilcox.test(combined_data$peptidase_UCell, 
                                           combined_data$plot_cluster, 
                                           p.adjust.method = "holm") # Apply Holm correction [11](@ref)
  
  print(pairwise_results)
  
  cat("\nResult interpretation: The table shows p-values after Holm correction.\n")
  cat("Significance is typically based on an adjusted p-value < 0.05.\n")
  
} else {
  cat("\nKruskal-Wallis test result is not significant (p =", round(kruskal_test_result$p.value, 4), ").\n")
  cat("Based on standard statistical practice, post-hoc pairwise comparisons are not typically performed if the overall test is not significant.\n")
}

