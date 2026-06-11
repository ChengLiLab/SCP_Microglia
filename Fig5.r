##---- 0. load packages --------
library(tidyverse)
library(ggplot2)
library(Seurat)
library(UCell)
library(plotthis)
library(ggsignif)
library(pheatmap)
source("scp_utils.R")

# Figure 5A
ex.genes_2 <- read.csv("20250121_exo_365genes_aging_0.7SCPcutoff.csv")[,2] # 365 Phagoproteome proteins


ttt <- tmp.seu[, !tmp.seu$plot_cluster == "C5"] # remove BAMs
ttt[["ALL"]] <- as(ttt[["ALL"]], "Assay")
DefaultAssay(ttt) <- "ALL"

ht3 <- SCP::GroupHeatmap(ttt,
                         features = c(ex.genes_2), border = F,
                         cluster_rows = T, show_row_names = F,
                         width = 2,
                         species = "Mus_musculus", db = "GO_CC", anno_terms = TRUE,
                         anno_keys = F, anno_features = F, slot = "data",
                         features_label = c("Slc7a11", "Dmd", "Nptx1", "Cbln4",
                                            "Bsn", "Shank3",
                                            "Syp", "Gabrb2", "Gabra1", "Gria2", "Aqp4",
                                            "Syt1", "Mapt", "Gap43", "Thy1", "Dlg3", "Nefl", "Snap25", "Gfap", "Gad2", "Homer1", "Stx1a", "Camk2a"
                         ), # C5
                         split.by = "plot_cluster",
                         cell_split_palcolor = c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F"),
                         n_split = 4, group.by = c("plot_cluster"))
ht3$plot


# Figure 5B
SCP::FeatureDimPlot(tmp.seu, assay = "ALL",
                    features = c("Bsn", "Syp", "Syt1"),
                    reduction = "bbknnUMAP", pt.size = 0.5, ncol = 4,
                    xlab = "UMAP 1", ylab = "UMAP 2",
                    palcolor = sc.hic.orange(100))

# Figure 5C
# Create a reusable function
calculate_proteinset_intensity <- function(seurat_obj, geneset, 
                                           assay = "ALL", slot = "counts") {
  mat <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  existing_genes <- intersect(geneset, rownames(mat))
  if(length(existing_genes) == 0) {
    warning("Specified proteins not found")
    return(rep(0, ncol(mat)))
  }
  colSums(mat[existing_genes, , drop = FALSE])
}

protein_sets <- list(
  Phago = ex.genes_2
)

ttt <- tmp.seu[, !tmp.seu$plot_cluster == "C5"] # remove BAMs

# Batch calculation
for(set_name in names(protein_sets)) {
  ttt[[paste0(set_name, "_total")]] <- calculate_proteinset_intensity(ttt, protein_sets[[set_name]])
}

df <- ttt@meta.data
df <- df %>% 
  mutate(log_phago_total = log10(1 + Phago_total))

cluster.color <- c("#ee8172", "#7DBFA7", "#AECDE1", "#EE934E", "#a5917f")

p1 <- plotthis::ViolinPlot(df, x = "plot_cluster",
                           y = "log_phago_total", 
                           sig_label = "p.signif",
                           sig_labelsize = 4,
                           step_increase = 0.1,
                           add_point = F,
                           add_box = T,
                           pt_alpha = 0.1,
                           y_max = 0.05,
                           legend.position = "none",
                           xlab = "",
                           ylab = "log10 (Intensity + 1)",
                           box_color = "black",
                           box_width = 0.2,
                           palcolor = cluster.color) +
  ggsignif::geom_signif(
    comparisons = list(c("C3", "C4"),
                       c("C3", "C2"),
                       c("C1", "C2"),
                       c("C1", "C3")),
    test = "t.test", color = "black",
    map_signif_level = TRUE, step_increase = 0.1
  ) +
  scale_y_continuous(limits = c(0, 10),
                     n.breaks = 3,
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off")
p1


shapiro_result <- shapiro.test(df$log_phago_total)
print(shapiro_result)

kruskal_test_result <- kruskal.test(log_phago_total ~ cluster, data = df)

# Print test results
cat("Kruskal-Wallis test results:\n")
print(kruskal_test_result)

if (kruskal_test_result$p.value < 0.05) {
  cat("\nKruskal-Wallis test is significant (p < 0.05). Performing post-hoc pairwise comparisons (Wilcoxon rank-sum test with Holm correction):\n")
  
  pairwise_results <- pairwise.wilcox.test(df$log_phago_total, 
                                           df$cluster, 
                                           p.adjust.method = "holm") # Apply Holm correction [11](@ref)
  
  print(pairwise_results)
  
  cat("\nResult interpretation: The table shows p-values after Holm correction.\n")
  cat("Significance is typically based on an adjusted p-value < 0.05.\n")
  
} else {
  cat("\nKruskal-Wallis test result is not significant (p =", round(kruskal_test_result$p.value, 4), ").\n")
  cat("Based on standard statistical practice, post-hoc pairwise comparisons are not typically performed if the overall test is not significant.\n")
}



# Figure 5E

df <- ht3$enrichment$input
P1 <- df %>% filter(geneID_groups == "C1") %>% 
  dplyr::select("geneID") %>% unlist()
P2 <- df %>% filter(geneID_groups == "C2") %>% 
  dplyr::select("geneID") %>% unlist()
P3 <- df %>% filter(geneID_groups == "C3") %>% 
  dplyr::select("geneID") %>% unlist()
P4 <- df %>% filter(geneID_groups == "C4") %>% 
  dplyr::select("geneID") %>% unlist()


gene_sets <- list(
  P1 = P1,
  P2 = P2,
  P3 = P3,
  P4 = P4)


tmp.seu <- AddModuleScore_UCell(
  tmp.seu,
  features = gene_sets,
  name = "_UCell",
  assay = "ALL",
  slot = "counts")

SCP::FeatureDimPlot(tmp.seu, 
                    features = c("P1_UCell", "P2_UCell", "P3_UCell", "P4_UCell"),
                    reduction = "bbknnUMAP", pt.size = 0.5,
                    xlab = "UMAP 1", ylab = "UMAP 2", ncol = 4,
                    palcolor = sc.hic.orange(100))

# Figure 5F
exo_module <- c("P1_UCell", "P2_UCell", "P3_UCell", "P4_UCell")
data <- GetAssayData(tmp.seu, slot = "data", assay = "SCP") %>%
  t() %>% as.data.frame()

receptors <- read.table("public data/mouse_lr_pair.txt", header = T) %>% 
  dplyr::select("receptor_gene_symbol") %>% as.vector() %>% unlist() %>% unique()
expressed_receptors <- receptors[receptors %in% colnames(data)]

df <- cbind(data[, c(colnames(data) %in% expressed_receptors)],
            tmp.seu@meta.data[, c(colnames(tmp.seu[[]]) %in% exo_module)]) %>% 
  t()


cor_results <- list()
for(Phago in expressed_receptors) {
  for(exo in exo_module) {
    # Get expression vectors
    phago_expr <- df[Phago, ]
    exo_expr <- df[exo, ]
    
    # Remove zero values, find non-zero cells where both genes are expressed
    non_zero_idx <- which(phago_expr > 0 & exo_expr > 0)
    
    # Check if there are enough data points
    if(length(non_zero_idx) >= 5) {
      # Perform correlation analysis using non-zero values
      phago_nonzero <- phago_expr[non_zero_idx]
      exo_nonzero <- exo_expr[non_zero_idx]
      
      cor_val <- cor(phago_nonzero, exo_nonzero, method = "spearman")
      p_val <- cor.test(phago_nonzero, exo_nonzero, method = "spearman")$p.value
      
      # Save results
      cor_results[[paste(Phago, exo, sep = "_")]] <- list(
        Phago = Phago,
        exo = exo,
        correlation = cor_val,
        p_value = p_val,
        n_nonzero = length(non_zero_idx)  # Add number of non-zero cells
      )
    } else {
      # Insufficient data points, save as NA
      cor_results[[paste(Phago, exo, sep = "_")]] <- list(
        Phago = Phago,
        exo = exo,
        correlation = NA,
        p_value = NA,
        n_nonzero = length(non_zero_idx)  # Record actual number of non-zero cells
      )
    }
  }
}

cor_df <- do.call(rbind, lapply(cor_results, as.data.frame))
cor_df$fdr <- p.adjust(cor_df$p_value, method = "fdr")
cor_df <- cor_df %>% 
  mutate(significant = case_when(
    fdr < 0.05 & abs(correlation) > 0.1 ~ "Sig",
    .default = "Not Sig")
  )
cor_df$correlation <- ifelse(cor_df$significant == "Not Sig", NA, cor_df$correlation)
significant_receptor <- cor_df %>% filter(significant == "Sig")

cor_mt <- cor_df[, c(1:3)] %>% pivot_wider(
  names_from = exo,      # Column names from the exo column
  values_from = correlation,  # Values from the correlation column
  values_fill = NA       # Fill missing values with NA
) %>%
  as.data.frame() %>% 
  column_to_rownames("Phago") %>% 
  as.matrix()

cor_mt <- cor_mt[rowSums(is.na(cor_mt)) < ncol(cor_mt), ]
gene_ranking <- data.frame(
  gene = rownames(cor_mt),
  P1_cor = cor_mt[, "P1_UCell"],
  P2_cor = cor_mt[, "P2_UCell"],
  P3_cor = cor_mt[, "P3_UCell"],
  P4_cor = cor_mt[, "P4_UCell"],
  stringsAsFactors = FALSE
)

# Determine which sample has the highest correlation for each gene
gene_ranking$max_sample <- apply(gene_ranking[, 2:5], 1, function(x) {
  colnames(gene_ranking[, 2:5])[which.max(x)]
})

# Sort as requested: P1 highest > P2 highest > P3 highest > P4 highest
gene_ranking <- gene_ranking %>%
  arrange(
    # First sort by the sample with the maximum correlation
    factor(max_sample, levels = c("P4_cor", "P3_cor", "P2_cor", "P1_cor")),
    # Then sort in descending order of correlation within each group
    case_when(
      max_sample == "P1_cor" ~ P1_cor,
      max_sample == "P2_cor" ~ P2_cor,
      max_sample == "P3_cor" ~ P3_cor,
      max_sample == "P4_cor" ~ P4_cor
    )
  )

# Get the sorted gene order
sorted_genes <- gene_ranking$gene

# Reorder the matrix
sorted_matrix <- cor_mt[sorted_genes, ]
rownames(sorted_matrix) <- toupper(rownames(sorted_matrix))
sorted_matrix[is.na(sorted_matrix)] <- 0

cor_color <- colorRampPalette(c(
  "#2E86AB", "#6AB5D0", "#A5D9F0", "#E0F3F8", 
  "#FFFFFF", "#FFE8D6", "#FFB8A1", "#FF7F6F", 
  "#D43D51"
))(100)

breaks <- seq(-0.7, 0.7, length.out = 101)

## All receptor
pheatmap(sorted_matrix, cluster_cols = F,
         cluster_rows = F,
         color = cor_color,
         breaks = breaks,
         main = "Phagoproteome vs Receptors (Spearman)")


# Figure 5G
SCP::FeatureDimPlot(tmp.seu, 
                    features = c("Itgb3", "Cd44", "Itgam", "Atp1a3"),
                    reduction = "bbknnUMAP", pt.size = 0.5, ncol = 4,
                    xlab = "UMAP 1", ylab = "UMAP 2",
                    palcolor = sc.hic.orange(100))
