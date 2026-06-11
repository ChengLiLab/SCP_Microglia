# Fig3.R
# Purpose: Dimensionality Reduction, Clustering, and Visualization of Expressed Protein Matrix

# 0.Load Packages -------------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(bbknnR)
library(patchwork)
library(ggthemes)
library(cluster)
library(vegan)
library(igraph) 
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(SCP)

#-----------------1.Load Data----------------------------------------------------
tmp.seu <- readRDS("res/SCP_MG_3085_seurat_post_QC.rds")

ex_df <- read.csv("data/exo_genes_list.csv", header = F)
ex_genes <- ex_df$V1

seu <- tmp.seu[setdiff(rownames(tmp.seu), ex_genes), ]


#----------2.Remove batch effect and clustering ------------------------------
tmp.run.seurat <- function(seu,tmp.dims=20){
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu,npcs = tmp.dims)
  seu <- RunUMAP(seu,
                 dims = 1:tmp.dims,
                 reduction = "pca")
  return(seu)
}
tmp.process.harmony <- function(seu){
  seu[["SCP"]] <- split(seu[["SCP"]], f = seu$column_id)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu,npcs = 20)
  seu <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE)
  seu[["SCP"]] <- JoinLayers(seu[["SCP"]])
  seu <- RunUMAP(seu,dims = 1:20,
                 reduction = "harmony",
                 reduction.name = "harmonyUMAP")
  return(seu)
}
tmp.process.CCA <- function(seu){
  seu[["SCP"]] <- split(seu[["SCP"]], f = seu$column_id)
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu,npcs = 20)
  seu <- IntegrateLayers(
    object = seu,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "cca",
    verbose = FALSE)
  seu[["SCP"]] <- JoinLayers(seu[["SCP"]])
  seu <- RunUMAP(seu,dims = 1:20,
                 reduction = "cca",
                 reduction.name = "CCAUMAP")
  return(seu)
}
seu <- tmp.run.seurat(seu,tmp.dims = 20)
seu <- tmp.process.harmony(seu = seu)
seu <- tmp.process.CCA(seu = seu)
seu <- RunBBKNN(object = seu,
                batch_key = "column_id",
                reduction = "pca",
                n_pcs = 20,
                run_TSNE = F,
                run_UMAP = T,
                UMAP_key = "bbknnUMAP_",
                UMAP_name = "bbknnUMAP")
Reductions(seu)


##### clustering
run_clustering <- function(seu, 
                           reductions = c("pca", "harmony", "cca"),
                           dims = 1:20,
                           resolution = 0.3) {
  
  for (reduc in reductions) {
    if (!reduc %in% Reductions(seu)) {
      warning(paste("Reduction", reduc, "not found, skipping..."))
      next
    }
    
    cluster_name <- paste0(reduc, "_cluster")
    
    seu <- FindNeighbors(seu, 
                         reduction = reduc, 
                         dims = dims) %>%
      FindClusters(resolution = resolution, 
                   cluster.name = cluster_name)
  }
  
  return(seu)
}
seu <- run_clustering(seu,
                      reductions = c("pca", "harmony", "cca"),
                      dims = 1:20,
                      resolution = 0.3)
seu <- FindClusters(seu,
                    resolution = 0.3,
                    graph.name = "SCP_bbknn",
                    cluster.name = "bbknn_cluster")

plot_clusters <- function(seu) {
  cluster_cols <- grep("_cluster$", colnames(seu@meta.data), value = TRUE)
  plots <- lapply(cluster_cols, function(col) {
    seu@meta.data[[col]] <- as.numeric(seu@meta.data[[col]])
    
    reduc_method <- strsplit(col, "_")[[1]][1]
    
    reduc_umap <- switch(reduc_method,
                         "pca" = "umap",
                         "harmony" = "harmonyUMAP",
                         "cca" = "CCAUMAP",
                         "bbknn" = "bbknnUMAP",
                         "umap")  
    
    if (!reduc_umap %in% Reductions(seu)) {
      warning(paste("Reduction", reduc_umap, "not found for", col))
      return(NULL)
    }
    
    set.seed(42)
    DimPlot(seu, 
            reduction = reduc_umap, 
            group.by = col,
            label = T,
            label.size = 7) + 
      scale_color_manual(values = randomcoloR::distinctColorPalette(20)) + 
      labs(title=paste(reduc_umap,"\n(0.3 resolution)"),
           x=paste(reduc_umap,"1",
                   sep=" "),
           y=paste(reduc_umap,"2",
                   sep=" ")) +
      theme(legend.position = "none")
  })
  wrap_plots(plots, ncol = 2)
}
plot_clusters(seu)


## *Optimal resolution ------------------------------------------------------
resolutions <- seq(0.1, 1.5, by = 0.1)
seu1 <- FindClusters(seu,
                     graph.name = "SCP_bbknn",
                     resolution = resolutions,
                     verbose = FALSE)

calculate_clustering_metrics <- function(seu_obj, resolution, reduction = "bbknnUMAP") {
  
  cluster_col <- paste0("SCP_bbknn_res.", resolution)
  clusters <- seu_obj@meta.data[[cluster_col]]
  
  reduction_data <- Embeddings(seu_obj, reduction)[,1:2]
  
  clusters_numeric <- as.numeric(as.factor(clusters))
  
  dist_matrix <- dist(reduction_data)
  
  # 1. Silhouette Index
  sil_result <- cluster::silhouette(clusters_numeric, dist_matrix)
  silhouette_score <- mean(sil_result[, "sil_width"])
  
  # 2. Calinski-Harabasz Index
  variance_ratio <- fpc::calinhara(reduction_data, clusters_numeric)
  
  # 3. Davies-Bouldin Index
  db_index <- clusterSim::index.DB(reduction_data, clusters_numeric)$DB
  
  return(list(
    resolution = resolution,
    silhouette = silhouette_score,
    variance_ratio = variance_ratio,
    db_index = db_index
  ))
}

metrics_list <- lapply(resolutions, function(r) {
  tryCatch({
    calculate_clustering_metrics(seu1, r)
  }, error = function(e) {
    warning("Resolution ", r, " calculation failed: ", e$message)
    return(list(resolution = r, silhouette = NA, variance_ratio = NA, db_index = NA))
  })
})

metrics_df <- do.call(rbind, lapply(metrics_list, as.data.frame))

# normalization
normalize_metrics <- function(df) {
  df_normalized <- df
  
  df_normalized$silhouette_norm <- (df$silhouette - min(df$silhouette, na.rm = TRUE)) / 
    (max(df$silhouette, na.rm = TRUE) - min(df$silhouette, na.rm = TRUE))
  
  df_normalized$variance_ratio_norm <- (df$variance_ratio - min(df$variance_ratio, na.rm = TRUE)) / 
    (max(df$variance_ratio, na.rm = TRUE) - min(df$variance_ratio, na.rm = TRUE))
  
  df_normalized$db_index_norm <- 1 - ((df$db_index - min(df$db_index, na.rm = TRUE)) / 
                                        (max(df$db_index, na.rm = TRUE) - min(df$db_index, na.rm = TRUE)))
  
  return(df_normalized)
}

# 4. Composite Score
calculate_composite_score <- function(normalized_df, weights = c(0.4, 0.4, 0.2)) {
  normalized_df$composite_score <- 
    weights[1] * normalized_df$silhouette_norm +
    weights[2] * normalized_df$variance_ratio_norm + 
    weights[3] * normalized_df$db_index_norm
  
  return(normalized_df)
}

metrics_normalized <- metrics_df %>%
  normalize_metrics() %>%
  calculate_composite_score()

optimal_res_composite <- metrics_normalized$resolution[which.max(metrics_normalized$composite_score)]


plot_data_long <- metrics_normalized %>%
  dplyr::select(resolution, silhouette_norm, variance_ratio_norm, db_index_norm, composite_score) %>%
  pivot_longer(cols = -resolution, names_to = "metric", values_to = "normalized_score") %>%
  mutate(metric = case_when(
    metric == "silhouette_norm" ~ "Silhouette (Norm)",
    metric == "variance_ratio_norm" ~ "Calinski-Harabasz (Norm)",
    metric == "db_index_norm" ~ "Davies-Bouldin (Norm)",
    metric == "composite_score" ~ "Composite Score"
  )) %>%
  mutate(metric = factor(metric, levels = c("Silhouette (Norm)",
                                            "Calinski-Harabasz (Norm)",
                                            "Davies-Bouldin (Norm)",
                                            "Composite Score")))

# Create custom color scale
custom_colors <- c(
  "Silhouette (Norm)" = "#488cbc",       
  "Calinski-Harabasz (Norm)" = "#fc8f45", 
  "Davies-Bouldin (Norm)" = "#44a64c",    
  "Composite Score" = "black"      
)

# Create line type scale
custom_linetypes <- c(
  "Silhouette (Norm)" = "solid",
  "Calinski-Harabasz (Norm)" = "solid",
  "Davies-Bouldin (Norm)" = "solid",
  "Composite Score" = "31"
)

# Create improved plot
p <- ggplot(plot_data_long, aes(x = resolution, y = normalized_score, 
                                color = metric, linetype = metric)) +
  geom_line(size = 1, alpha = 0.9) +
  geom_point(size = 1.5, alpha = 0.9) +
  # Add red dotted vertical line for best composite score
  geom_vline(xintercept = optimal_res_composite, 
             linetype = "dotted", color = "red", size = 1) +
  # Add label for best resolution
  geom_text(data = data.frame(x = 0.5, y = 0.1, 
                              label = paste("Best Res =", optimal_res_composite)),
            aes(x = x, y = y, label = label),
            color = "red", size = 4, fontface = "bold", inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(0.1, 1.5, by = 0.2)) +
  scale_y_continuous(limits = c(0, 1.1)) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = custom_linetypes) +
  labs(title = "Clustering Performance vs. Resolution",
       x = "Resolution",
       y = "Normalized Score",
       color = "",
       linetype = "") +
  theme_few(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = c(0.75,0.8),
    legend.key.width = unit(1, "cm")
  )

## *Cluster stability analysis ------------------------------------------------------
## Bootstrap
# parameter setting
n_bootstraps <- 50  
subsample_ratio <- 0.8  

consensus_matrix <- matrix(0,
                           nrow = ncol(seu),
                           ncol = ncol(seu),
                           dimnames = list(colnames(seu), colnames(seu)))

pb <- txtProgressBar(min = 0, max = n_bootstraps, style = 3)
for (i in 1:n_bootstraps) {
  setTxtProgressBar(pb, i)
  
  set.seed(i)
  subsample_cells <- sample(colnames(seu), size = round(subsample_ratio * ncol(seu)), replace = F)
  seu_sub <- seu[, subsample_cells]
  
  seu_sub <- RunBBKNN(seu_sub,
                      batch_key = "column_id",
                      reduction = "pca",
                      n_pcs = 20)
  
  seu_sub <- FindClusters(seu_sub,
                          resolution = 0.3,
                          graph.name = "SCP_bbknn",
                          cluster.name = "boot_cluster")
  
  cluster_labels <- seu_sub$boot_cluster
  cell_indices <- match(names(cluster_labels), colnames(consensus_matrix))
  consensus_matrix[cell_indices, cell_indices] <- consensus_matrix[cell_indices, cell_indices] +
    outer(cluster_labels, cluster_labels, "==")
}
close(pb)

consensus_matrix <- consensus_matrix / n_bootstraps

cell_order <- order(seu$bbknn_cluster)
p <- ComplexHeatmap::Heatmap(consensus_matrix[cell_order, cell_order],
                             name = "Co-clustering\nFrequency",
                             col = colorRampPalette(c("white", "blue"))(100),
                             cluster_rows = FALSE,
                             cluster_columns = FALSE,
                             show_row_names = FALSE,
                             show_column_names = FALSE,
                             row_split = seu$bbknn_cluster[cell_order],
                             column_split = seu$bbknn_cluster[cell_order],
                             row_title = paste0("C", unique(as.numeric(seu$bbknn_cluster[cell_order]))),
                             column_title = paste0("C", unique(as.numeric(seu$bbknn_cluster[cell_order]))),
                             border = TRUE)


#---------- 3.Visualization ------------------------------------------------------
##-------------3.1 UMAP plot ----------------------------------------------
seu@meta.data$plot_cluster <- plyr::mapvalues(x = as.character(as.numeric(seu@meta.data$bbknn_cluster)),
                                              from = c("1","2","3","4","5"),
                                              to = c("C1","C2","C3","C4","C5")) %>%
  factor(levels = c("C1","C2","C3","C4","C5"))
table(seu@meta.data$plot_cluster)
cluster_labels <- c("C1"="C1 (1,132)",
                    "C2"="C2 (800)",
                    "C3"="C3 (631)",
                    "C4"="C4 (392)",
                    "C5"="C5-BAM (130)")
DimPlot(seu,
        reduction = "bbknnUMAP",
        group.by = "plot_cluster",
        pt.size = 0.5,
        label = T,
        label.size = 5) +
  scale_color_manual(values = c("#ee8172","#83c1a7","#AECDE1","#EE934E","#a5917f"),
                     labels = cluster_labels, name = "") +
  labs(title = "",
       x = "bbknnUMAP 1",
       y = "bbknnUMAP 2") +
  theme(legend.text = element_text(face = "bold", color = "black"))

### Add "ALL" assay
dim(GetAssayData(tmp.seu, layer ="counts"))
scp_count_all <- GetAssayData(tmp.seu, layer = "counts")
seu[["ALL"]] <- CreateAssayObject(counts = scp_count_all)
seu <- NormalizeData(seu, assay = "ALL")
seu <- ScaleData(seu, assay = "ALL")

saveRDS(seu, file = "res/SCP_MG_3085_seurat_clustering.rds")

##-----------------3.2 stack plot ----------------------------------------------
## Age Proportion
tmp.cluster <- "plot_cluster"
tmp.data.plot <- seu[[]] %>%
  mutate(celltype = as.factor(as.numeric(!!sym(tmp.cluster)))) %>%
  group_by(celltype,age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

tmp.data.all <- seu[[]] %>%
  group_by(age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  mutate(celltype = "all") # add "all"

tmp.data.plot <- bind_rows(tmp.data.plot, tmp.data.all) %>%
  mutate(age = factor(age, levels = c("2M","14M","24M"))) %>%
  mutate(celltype_num = as.numeric(gsub("all", "", celltype)),
         celltype = factor(celltype,levels = c(sort(unique(na.omit(celltype_num))),"all"))) %>%
  mutate(celltype = ifelse(is.na(celltype_num), 
                           as.character(celltype), 
                           paste0("C", celltype))) %>%
  mutate(celltype = fct_relevel(celltype, "all", after = Inf))

p <- ggplot(tmp.data.plot,
            aes(celltype,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Cell cluster")+
  scale_fill_manual(values = c("#CEEB9C","#5BB6A9","#5E4FA2"),
                    name = "Age")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 20)+
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5))


## Cluster Proportion
tmp.data.plot <- seu[[]] %>%
  mutate(celltype = as.factor(as.numeric(!!sym(tmp.cluster)))) %>%
  group_by(age,celltype) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

tmp.data.plot <- tmp.data.plot %>%
  mutate(age = factor(age, levels = c("2M","14M","24M"))) %>%
  mutate(celltype = factor(celltype,levels = c(sort(unique(celltype))))) %>%
  mutate(celltype = paste0("C", celltype))

p <- ggplot(tmp.data.plot,
            aes(age,pct,fill=celltype))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Age")+
  scale_fill_manual(values = c("#ee8172","#83c1a7","#AECDE1","#EE934E","#a5917f"),
                    name="Cluster")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 20)+
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5))


## Brain Region Proportion
tmp.data.plot <- seu[[]] %>%
  mutate(celltype = as.factor(as.numeric(!!sym(tmp.cluster)))) %>%
  group_by(celltype,brain_region) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

tmp.data.all <- seu[[]] %>%
  group_by(brain_region) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  mutate(celltype = "all") # add "all"

tmp.data.plot <- bind_rows(tmp.data.plot, tmp.data.all) %>%
  mutate(celltype_num = as.numeric(gsub("all", "", celltype)),
         celltype = factor(celltype,levels = c(sort(unique(na.omit(celltype_num))),"all"))) %>%
  mutate(celltype = ifelse(is.na(celltype_num), 
                           as.character(celltype), 
                           paste0("C", celltype))) %>%
  mutate(celltype = fct_relevel(celltype, "all", after = Inf))

p <- ggplot(tmp.data.plot,
            aes(celltype,pct,fill=brain_region))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("")+
  ylab("Proportion (%)")+
  xlab("Cell cluster")+
  scale_fill_manual(values = c("#00758b","#9cca62"),
                    name="Brain Region")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 20)+
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5))

##-----------------3.3 find markers----------------------------------------------
Idents(seu) <- seu@meta.data$plot_cluster

all_markers <- FindAllMarkers(seu, 
                              logfc.threshold = 0.5, 
                              only.pos = T,
                              min.pct = 0.3,
                              test.use = 'wilcox') %>%
  mutate(pct.diff = pct.1-pct.2)

DEGs <- all_markers[with(all_markers, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
table(DEGs$cluster)


# GO enrichment analysis
seu[[]]$plot_cluster <- factor(seu[[]]$plot_cluster,levels = c("C1","C2","C3","C4","C5"))
ego_df.list <- lapply(unique(seu[[]]$plot_cluster),function(xx){
  print(xx)
  sig_deg_up <- subset(DEGs, avg_log2FC > 0.5 & cluster== xx)
  ego <- enrichGO(gene          = row.names(sig_deg_up),
                  OrgDb         = 'org.Mm.eg.db',
                  keyType       = 'SYMBOL',
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
  ego_simplify <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)  
  ego_df <- data.frame(ego_simplify) %>% 
    mutate(cluster=xx)
  return(ego_df)
})
go_results <- Reduce(rbind,ego_df.list)


#------------------------4.Diversity analysis----------------------------------------------
age_abundance <- seu[[]] %>%
  group_by(plot_cluster,age) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  pivot_wider(id_cols = age,
              names_from = plot_cluster,
              values_from = n) %>% 
  column_to_rownames("age")

# Shannon-Wiener index
shannon_index <- diversity(age_abundance, index = "shannon")
# Simpson index
simpson_index <- diversity(age_abundance, index = "simpson")

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

p <- ggplot(tmp.data.plot,
            aes(x=age,y=value,
                color=diversity_index,
                group=diversity_index))+
  geom_point(size=3)+
  geom_line(size=2)+
  scale_color_manual(values = c("#CEEB9C","#5BB6A9"),,
                     name = "Diversity Index")+
  scale_y_continuous(limits = c(0.65,1.55))+
  labs(title = "Diversity analysis",x="Age",y="Index Value")+
  theme_bw()


#--------- 5.Age-dependent differential protein analysis by cluster ----------------------------------------------
compute_pairwise_dep_by_cluster <- function(seu_obj,
                                            comparisons = list(
                                              c("24M", "2M"),
                                              c("24M", "14M"),
                                              c("14M", "2M")
                                            ),
                                            cluster_col = "plot_cluster",
                                            age_col = "age",
                                            min_cells = 10,
                                            logfc_threshold = 0.5,
                                            pval_threshold = 0.05) {
  
  clusters <- unique(seu_obj@meta.data[[cluster_col]])
  all_results <- list()
  
  for (comparison_idx in seq_along(comparisons)) {
    comp <- comparisons[[comparison_idx]]
    ident1 <- comp[1]
    ident2 <- comp[2]
    comp_name <- paste0(ident1, " vs. ", ident2)
    
    cat("=== Comparison:", comp_name, "===\n")
    
    for (cluster in clusters) {
      cat("  Cluster:", cluster, "\n")
      
      cluster_cells <- colnames(seu_obj)[seu_obj@meta.data[[cluster_col]] == cluster]
      
      if (length(cluster_cells) < min_cells) {
        cat("    Too few cells; skipped.\n")
        next
      }
      
      seu_subset <- subset(seu_obj, cells = cluster_cells)
      Idents(seu_subset) <- seu_subset@meta.data[[age_col]]
      
      ages_present <- unique(Idents(seu_subset))
      if (!(ident1 %in% ages_present && ident2 %in% ages_present)) {
        cat("    Missing", ident1, "or", ident2, "cells; skipped.\n")
        next
      }
      
      markers <- FindMarkers(
        seu_subset,
        ident.1 = ident1,
        ident.2 = ident2,
        min.pct = 0.1,
        logfc.threshold = logfc_threshold,
        test.use = "wilcox"
      )
      
      if (nrow(markers) > 0) {
        markers$cluster <- cluster
        markers$protein <- rownames(markers)
        markers$comparison <- comp_name
        markers$ident1 <- ident1
        markers$ident2 <- ident2
        
        all_results[[paste0(comp_name, "_", cluster)]] <- markers
        
        cat("    Differential proteins identified:", nrow(markers), "\n")
      } else {
        cat("    No differential proteins found.\n")
      }
    }
  }
  
  if (length(all_results) > 0) {
    combined_results <- bind_rows(all_results)
    return(combined_results)
  } else {
    warning("No differential proteins were found.")
    return(NULL)
  }
}

prepare_comparison_data <- function(dep_results,
                                    comparison_name,
                                    fc_threshold = 0.5,
                                    pval_threshold = 0.05) {
  
  comp_data <- dep_results %>%
    filter(comparison == comparison_name) %>%
    filter(abs(avg_log2FC) > fc_threshold, p_val_adj < pval_threshold) %>%
    mutate(
      type = ifelse(avg_log2FC > 0, "Up", "Down"),
      comparison_label = comparison_name
    )
  
  if (nrow(comp_data) == 0) {
    warning(paste("No significant differential proteins found for", comparison_name))
    return(NULL)
  }
  
  comp_data <- comp_data %>%
    group_by(cluster, type) %>%
    mutate(
      is_up_top = ifelse(
        type == "Up",
        avg_log2FC %in% head(sort(avg_log2FC[type == "Up"], decreasing = TRUE), 5),
        FALSE
      ),
      is_down_top = ifelse(
        type == "Down",
        avg_log2FC %in% head(sort(avg_log2FC[type == "Down"]), 5),
        FALSE
      )
    ) %>%
    ungroup() %>%
    mutate(
      type_detail = case_when(
        is_up_top ~ "Up-Top5",
        is_down_top ~ "Down-Top5",
        TRUE ~ type
      )
    )
  
  comp_data$type_detail <- factor(
    comp_data$type_detail,
    levels = c("Up-Top5", "Up", "Down", "Down-Top5")
  )
  
  comp_data$cluster <- factor(
    comp_data$cluster,
    levels = c("C1", "C2", "C3", "C4", "C5")
  )
  
  return(comp_data)
}


# Transcriptome-consistent aging signatures 
match_signature_to_features <- function(signature_genes, reference_features) {
  candidate_genes <- unique(c(
    signature_genes,
    str_to_title(tolower(signature_genes)),
    toupper(signature_genes)
  ))
  
  intersect(candidate_genes, reference_features)
}

aging_up_set1 <- c(
  "CXCL13", "SPP1", "APOE", "AXL", "CCL3", "MMP12", "S100A8", "TTR",
  "S100A9", "CD36", "CD74", "IFITM3", "H2-AA", "CYBB", "CXCR4", "CLEC7A",
  "RAI14", "LILRB4", "FXYD5", "5430435G22RIK", "IFITM2", "H2-EB1", "FGL2",
  "CST7", "LTF", "CD44", "CYFIP2", "CTSE", "CAMP", "LGALS3", "OASL2",
  "H2-AB1", "RSAD2", "ST8SIA6", "CD274", "LGALS3BP", "ANXA2", "IGF2R",
  "LCN2", "EGR2", "LY9", "H2-M2", "UBE2E2", "NGP", "IQGAP1", "IL1R1",
  "DUSP1", "PCDHB22", "CXCL16", "USP18"
)

aging_up_set2 <- c(
  "Ccl4", "Lgals3", "Il1b", "Lpl", "Fam20c", "Cst7", "Csf1",
  "Ifitm3", "Ifit3", "Rtp4", "Irf7", "Isg15", "Oasl2",
  "F13a1", "Mrc1", "Pf4", "Clec12a", "Ms4a7", "Mgl2",
  "Ifit2", "Bcl2a1a", "Bcl2a1d", "Cdkn1a"
)

aging_up_set3 <- c(
  "Pdcd1", "Lyz2", "H2-K1", "Cd274", "Ifit1", "Slc2a4", "Tnf", "C3",
  "Neat1", "Cd52", "Pmch", "Ccl4", "Bst2", "Tac2", "Ifitm3", "Cdkn2a",
  "Cd79a", "Il2ra", "Stat1", "Prox1", "Gzmb", "Cxcl16", "Tpm4", "H2-D1",
  "B2m", "Cldn11", "Mcm2", "Crh", "S100a6", "Fgf2", "Klrd1", "Cxcl1",
  "Klk6", "Ccr7", "Pisd", "Cdkn1a", "Ly6c1", "Ifi27"
)

aging_down_set1 <- c(
  "PNMA2", "SLC15A2", "PPM1L", "PECAM1", "SNTA1", "STAU2", "ABHD6",
  "CDKN2AIP", "NRARP", "5730409E04RIK", "HIST1H2BM", "RGL1", "CLPTM1",
  "NCK1", "KCNRG", "STARD8", "FBXO30", "SUGT1", "DNAJB1", "ANKS6",
  "SLC24A3", "GADD45G", "MOCS2", "VEZF1", "IL7R", "TMC7", "NUP88",
  "CCDC90B", "USP38", "SYT12", "CD276", "CMTM7", "SLC35A2", "NRIP1",
  "GALNT10", "TSPAN33", "VPS26B", "HIST1H2AF", "DLGAP4", "TUBB2A",
  "TNFRSF17", "BBS9", "HIST1H2AN", "WAS", "CCDC12", "THTPA", "MSRB2",
  "0610037L13RIK", "ACADS"
)

aging_down_set2 <- c(
  "Gamt", "Ifngr1", "H2afx", "Dlx6os1", "Itgam", "Gpr17", "Snhg11",
  "Il7r", "Lsamp", "Insm1", "Ppargc1a", "Alcam", "Hexb", "Cd24a",
  "Hes5", "Cx3cr1", "Nts", "Pecam1", "Dcx", "Apln", "Fcrls", "P2ry12"
)

aging_down_set3 <- c(
  "P2ry12", "Tmem119", "Aif1", "Mertk", "Cx3cr1", "Hexb", "Csf1r", "P2ry13"
)

aging_up_merge <- unique(c(aging_up_set1, aging_up_set2, aging_up_set3))
aging_down_merge <- unique(c(aging_down_set1, aging_down_set2, aging_down_set3))

aging_up_protein <- match_signature_to_features(aging_up_merge, rownames(seu))
aging_down_protein <- match_signature_to_features(aging_down_merge, rownames(seu))


# Add transcriptome-consistency annotation for plotting

annotate_transcriptome_consistency <- function(comp_data,
                                               aging_up_protein,
                                               aging_down_protein) {
  comp_data %>%
    mutate(
      transcriptome_consistency = case_when(
        protein %in% aging_up_protein & avg_log2FC > 0 ~ "mRNA-level Up",
        protein %in% aging_down_protein & avg_log2FC < 0 ~ "mRNA-level Down",
        protein %in% aging_up_protein & avg_log2FC < 0 ~ "RNA-protein discordant",
        protein %in% aging_down_protein & avg_log2FC > 0 ~ "RNA-protein discordant",
        TRUE ~ "No consistency"
      ),
      plot_group = case_when(
        transcriptome_consistency == "mRNA-level Up" ~ "mRNA-level Up",
        transcriptome_consistency == "mRNA-level Down" ~ "mRNA-level Down",
        type == "Up" ~ "protein-level Up",
        type == "Down" ~ "protein-level Down",
        TRUE ~ "Other"
      ),
      label_color = case_when(
        transcriptome_consistency == "mRNA-level Up" ~ "#cd3333",
        transcriptome_consistency == "mRNA-level Down" ~ "#000080",
        TRUE ~ "black"
      )
    )
}


# Scatter plot with transcriptome-consistent labels 
create_transcriptome_consistency_scatter <- function(markers_data,
                                                     comparison_label,
                                                     cluster_colors = c(
                                                       "C1" = "#EE8172",
                                                       "C2" = "#7DBFA7",
                                                       "C3" = "#AECDE1",
                                                       "C4" = "#EE934E",
                                                       "C5" = "#A5917F"
                                                     ),
                                                     jitter_width = 0.35,
                                                     seed = 42) {
  markers_data <- markers_data %>%
    mutate(
      cluster = factor(cluster, levels = c("C1", "C2", "C3", "C4", "C5"))
    ) %>%
    arrange(cluster, protein)
  
  set.seed(seed)
  markers_data <- markers_data %>%
    mutate(
      cluster_x = as.numeric(cluster),
      x_jitter = cluster_x + runif(n(), -jitter_width, jitter_width)
    )
  
  cluster_levels <- levels(markers_data$cluster)
  cluster_levels <- cluster_levels[cluster_levels %in% unique(as.character(markers_data$cluster))]
  
  bg_df <- tibble(
    cluster = factor(cluster_levels, levels = levels(markers_data$cluster)),
    cluster_x = as.numeric(factor(cluster_levels, levels = levels(markers_data$cluster))),
    xmin = cluster_x - 0.37,
    xmax = cluster_x + 0.37,
    ymin = -0.5,
    ymax = 0.5
  )
  
  label_df <- markers_data %>%
    filter(transcriptome_consistency %in% c("mRNA-level Up", "mRNA-level Down"))
  
  p <- ggplot(markers_data, aes(x = x_jitter, y = avg_log2FC)) +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = bg_df,
      aes(x = cluster_x, y = 0, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      color = "black",
      show.legend = FALSE
    ) +
    geom_point(
      aes(color = plot_group, size = plot_group, alpha = plot_group),
      stroke = 0
    ) +
    geom_text_repel(
      data = label_df,
      aes(
        x = x_jitter,
        y = avg_log2FC,
        label = toupper(protein),
        color = plot_group
      ),
      size = 2.5,
      segment.color = "gray50",
      segment.size = 0.4,
      segment.alpha = 0.7,
      min.segment.length = 0.1,
      box.padding = 0.4,
      point.padding = 0.3,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    scale_color_manual(
      values = c(
        "protein-level Up" = "#eebbbb",
        "mRNA-level Up" = "#cd3333",
        "protein-level Down" = "#aaaad4",
        "mRNA-level Down" = "#000080",
        "Other" = "gray70"
      ),
      breaks = c(
        "protein-level Up",
        "mRNA-level Up",
        "protein-level Down",
        "mRNA-level Down"
      ),
      labels = c(
        "protein-level Up",
        "mRNA-level Up",
        "protein-level Down",
        "mRNA-level Down"
      ),
      name = ""
    ) +
    scale_size_manual(
      values = c(
        "protein-level Up" = 1,
        "mRNA-level Up" = 1.2,
        "protein-level Down" = 1,
        "mRNA-level Down" = 1.2,
        "Other" = 1
      ),
      guide = "none"
    ) +
    scale_alpha_manual(
      values = c(
        "protein-level Up" = 0.6,
        "mRNA-level Up" = 0.9,
        "protein-level Down" = 0.6,
        "mRNA-level Down" = 0.9,
        "Other" = 0.6
      ),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq_along(levels(markers_data$cluster)),
      labels = levels(markers_data$cluster),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_y_continuous(
      limits = c(
        min(markers_data$avg_log2FC, na.rm = TRUE) - 1,
        max(markers_data$avg_log2FC, na.rm = TRUE) + 1
      )
    ) +
    labs(
      x = "",
      y = paste0("log2FC (", comparison_label, ")"),
      title = paste("Age-differentially Expressed Proteins:", comparison_label)
    ) +
    guides(color = guide_legend(override.aes = list(size = 2, shape = 19))) +
    theme_few() +
    theme(
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      legend.position = c(0.82, 0.93),
      legend.key = element_blank(),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.text = element_text(size = 9, margin = margin(l = -1.5)),
      plot.title = element_text(size = 11, hjust = 0.5)
    )
  
  return(p)
}

# Bar plot for significant DEP counts 
create_dep_count_barplot <- function(markers_data, comparison_label) {
  
  count_data <- markers_data %>%
    group_by(cluster, type) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(cluster, type = c("Up", "Down"), fill = list(n = 0)) %>%
    mutate(
      type = factor(type, levels = c("Up", "Down")),
      y = ifelse(type == "Up", n, -n)
    )
  
  p <- ggplot(count_data, aes(x = cluster, y = y, fill = type)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
    scale_fill_manual(values = c("Up" = "#eebbbb", "Down" = "#aaaad4"), name = "") +
    scale_y_continuous(labels = abs) +
    labs(
      title = paste("DEP Counts:", comparison_label),
      x = "",
      y = "DEP counts"
    ) +
    theme_few() +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5),
      legend.position = c(0.85, 0.9),
      legend.background = element_blank(),
      legend.key.size = unit(0.3, "cm")
    )
  
  return(p)
}

# Run analysis and generate plots 
run_age_dep_analysis <- function(seu_obj,
                                 comparisons = list(
                                   c("24M", "2M"),
                                   c("24M", "14M"),
                                   c("14M", "2M")
                                 ),
                                 aging_up_protein,
                                 aging_down_protein) {
  
  cat("=== Age-dependent DEP analysis ===\n")
  
  all_dep_results <- compute_pairwise_dep_by_cluster(
    seu_obj = seu_obj,
    comparisons = comparisons
  )
  
  if (is.null(all_dep_results)) {
    stop("No differential proteins were found. Analysis stopped.")
  }
  
  comparison_names <- sapply(comparisons, function(x) paste0(x[1], " vs. ", x[2]))
  
  scatter_plots <- list()
  count_plots <- list()
  plot_data <- list()
  
  for (comp_name in comparison_names) {
    cat("\n--- Processing:", comp_name, "---\n")
    
    comp_data <- prepare_comparison_data(all_dep_results, comp_name)
    
    if (is.null(comp_data) || nrow(comp_data) == 0) {
      cat("  No significant differential proteins; skipped.\n")
      next
    }
    
    comp_data <- annotate_transcriptome_consistency(
      comp_data = comp_data,
      aging_up_protein = aging_up_protein,
      aging_down_protein = aging_down_protein
    )
    
    consistent_stats <- comp_data %>%
      filter(transcriptome_consistency %in% c("mRNA-level Up", "mRNA-level Down")) %>%
      group_by(cluster, transcriptome_consistency) %>%
      summarise(n = n(), .groups = "drop")
    
    cat("  Transcriptome-consistent proteins:\n")
    print(consistent_stats)
    
    scatter_plots[[comp_name]] <- create_transcriptome_consistency_scatter(
      markers_data = comp_data,
      comparison_label = comp_name
    )
    
    count_plots[[comp_name]] <- create_dep_count_barplot(
      markers_data = comp_data,
      comparison_label = comp_name
    )
    
    plot_data[[comp_name]] <- comp_data
  }
  
  return(list(
    all_results = all_dep_results,
    plot_data = plot_data,
    scatter_plots = scatter_plots,
    count_plots = count_plots
  ))
}

### Execute analysis 
age_comparisons <- list(
  c("24M", "2M"),
  c("24M", "14M"),
  c("14M", "2M")
)

age_dep_output <- run_age_dep_analysis(
  seu_obj = seu,
  comparisons = age_comparisons,
  aging_up_protein = aging_up_protein,
  aging_down_protein = aging_down_protein
)

all_results <- age_dep_output$all_results

# Plot objects
scatter_plots <- age_dep_output$scatter_plots
count_plots <- age_dep_output$count_plots

# Example:
scatter_plots[["24M vs. 2M"]]
count_plots[["24M vs. 2M"]]


