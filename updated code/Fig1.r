# Fig1.R
# Purpose: Comparative analysis with published scRNA-seq datasets.

# 0.Load Packages -------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(Seurat)
library(qs)
library(Matrix)
library(ggthemes)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(eulerr)

# 0.Load Files -------------------------------------------------------------------------------
# SCP
scp_seu <- readRDS("res/SCP_MG_3085_seurat_post_QC.rds")
DefaultAssay(scp_seu) <- "SCP"
scp_count <- as.data.frame(GetAssayData(scp_seu, layer = "counts"))
scp_data <- as.data.frame(LogNormalize(data = scp_count))

# scRNA-seq datasets
# Hammond et al.
hammond_all_samples <- qread("data/public_RNA/Hammond/Hammond_et-al-2019_Seurat_Converted_v4.qs")
rna1_seu <- hammond_all_samples[rowSums(hammond_all_samples@assays$RNA@counts) > 0, ] 
rna1_count <- as.data.frame(GetAssayData(hammond_all_samples, layer = "counts"))

# Masuda et al.
rna2_count <- fread("data/public_RNA/Masuda/GSE120745_matrix.txt") %>%
  column_to_rownames("V1")
rna2_count <- rna2_count[rowSums(rna2_count == 0) != ncol(rna2_count), ] 

# Li et al.
rna3_count <- fread("data/public_RNA/Li/GSE123025_matrix.txt") %>%
  column_to_rownames("V1") %>%  
  mutate(across(everything(), ~ replace_na(., 0)))
rna3_count <- rna3_count[rowSums(rna3_count == 0) != ncol(rna3_count), ] 

# Jin et al.
mtx <- readMM("data/public_RNA/Jin/nature2025_mouse_microglia_BAM_7w.mtx")
cl <- fread("data/public_RNA/Jin/nature2025_mouse_microglia_BAM_cellinfo.csv", header = T, data.table = F) 
rl <- fread("data/public_RNA/Jin/nature2025_mouse_microglia_BAM_geneinfo.csv", header = T, data.table = F) 
rownames(mtx) <- cl$sample_id
colnames(mtx) <- rl$V1
rna4_count <- as.data.frame(t(mtx))


coding_gene <- fread("data/public_RNA/genecode.mm10_protein.coding.txt")
 
full_count <- list( 
  `Hammond et al.` = rna1_count,
  `Masuda et al.` = rna2_count,  
  `Li et al.` = rna3_count,
  `Jin et al.` = rna4_count
)  
# set the threshold 
min_genes_per_cell <- 650
min_expression_cells <- 1

full_data_filterd <- list()
for (dataset_name in names(full_count)) {
  adata <- full_count[[dataset_name]]
  
  gene_expression_matrix <- adata %>%
    filter(rownames(adata) %in% coding_gene$V5)
  
  # cell-level QC
  genes_detected_per_cell <- colSums(gene_expression_matrix != 0)
  filtered_cells_indices <- which(genes_detected_per_cell >= min_genes_per_cell)
  filtered_cells <- gene_expression_matrix[, filtered_cells_indices]
  
  # gene-level QC
  cells_detected_per_gene <- rowSums(filtered_cells != 0)
  filtered_genes <- filtered_cells[which(cells_detected_per_gene >= min_expression_cells), ]
  
  print(dim(filtered_genes))
  
  rna_data <- as.data.frame(LogNormalize(data = filtered_genes))
  full_data_filterd[[dataset_name]] <- rna_data
}

full_data_filterd$SCP <- scp_data
matrix_list <- full_data_filterd


# 1.SCP Saturation Curve -----------------------------------------------------------------------------------
set.seed(42)   
expression_matrix <- scp_count[, sample(ncol(scp_count))] 
max_cells <- ncol(expression_matrix)  

saturation_results <- data.frame(Cell_Number = integer(), 
                                 Identified_Proteins = integer())  

for (cell_count in 1:max_cells) {  
  sampled_cells <- data.frame(expression_matrix[, 1:cell_count]) 
  identified_proteins <- sum(rowSums(sampled_cells) > 0)  
  saturation_results <- rbind(saturation_results, data.frame(Cell_Number = cell_count, 
                                                             Identified_Proteins = identified_proteins))  
}  

index <- head(which(saturation_results$Identified_Proteins > nrow(scp_count)*0.95), n = 1)  
point <- saturation_results$Cell_Number[index]  
print(point)
# 661
p_curve <- ggplot(saturation_results, aes(x = Cell_Number, y = Identified_Proteins)) +  
  geom_line(color = "gray", linewidth = 0.5) +  
  geom_point(color = "#da6968", size = 0.5) +  
  geom_vline(xintercept = saturation_results$Cell_Number[point], linetype = "dashed", color = "black", alpha = 0.6) +  
  labs(title = "Protein Saturation Curve",  
       x = "Number of Cells",  
       y = "Number of Proteins") +  
  theme_few()  


# 2.Protein completeness -----------------------------------------------------------------------------------
completeness_df <- data.frame(completeness = rowSums(scp_count > 0) / ncol(scp_count) * 100) 
p_com <- ggplot(completeness_df, aes(x = completeness)) +  
  geom_histogram(alpha = 0.8,linewidth=0.2,color="black",binwidth = 5,fill = colorRampPalette(c("#f8e6e5", "#da6968"))(21)) + 
  labs(title = "Data completeness across 3085 cells",  
       x = "Completeness (%)",  
       y = "Number of proteins") +  
  theme_few() +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
sum(completeness_df$completeness > 75)
# 588


# 3.Comparative analysis with Transcriptome -------------------------------------------------------------------------------
## 3.1 Overlap between SCP and RNA -------------------------------------------------------------------------------
proteins <- rownames(matrix_list[[5]])
gene_list <- list(rownames(matrix_list[[1]]), 
                  rownames(matrix_list[[2]]), 
                  rownames(matrix_list[[3]]), 
                  rownames(matrix_list[[4]]))  

overlap_df <- data.frame(  
  source = c("Hammond et al.","Masuda et al.","Li et al.","Jin et al."),  
  shared = c(length(intersect(proteins, gene_list[[1]])),  
             length(intersect(proteins, gene_list[[2]])),  
             length(intersect(proteins, gene_list[[3]])),  
             length(intersect(proteins, gene_list[[4]]))) 
)  

overlap_df <- overlap_df %>%  
  mutate(unique = length(proteins) - shared) %>%  
  pivot_longer(cols = c(shared, unique),   
               names_to = "Overlap",   
               values_to = "Count") %>%  
  mutate(source = factor(source, levels = c("Masuda et al.","Li et al.","Hammond et al.","Jin et al."))) %>%  
  mutate(Overlap = factor(Overlap, levels = c("unique","shared")))
 
p_overlap <- ggplot(overlap_df, aes(x = source, y = Count, fill = Overlap)) +  
  geom_bar(stat="identity",position="fill",width=0.6) +
  scale_fill_manual(values = c("gray","black")) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Overlap between SCP\nand Transcriptome",  
       x = "",  
       y = "Overlap ratio") +   
  theme_cowplot(font_size = 22) +
  theme(plot.title = element_text(hjust = 0.5,size = 13),
        axis.title = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length = unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "none",
        legend.position.inside = c(0,1.1),
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3,'cm'),
        legend.text = element_text(size = 10, margin = margin(t = -2, r = 5, b = 0, l = 0)))

for (dataset_name in names(matrix_list)[1:4]) {  
  adata <- matrix_list[[dataset_name]]
  genes <- rownames(adata)
  overlap <- intersect(proteins,genes)   
  df <- c("Proteome" = length(proteins)-length(overlap),
          "Transcriptome" = length(genes)-length(overlap), 
          "Proteome&Transcriptome" = length(overlap))
  p <- plot(euler(df), 
            quantities = list(type = "counts",cex=1),          
            edges = list(col = "lightgray", lex = 2, lwd=0.6), 
            fills = list(fill = c("#2e3c65", "#77a8bc"),alpha=0.7),
            main = list(label=dataset_name,cex=1.5))
  print(p)
}


##---3.2 Protein/gene identification depth (number of genes detected in each cell) ---------------------------------
result_count <- list()
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]
  result_df <- data.frame(cell = colnames(adata), 
                          gene_count = colSums(adata > 0),
                          source = dataset_name)
  result_count[[dataset_name]] <- result_df  
}

combined_count <- do.call(rbind, result_count) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "Masuda et al.",
                                           "Li et al.",
                                           "Hammond et al.",
                                           "Jin et al.")))

p_count <- ggplot(combined_count, aes(x = source, y = gene_count, fill=source)) +
  geom_violin(linewidth = 0.4) +
  geom_boxplot(fill = "white",
               size=0.35,outlier.size=0.8,width=0.15,outlier.shape = NA, median.linewidth = 0.4) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#B2B6C1","#9196A0","#71757F")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Gene/Protein count per cell", title = "") +  
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank())
median_by_source <- combined_count %>%
  group_by(source) %>%
  summarise(median_gene_count = median(gene_count, na.rm = TRUE),
            mean_gene_count = mean(gene_count, na.rm = TRUE))


## 3.3 Expression completeness per cell -------------------------------------------------------------------------------
result_completeness <- list()
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]
  result_df <- data.frame(cell = colnames(adata), 
                          completeness = apply(adata, MARGIN = 2,   
                                               FUN = function(x) sum(x != 0, na.rm = TRUE) / length(x[!is.na(x)])),
                          source = dataset_name)
  result_completeness[[dataset_name]] <- result_df  
}

combined_completeness <- do.call(rbind, result_completeness) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "Masuda et al.",
                                           "Li et al.",
                                           "Hammond et al.",
                                           "Jin et al.")))

p_completeness <- ggplot(combined_completeness, aes(x = source, y = completeness, fill=source)) +
  stat_boxplot(geom ="errorbar",size = 0.5,width=0.25,position = position_dodge(width=0.5)) +
  geom_boxplot(size=0.4,outlier.size=0.8,width=0.5,outlier.shape = NA, median.linewidth = 0.6) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#B2B6C1","#9196A0","#71757F")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Data completeness", title = "") +  
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank())


## random 4349 genes
result_completeness_sampled <- list()

for (dataset_name in names(matrix_list)) {
  adata <- matrix_list[[dataset_name]]
  
  set.seed(42)  
  sampled_genes <- sample(rownames(adata), size = 4349, replace = FALSE)
  
  adata_sampled <- adata[rownames(adata) %in% sampled_genes, ]
  
  result_df <- data.frame(
    cell = colnames(adata_sampled),
    completeness = apply(adata_sampled, 2, 
                         function(x) sum(x != 0, na.rm = TRUE) / length(x[!is.na(x)])),
    source = dataset_name
  )
  
  result_completeness_sampled[[dataset_name]] <- result_df
}

combined_completeness_sampled <- do.call(rbind, result_completeness_sampled) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "Masuda et al.",
                                           "Li et al.",
                                           "Hammond et al.",
                                           "Jin et al.")))

p_completeness_sampled <- ggplot(combined_completeness_sampled, aes(x = source, y = completeness, fill = source)) +
  stat_boxplot(geom = "errorbar", size = 0.5, width = 0.25, position = position_dodge(width = 0.5)) +
  geom_boxplot(size = 0.4, outlier.size = 0.8, width = 0.5, outlier.shape = NA, median.linewidth = 0.6) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#B2B6C1","#9196A0","#71757F")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Data completeness\n(Random 4,349 genes)", title = "") +  
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11),
        axis.ticks.length = unit(0.1, 'cm'),
        legend.title = element_blank())


## top 1000 proteins/genes
result_completeness_top <- list()
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]
  
  gene_means <- rowMeans(expm1(adata))
  
  top_genes <- names(sort(gene_means, decreasing = TRUE))[1:1000]
  
  adata_top <- adata[rownames(adata) %in% top_genes, ]
  
  result_df <- data.frame(cell = colnames(adata_top), 
                          completeness = apply(adata_top, 2,   
                                               function(x) sum(x != 0, na.rm = TRUE) / length(x[!is.na(x)])),
                          source = dataset_name)
  result_completeness_top[[dataset_name]] <- result_df  
}
 
combined_completeness_top <- do.call(rbind, result_completeness_top) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "Masuda et al.",
                                           "Li et al.",
                                           "Hammond et al.",
                                           "Jin et al.")))

p_completeness_top <- ggplot(combined_completeness_top, aes(x = source, y = completeness, fill=source)) +
  stat_boxplot(geom ="errorbar",size = 0.5,width=0.25,position = position_dodge(width=0.5)) +
  geom_boxplot(size=0.4,outlier.size=0.8,width=0.5,outlier.shape = NA, median.linewidth = 0.6) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#B2B6C1","#9196A0","#71757F")) +
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Data completeness\n(Top 1,000 proteins/genes)", title = "") +  
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank())


## shared high-abundance genes
get_top_genes <- function(expr_matrix, quantile_thresh = 0.5) {
  gene_means <- rowMeans(expm1(expr_matrix))
  cutoff <- quantile(gene_means, probs = 1 - quantile_thresh, na.rm = TRUE)
  return(names(gene_means)[gene_means >= cutoff])
}

protein_genes <- get_top_genes(matrix_list[[5]])

rna_gene_sets <- list(
  `Hammond et al.` = get_top_genes(matrix_list[[1]]),
  `Masuda et al.` = get_top_genes(matrix_list[[2]]),
  `Li et al.` = get_top_genes(matrix_list[[3]]),
  `Jin et al.` = get_top_genes(matrix_list[[4]])
)
gene_overlap <- reshape2::melt(rna_gene_sets) %>% 
  group_by(value) %>% 
  summarise(n_datasets = n_distinct(L1))
required_datasets <- 4  
shared_rna_genes <- gene_overlap %>% 
  filter(n_datasets >= required_datasets) %>% 
  pull(value)

shared_genes <- intersect(shared_rna_genes, protein_genes)


result_completeness_shared <- list()
for (dataset_name in names(matrix_list)) {
  adata <- matrix_list[[dataset_name]]
  
  adata_shared <- adata[rownames(adata) %in% shared_genes, ]
  
  result_df <- data.frame(
    cell = colnames(adata_shared),
    completeness = apply(adata_shared, 2, 
                         function(x) sum(x != 0, na.rm = TRUE) / nrow(adata_shared)),
    source = dataset_name
  )
  result_completeness_shared[[dataset_name]] <- result_df
}
 
combined_completeness_shared <- do.call(rbind, result_completeness_shared) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "Masuda et al.",
                                           "Li et al.",
                                           "Hammond et al.",
                                           "Jin et al.")))

p_completeness_shared <- ggplot(combined_completeness_shared, aes(x = source, y = completeness, fill=source)) +
  stat_boxplot(geom ="errorbar",size=0.5,width=0.25,position = position_dodge(width=0.5)) +
  geom_boxplot(size=0.4,outlier.size=0.8,width=0.5,outlier.shape = NA, median.linewidth = 0.6) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#B2B6C1","#9196A0","#71757F")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Data completeness\n(Shared high-abundance genes)", title = "") +  
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank())


##---3.4 Heat map of cell–cell correlations ---------------------------------
shared_genes <- Reduce(intersect, list(rownames(matrix_list[[1]]), 
                                       rownames(matrix_list[[2]]), 
                                       rownames(matrix_list[[3]]),
                                       rownames(matrix_list[[4]]),
                                       rownames(matrix_list[[5]])))  
data_filtered <- list(
  SCP = matrix_list[[5]][shared_genes, ],
  RNA1 = matrix_list[[2]][shared_genes, ],
  RNA2 = matrix_list[[3]][shared_genes, ],  
  RNA3 = matrix_list[[1]][shared_genes, ],  
  RNA4 = matrix_list[[4]][shared_genes, ]
)

sample_size <- 500
sampled_data <- list()  
for (name in names(data_filtered)) {  
  adata <- data_filtered[[name]]  
  set.seed(42)
  sampled_columns <- sample(ncol(adata), sample_size)  
  sampled_data[[name]] <- adata[, sampled_columns]  
}

cor_results <- lapply(sampled_data, function(sampled) {
  cor(sampled, method = "pearson") 
})
 
mean_values <- lapply(cor_results, mean)

joined <- bind_cols(sampled_data)  
cor <- cor(joined, method = "pearson") 

p_cor <- ComplexHeatmap::pheatmap(cor,   
                                  show_rownames = F,   
                                  show_colnames = F,  
                                  cluster_rows = F, 
                                  cluster_cols = F,
                                  gaps_col = c(500,1000,1500,2000),
                                  gaps_row = c(500,1000,1500,2000),
                                  border_gp = gpar(col = "black"),
                                  color = YlGnBu_gradient(100),  
                                  use_raster = FALSE,
                                  fontsize = 10,
                                  heatmap_legend_param = list(title = "Pearson correlation"))


# 4.Protein-to-mRNA abundance plot -----------------------------------------------------------------------------------
scp_seu <- readRDS("res/SCP_MG_3085_seurat_post_QC.rds")
DefaultAssay(scp_seu) <- "SCP"
scp_count <- as.data.frame(GetAssayData(scp_seu, layer = "counts"))

hammond_all_samples <- qread("data/public_RNA/Hammond/Hammond_et-al-2019_Seurat_Converted_v4.qs")
rna1_seu_age <- hammond_all_samples %>% 
  subset(Age %in% c("P30", "P100", "P540")) %>%  
  subset(Sex == "Male")
rna1_seu_age <- rna1_seu_age[rowSums(rna1_seu_age@assays$RNA@counts) > 0, ] 
rna1_count <- as.data.frame(GetAssayData(rna1_seu_age, layer = "counts"))

common_genes <- intersect(rownames(scp_count), rownames(rna1_count))
scp_common <- scp_count[common_genes, ]
rna1_common <- rna1_count[common_genes, ]
identical(rownames(scp_common), rownames(rna1_common))
seurat_obj <- CreateSeuratObject(counts = cbind(scp_common, rna1_common))

data_types <- c(rep("SCP", ncol(scp_common)), rep("RNA", ncol(rna1_common)))  
seurat_obj <- AddMetaData(seurat_obj, metadata = data_types, col.name = "type")  
head(seurat_obj@meta.data) 

seurat_obj <- NormalizeData(seurat_obj)

exp_data <- as.data.frame(AverageExpression(seurat_obj, group.by = c("type"), slot = "data"))
head(exp_data)  
colnames(exp_data) <- c("RNA","SCP")

correlation <- cor(exp_data$RNA, 
                   exp_data$SCP, 
                   method = "spearman")  

p_cor <- ggplot(exp_data, aes(x = log(RNA+1), y = log(SCP+1)))+
  geom_point(size=0.8, alpha=0.6)+
  ggpubr::stat_cor(method="spearman",label.x = 0.1)+
  labs(x = "RNA Expression", y = "Protein Abundance") +
  theme_few()


## Calculation of fold-change and p-value
rna1_common <- as.matrix(rna1_common)
scp_common <- as.matrix(scp_common)
result_fc <- exp_data %>%
  mutate(
    fold_change = (RNA+1)/(SCP+1),  
    log2_fold_change = log2(fold_change),
    p_value = sapply(1:nrow(scp_common), function(i) {
      wilcox.test(rna1_common[i, ], scp_common[i, ])$p.value
    }),
    p_adj = p.adjust(p_value, method = "fdr")
  ) 

plot <- result_fc %>%
  mutate(significance = case_when(
    p_adj < 0.05 & 
      log2_fold_change > 2 ~ "Lower in protein",
    p_adj < 0.05 & 
      log2_fold_change < -2 ~ "Higher in protein",
    TRUE ~ "NoSig"
  ))
table(plot$significance)
plot$Protein <- rownames(plot)

p_cor1 <- ggplot(plot, aes(x = log(RNA+1), y = log(SCP+1))) +
  geom_point(aes(color = significance), size = 1, alpha=0.5) +
  scale_color_manual(values = c("Higher in protein" = "#d35230",
                                "Lower in protein" = "#2b7cd3",
                                "NoSig" = "gray")) +
  ggrepel::geom_text_repel(data = subset(plot,significance!="NoSig"),
                           aes(label = Protein),
                           size = 3, 
                           box.padding = 0.5, 
                           point.padding = 0.8,
                           min.segment.length = 0.5, 
                           segment.color = "black",
                           show.legend = F) +
  ggpubr::stat_cor(method="spearman", label.x = 4) +
  geom_abline(intercept = c(-log(4), 0, log(4)), slope = 1,   # |log2FC|>2
              linetype = c("dotted", "dashed", "dotted"),
              color = c("gray", "black", "gray") ,linewidth = 0.5) +  
  scale_x_continuous(limits = c(0,6.6)) +
  scale_y_continuous(limits = c(0,6.6)) +
  labs(x = "Mean log normalized mRNA expression", 
       y = "Mean log normalized protein abundance") +
  theme_few() +
  theme(legend.position = 'none')


