# 0.Load Packages -------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(Seurat)
library(scCustomize)
library(qs)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(cowplot)
library(eulerr)
source("code/scp_utils.R")

# 0.Load Files -------------------------------------------------------------------------------
# SCP
scp_seu <- readRDS("data/SCP_MG_Seurat.rds")
DefaultAssay(scp_seu) <- "SCP"
scp_count <- as.data.frame(GetAssayData(scp_seu, layer = "counts"))
scp_data <- as.data.frame(LogNormalize(data = scp_count))

# scRNA-seq1
hammond_all_samples <- qread("data/scRNA-seq/Hammond/Hammond_et-al-2019_Seurat_Converted_v4.qs")
rna1_seu <- hammond_all_samples[rowSums(hammond_all_samples@assays$RNA@counts) > 0, ] 
rna1_count <- as.data.frame(GetAssayData(rna1_seu, layer = "counts"))

# scRNA-seq2
rna2_count <- fread("data/scRNA-seq/Masuda/GSE120745_matrix.txt") %>%
  column_to_rownames("V1")
rna2_count <- rna2_count[rowSums(rna2_count == 0) != ncol(rna2_count), ] 

# scRNA-seq3
rna3_count <- fread("data/scRNA-seq/Li/GSE123025_matrix.txt") %>%
  column_to_rownames("V1") %>%  
  mutate(across(everything(), ~ replace_na(., 0)))
rna3_count <- rna3_count[rowSums(rna3_count == 0) != ncol(rna3_count), ] 

full_count <- list( 
  `2019 Immunity` = rna1_count,
  `2019 Nature` = rna2_count,  
  `2019 Neuron` = rna3_count 
)  

coding_gene <- fread("data/scRNA-seq/genecode.mm10.protein.coding.txt")

min_genes_per_cell <- 650
min_expression_cells <- 1

full_data_filterd <- list()
for (dataset_name in names(full_count)) {
  adata <- full_count[[dataset_name]]
  
  gene_expression_matrix <- adata %>%
    filter(rownames(adata) %in% coding_gene$V5)
  
  genes_detected_per_cell <- colSums(gene_expression_matrix != 0)
  filtered_cells_indices <- which(genes_detected_per_cell >= min_genes_per_cell)
  filtered_cells <- gene_expression_matrix[, filtered_cells_indices]
  
  cells_detected_per_gene <- rowSums(filtered_cells != 0)
  filtered_genes <- filtered_cells[which(cells_detected_per_gene >= min_expression_cells), ]
  
  print(dim(filtered_genes))
  
  rna_data <- as.data.frame(LogNormalize(data = filtered_genes))
  full_data_filterd[[dataset_name]] <- rna_data
}
data_list <- full_data_filterd
data_list$SCP <- scp_data
matrix_list <- data_list

# 1.overlap between SCP and RNA -------------------------------------------------------------------------------
proteins <- rownames(matrix_list[[4]])
gene_list <- list(rownames(matrix_list[[1]]), 
                  rownames(matrix_list[[2]]), 
                  rownames(matrix_list[[3]]))  

overlap_df <- data.frame(  
  source = c("10x Genomics","Smart-seq2 ①","Smart-seq2 ②"),  
  shared = c(length(intersect(proteins, gene_list[[1]])),  
             length(intersect(proteins, gene_list[[2]])),  
             length(intersect(proteins, gene_list[[3]]))) 
)  

overlap_df <- overlap_df %>%  
  mutate(unique = length(proteins) - shared) %>%  
  pivot_longer(cols = c(shared, unique),   
               names_to = "Overlap",   
               values_to = "Count") %>%  
  mutate(Overlap=factor(Overlap,levels=c("unique","shared")))

p_overlap <- ggplot(overlap_df, aes(x = source, y = Count, fill = Overlap)) +  
  geom_bar(stat="identity",position="fill",width=0.6) +
  scale_fill_manual(values = c("gray","black")) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "",  
       x = "",  
       y = "Overlap ratio between proteins\nand three gene datasets") +   
  theme_cowplot(font_size = 22) +
  theme(axis.title = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length = unit(0.1,'cm'),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0,1.1),
        legend.direction = "horizontal", 
        legend.key.size = unit(0.3,'cm'),
        legend.text = element_text(size = 10, margin = margin(t = -2, r = 5, b = 0, l = 0)))
  

# 2.expression completeness per cell -------------------------------------------------------------------------------
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
  mutate(source = plyr::mapvalues(from = c("2019 Immunity",
                                           "2019 Nature",
                                           "2019 Neuron"),
                                  to = c("10x Genomics",
                                         "Smart-seq2 ①",
                                         "Smart-seq2 ②"),
                                  source)) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "10x Genomics",
                                           "Smart-seq2 ①",
                                           "Smart-seq2 ②")))

p_completeness <- ggplot(combined_completeness, aes(x = source, y = completeness, fill=source)) +
  stat_boxplot(geom ="errorbar",width=0.25,position = position_dodge(width=0.5)) +
  geom_boxplot(size=0.4,outlier.size=0.8,width=0.5,outlier.shape = NA,fatten=1.2) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#a8aebc","#7e828d")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Gene/Protein expression\ncompleteness per cell", title = "") +  
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'),
        legend.title = element_blank())


# Protein/gene identification depth (number of genes detected in each cell)
result_count <- list()
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]
  result_df <- data.frame(cell = colnames(adata), 
                          gene_count = colSums(adata > 0),
                          source = dataset_name)
  result_count[[dataset_name]] <- result_df  
}

combined_count <- do.call(rbind, result_count) %>%
  mutate(source = plyr::mapvalues(from = c("2019 Immunity",
                                           "2019 Nature",
                                           "2019 Neuron"),
                                  to = c("10x Genomics",
                                         "Smart-seq2 ①",
                                         "Smart-seq2 ②"),
                                  source)) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "10x Genomics",
                                           "Smart-seq2 ①",
                                           "Smart-seq2 ②")))

p_count <- ggplot(combined_count, aes(x = source, y = gene_count, fill=source)) +
  geom_violin() +
  geom_boxplot(fill = "white",size=0.6,outlier.size=0.8,width=0.15,outlier.shape = NA,lwd= 0.5,fatten=1) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#a8aebc","#7e828d")) + 
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
  summarise(median_gene_count = median(gene_count, na.rm = TRUE))


#---3.Comparison of CV between SCP and scRNA-seq ---------------------------------
cal_cv <- function(x){  
  y <- na.omit(x)
  return(sd(y)/mean(y))
}

full_data_cv <- list() 
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]    
  adata.cv <- data.frame(cv = apply(adata, 1, cal_cv)) %>%   
    rownames_to_column() %>%  
    mutate(source = dataset_name)  
  full_data_cv[[dataset_name]] <- adata.cv   
}

combined.cv <- do.call(rbind, full_data_cv) %>%
  mutate(source = plyr::mapvalues(from = c("2019 Immunity",
                                           "2019 Nature",
                                           "2019 Neuron"),
                                  to = c("10x Genomics",
                                         "Smart-seq2 ①",
                                         "Smart-seq2 ②"),
                                  source)) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "10x Genomics",
                                           "Smart-seq2 ①",
                                           "Smart-seq2 ②")))

p_cv_all <- ggplot(combined.cv, aes(x=source, y=cv,fill=source)) +
  stat_boxplot(geom ="errorbar",width=0.25,position = position_dodge(width=0.5)) +
  geom_boxplot(size=0.4,outlier.size=0.8,width=0.5,outlier.shape = NA,lwd= 0.5,fatten=1) +
  scale_fill_manual(values = c("#da6968","#d3d7e2","#a8aebc","#7e828d")) + 
  theme_cowplot(font_size = 22) +
  theme(legend.position = 'none') +
  labs(x = "", y = "Coefficient of variation", title = "") +  
  coord_cartesian(ylim = c(0,150)) +
  theme(axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11),
        axis.ticks.length=unit(0.1,'cm'))


# 4.Downsampling analysis -----------------------------------------------------------------------------------
result_saturation <- data.frame(Cell_Number = integer(), 
                                Identified_Proteins = integer(),
                                Percentage = numeric())
for (dataset_name in names(matrix_list)) {  
  adata <- matrix_list[[dataset_name]]
  set.seed(111)   
  expression_matrix <- adata[, sample(ncol(adata))] 
  max_cells <- ncol(expression_matrix)  
  saturation_results <- data.frame(Cell_Number = integer(), 
                                   Identified_Proteins = integer(),
                                   Percentage = numeric())  
  for (cell_count in 1:max_cells) {  
    sampled_cells <- data.frame(expression_matrix[, 1:cell_count])  
    identified_proteins <- sum(rowSums(sampled_cells) > 0)  
    saturation_results <- rbind(saturation_results, data.frame(Cell_Number = cell_count, 
                                                               Identified_Proteins = identified_proteins,
                                                               Percentage = identified_proteins/nrow(expression_matrix)*100))
  }  
  saturation_results$source <- dataset_name
  result_saturation <- rbind(result_saturation,saturation_results)
}

result_saturation <- result_saturation %>%
  mutate(source = plyr::mapvalues(from = c("2019 Immunity",
                                           "2019 Nature",
                                           "2019 Neuron"),
                                  to = c("10x Genomics",
                                         "Smart-seq2 ①",
                                         "Smart-seq2 ②"),
                                  source)) %>%
  mutate(source = factor(source,levels = c("SCP",
                                           "10x Genomics",
                                           "Smart-seq2 ①",
                                           "Smart-seq2 ②")))
 
p_saturation <- ggplot(result_saturation, aes(x = Cell_Number, y = Percentage, color = source)) +  
  geom_line(linewidth = 0.3) +  
  geom_point(alpha = 0.6, size = 0.3) +
  scale_color_manual(values = c("#da6968","#d3d7e2","#a8aebc","#7e828d")) +  
  labs(
    title = "Protein/mRNA Saturation Curve", 
    x = "Number of Cells",  
    y = "Percentage of total proteins/genes (%)",
    color = "Dataset"
  ) +
  scale_x_sqrt() +
  theme_minimal() +  
  theme(legend.position = "inside",
        legend.position.inside = c(0.7,0.35),
        legend.text = element_text(size=10)) 


#---5.Heat map of cell–cell correlations ---------------------------------
shared_genes <- Reduce(intersect, list(rownames(matrix_list[[1]]), 
                                       rownames(matrix_list[[2]]), 
                                       rownames(matrix_list[[3]]),
                                       rownames(matrix_list[[4]])))  
data_filtered <- list(
  SCP = matrix_list[[4]][shared_genes, ],
  RNA1 = matrix_list[[1]][shared_genes, ],
  RNA2 = matrix_list[[2]][shared_genes, ],  
  RNA3 = matrix_list[[3]][shared_genes, ]
)

sample_size <- 500 
sampled_data <- list()  
for (name in names(data_filtered)) {  
  adata <- data_filtered[[name]] 
  set.seed(100)
  sampled_columns <- sample(ncol(adata), sample_size)  
  sampled_data[[name]] <- adata[, sampled_columns]  
}

cor_results <- lapply(sampled_data, function(sampled) {
  cor(sampled, method = "pearson") 
})
  
mean_values <- lapply(cor_results, mean)

joined <- bind_cols(sampled_data)   
c <- cor(joined, method = "pearson") 

p_cor <- ComplexHeatmap::pheatmap(c,    
                                  show_rownames = F,   
                                  show_colnames = F,  
                                  cluster_rows = F, 
                                  cluster_cols = F,
                                  gaps_col = c(500,1000,1500),
                                  gaps_row = c(500,1000,1500),
                                  border_gp = gpar(col = "black"),
                                  color = YlGnBu_gradient(100),  
                                  fontsize = 10,
                                  heatmap_legend_param = list(title = "Pearson correlation"))


# 6.Protein-to-mRNA abundance plot -----------------------------------------------------------------------------------
filtered_mdata <- rna1_seu@meta.data[rna1_seu@meta.data$Age %in% c("P30", "P100", "P540") & rna1_seu@meta.data$Sex == "Male", ]  
rna1_data_filter <- matrix_list[[1]][, colnames(matrix_list[[1]]) %in% rownames(filtered_mdata), drop = FALSE]  

common_genes <- intersect(rownames(matrix_list[[4]]), rownames(rna1_data_filter))  
scp_data_common <- matrix_list[[4]][common_genes, , drop = FALSE]
rna1_data_common <- rna1_data_filter[common_genes, , drop = FALSE]

average_expr_scp <- data.frame(Protein = rownames(scp_data_common),   
                               MeanExpression = rowMeans(scp_data_common))
average_expr_rna <- data.frame(Gene = rownames(rna1_data_common),   
                               MeanExpression = rowMeans(rna1_data_common))

merged_cor <- left_join(average_expr_scp, 
                        average_expr_rna,
                        by = c("Protein" = "Gene"),
                        suffix = c(".protein", ".rna"))

correlation <- cor(merged_cor$MeanExpression.protein, 
                   merged_cor$MeanExpression.rna, 
                   method = "spearman")  

#### Calculate fold change and p-value
rna1_data_common <- as.matrix(rna1_data_common)
scp_data_common <- as.matrix(scp_data_common)
set.seed(123) 
result_fc <- merged_cor %>%
  mutate(
    fold_change = exp(MeanExpression.rna - MeanExpression.protein),  
    log2_fold_change = log2(fold_change), 
    p_value = sapply(1:nrow(scp_data_common), function(i) {
      n_samples <- min(ncol(rna1_data_common), ncol(scp_data_common))
      rna_sample <- sample(rna1_data_common[i, ], n_samples, replace = FALSE)
      protein_sample <- sample(scp_data_common[i, ], n_samples, replace = FALSE)
      
      wilcox.test(rna_sample, protein_sample)$p.value
    }),
    p_adj = p.adjust(p_value, method = "fdr")
  ) 

plot <- result_fc %>%
  mutate(significance = case_when(
    p_adj < 0.05 & log2_fold_change >= 2 ~ "Lower in protein",
    p_adj < 0.05 & log2_fold_change <= -2 ~ "Higher in protein",
    TRUE ~ "NoSig"
  ))

p_cor <- ggplot(plot, aes(x = MeanExpression.rna, y = MeanExpression.protein)) +
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
  geom_abline(intercept = c(-log(4), 0, log(4)), slope = 1, 
              linetype = c("dotted", "dashed", "dotted"),
              color = c("gray", "black", "gray") ,linewidth = 0.5) + 
  scale_x_continuous(limits = c(0,6.6)) +
  scale_y_continuous(limits = c(0,6.6)) +
  labs(x = "Mean log normalized mRNA expression", 
       y = "Mean log normalized protein abundance") +
  theme_few() +
  theme(legend.position = 'none')


# 7.GO enrichment -----------------------------------------------------------------------------------
go_analysis_high_rna <- enrichGO(gene = subset(plot,significance=="Lower in protein")[,1],
                                 OrgDb = org.Mm.eg.db,
                                 keyType = "SYMBOL",
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)
go.rna <- simplify(go_analysis_high_rna,cutoff=0.7,by="p.adjust",select_fun=min)  

go_analysis_high_protein <- enrichGO(gene = subset(plot,significance=="Higher in protein")[,1],
                                     OrgDb = org.Mm.eg.db,
                                     keyType = "SYMBOL",
                                     ont = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05,
                                     readable = TRUE)

go.protein <- simplify(go_analysis_high_protein,cutoff=0.7,by="p.adjust",select_fun=min) 

go_1 <- go.protein@result %>%
  dplyr::select(c("Description","Count","p.adjust")) %>%
  mutate(logp.adjust = -log10(p.adjust)) %>%
  arrange(desc(logp.adjust)) %>%
  slice_head(n=10) %>%
  mutate(class="high protein abundance") %>%
  mutate(type=1)

go_2 <- go.rna@result %>%
  dplyr::select(c("Description","Count","p.adjust")) %>%
  mutate(logp.adjust = -log10(p.adjust)) %>%
  arrange(desc(logp.adjust)) %>%
  slice_head(n=10) %>%
  mutate(class="low protein abundance") %>%
  mutate(type=-1)

plot.data <- rbind(go_1,go_2) %>%
  mutate(x = -log10(p.adjust) * type,
         margin = ifelse(type == 1, -0.2, 0.2),
         hjust = ifelse(type == 1, 1, 0)) %>%
  arrange(x) %>%
  mutate(Description = factor(Description, levels = Description))

p_go <- ggplot(plot.data,aes(x, Description, colour = class)) +
  geom_segment(aes(xend = 0, yend = Description), linewidth = 1,
               show.legend = FALSE) +
  geom_point(aes(size = Count)) +
  geom_text(aes(x = margin, label = Description, hjust = hjust),
            show.legend = FALSE, colour = "black") +
  geom_hline(yintercept = 0) +
  labs(x = '-log10(p.adjust)', y = NULL) +
  scale_color_manual(values = c("high protein abundance" = "#d35230",
                                "low protein abundance" = "#2b7cd3")) +
  theme_bw() +
  theme(axis.ticks.x.top = element_blank(),  
        axis.title.x.top = element_blank(),  
        axis.text.x.top = element_blank(),  
        axis.line.x.top = element_blank(),
        axis.line.x = element_line(color = "black"), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 20, l = 20),
        panel.border = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    strip.placement = "outside",
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.4,"cm")
  )


# 8.Correlation Analysis of different cell-type scRNA-seq-----------------------------------------------------------------------------------
mm_rna_data <- readRDS("data/scRNA-seq/Tabula_Muris_Senis/mm.rds")
table(mm_rna_data@meta.data$age,mm_rna_data@meta.data$sex)

mm_rna1 <- subset(mm_rna_data, age %in% c("3m") & sex == "male") 
table(mm_rna1@meta.data$age,mm_rna1@meta.data$sex)
table(mm_rna1@meta.data$age,mm_rna1@meta.data$cell_type)
selected_cell_type <- c("fibroblast",
                        "mesenchymal stem cell", 
                        "chondrocyte",
                        "macrophage",
                        "B cell",
                        "keratinocyte",
                        "basal cell of epidermis",
                        "basal epithelial cell of tracheobronchial tree",
                        "fibroblast of lung",
                        "mesenchymal cell",
                        "bladder cell",
                        "bladder urothelial cell") 
mm_rna <- subset(mm_rna1, cell_type %in% selected_cell_type) 

scp_seu_2m <- scp_seu %>%
  subset(age == "2M")
rna_3m <- rna1_seu %>%
  subset(Age == "P100" & Sex == "Male") 

# Correlation calculation 
df <- as.data.frame(GetAssayData(scp_seu_2m, layer = "data"))  
average_expr_scp_2m <- data.frame(Protein = rownames(df),   
                                  MeanExpression = rowMeans(df))
average_expr_rna_3m <- data.frame(Gene = rownames(GetAssayData(rna_3m, layer = "data")),   
                                  MeanExpression = rowMeans(GetAssayData(rna_3m, layer = "data")))  
common_genes <- intersect(rownames(average_expr_scp_2m), rownames(average_expr_rna_3m))  
average_expr1_filtered <- average_expr_scp_2m[common_genes, , drop = FALSE]  
average_expr2_filtered <- average_expr_rna_3m[common_genes, , drop = FALSE]  
merged_cor <- left_join(average_expr1_filtered, 
                        average_expr2_filtered,
                        by = c("Protein" = "Gene"),
                        suffix = c(".protein", ".rna"))
correlation_3m <- cor(merged_cor$MeanExpression.protein, 
                      merged_cor$MeanExpression.rna, 
                      method = "spearman")  

correlation_results <- numeric(length(selected_cell_type))  
for (t in 1:length(selected_cell_type)){
  rna_t <- subset(mm_rna, cell_type == selected_cell_type[t])
  expr_rna <- GetAssayData(rna_t, layer = "data")
  average_expr_rna <- data.frame(Gene = rownames(expr_rna),   
                                 MeanExpression = rowMeans(expr_rna)) 
  name <- bitr(average_expr_rna$Gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Mm.eg.db")
  average_expr_rna <- inner_join(name,average_expr_rna,by=c("ENSEMBL"="Gene")) 
  
  common_genes <- intersect(rownames(average_expr_scp), average_expr_rna$SYMBOL)  
  average_expr1_filtered <- average_expr_scp[common_genes, , drop = FALSE]  
  
  merged_cor <- left_join(average_expr1_filtered, 
                          average_expr_rna[,2:3],
                          by = c("Protein" = "SYMBOL"),
                          suffix = c(".protein", ".rna"))
  
  correlation_results[t] <- cor(merged_cor$MeanExpression.protein,
                                merged_cor$MeanExpression.rna, 
                                method = "spearman")  
}
correlation_df <- data.frame(cell_type = c("Microglia",selected_cell_type), 
                             correlation = c(correlation_3m,correlation_results))%>%  
  mutate(cell_type = recode(cell_type,  
                            "basal epithelial cell of tracheobronchial tree" = "basal epithelial cell")) %>%
  mutate(cell_type = factor(cell_type,levels = c("Microglia",
                                                 "fibroblast",
                                                 "mesenchymal stem cell", 
                                                 "chondrocyte",
                                                 "macrophage",
                                                 "B cell",
                                                 "keratinocyte",
                                                 "basal cell of epidermis",
                                                 "basal epithelial cell",
                                                 "fibroblast of lung",
                                                 "mesenchymal cell",
                                                 "bladder cell",
                                                 "bladder urothelial cell")))
p_cor <- ggplot(correlation_df,aes(x=reorder(cell_type,-correlation),
                                   y=correlation,
                                   fill=cell_type))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(correlation,2)),vjust = -0.25,size=5) +
  scale_fill_manual(values = c("#fc4e07",rep("#595959", 12)))+ 
  theme_cowplot(font_size = 25)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.1,'cm')) +
  guides(fill = "none") +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5)) +
  xlab("Cell Type")







