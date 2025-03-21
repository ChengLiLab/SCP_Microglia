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
library(SCP)
library(harmony)
library(dplyr)
setwd('/Users/guangxinzhang/Documents/new_analyze')
source("code/scp_utils.R")

seurat_obj <- readRDS("seurat_object/SCP_MG_3184_seurat_zhr_step1_harmony_20250217_label_transferred_morpho_score_v5.rds")

####################### Fig S9C #######################

tmp.seu <- seurat_obj

tmp.seu[[]] <- tmp.seu[[]] %>% 
  mutate(age_region = paste0(age," ",brain_region))
tmp.seu[[]]$age_region <- factor(tmp.seu[[]]$age_region,
                                 levels = 
                                   # c("2M HIP","14M HIP","24M HIP",
                                   #          "2M PFC","14M PFC","24M PFC")
                                 c("2M HIP","2M PFC",
                                   "14M HIP","14M PFC",
                                   "24M HIP","24M PFC"))

ttt <- tmp.seu[,!tmp.seu$celltype == "BAM"]
DefaultAssay(ttt) <- "SCP"

current_assay_class <- class(ttt[["SCP"]])
print(paste("当前 Assay 的类为：", current_assay_class))

if (current_assay_class == "Assay5") {
  ttt[["SCP_v3"]] <- as(object = ttt[["SCP"]], Class = "Assay")
  DefaultAssay(ttt) <- "SCP_v3"
}

p <- FeatureStatPlot(ttt, stat.by = c("Diameter", "Elongation", "Intensity",
                                 "Circularity", "nFeature_SCP", "nCount_SCP"),
                group.by = "age_region", add_box = TRUE, plot_type = "box",
                ylab = "Value", xlab = "", bg.by = "age", ncol = 2,
                bg_palcolor = c("#1f78b4", "#33a02c", "#ff7f00"),
                comparisons = list(c("2M HIP", "2M PFC"),
                                   c("14M HIP", "14M PFC"),
                                   c("24M HIP", "24M PFC"),
                                   c("2M HIP", "24M HIP")))

p

####################### Fig S8C #######################

rna <- qread("seurat_object/Hammond_et-al-2019_Seurat_Converted_v4.qs")
rna <- subset(rna, subset = Sex == "Male")
rna <- subset(rna, subset = Condition == "WT")
rna <- subset(rna, subset = Age %in% c("P100", "P30","P540") & Batch == "8")

marker_lst <- read_excel("data/MG_marker_list.xlsx")
synapse_markers <- marker_lst[, c("Pre-synaptic marker", 
                                 "Post-synaptic marker",
                                 "GOCC_GABA_synapse(GO:0098982)",
                                 "GOCC_GLUTA_synapse(GO:0098978)", 
                                 "GOCC_Presynapse(GO:0098793)",
                                 "GOCC_Postsynapse(GO:0098794)",
                                 "GOCC_ExcitatorySynapse(GO:0060076)",
                                 "GOCC_InhibitorySynapse(GO:0060077)")]

synapse_markers_filled <- data.frame(
  lapply(synapse_markers, function(col_genes) {
    common_genes <- intersect(col_genes, rownames(seurat_obj))
    if (length(common_genes) < length(col_genes)) {
      c(common_genes, rep(NA, length(col_genes) - length(common_genes)))
    } else {
      common_genes
    }
  }),
  stringsAsFactors = FALSE
)

colnames(synapse_markers_filled) <- colnames(synapse_markers)
all_synapse_markers <- unname(unlist(synapse_markers_filled))
all_synapse_markers <- all_synapse_markers[!is.na(all_synapse_markers)]

count_matrix <- rna@assays$RNA@counts

synapse_expression <- sapply(all_synapse_markers, function(gene) {
  if(gene %in% rownames(count_matrix)) {
    cells_expressing <- sum(count_matrix[gene,] > 0)
    percentage <- (cells_expressing / ncol(count_matrix)) * 100
    return(percentage)
  } else {
    return(NA)
  }
})

names(synapse_expression) <- all_synapse_markers


data_matrix <- rna@assays$RNA@data

gene_means <- rowMeans(data_matrix)

gene_percentiles <- rank(gene_means) / length(gene_means) * 100

synapse_percentiles <- sapply(all_synapse_markers, function(gene) {
  if(gene %in% rownames(count_matrix)) {
    return(gene_percentiles[gene])
  } else {
    return(NA)
  }
})

names(synapse_percentiles) <- all_synapse_markers

synapse_stats <- data.frame(
  Gene = all_synapse_markers,
  Expression_Percentage = synapse_expression,
  Expression_Percentile = synapse_percentiles,
  row.names = NULL 
)

synapse_stats <- synapse_stats[!duplicated(synapse_stats$Gene), ]
synapse_stats <- synapse_stats[!is.na(synapse_stats$Expression_Percentile), ]

synapse_stats <- synapse_stats[order(synapse_stats$Expression_Percentile), ]

low_synapse_stats <- synapse_stats[synapse_stats$Expression_Percentile < 65, ]

low_synapse_stats <- low_synapse_stats[order(low_synapse_stats$Expression_Percentile), ]

synapse_genes <- low_synapse_stats$Gene

genes_to_keep <- setdiff(rownames(seurat_obj), synapse_genes)
genes_to_keep <- as.vector(genes_to_keep)

seurat_obj.tmp <- CreateSeuratObject(counts = seurat_obj[["SCP"]])
seurat_obj.tmp <- AddMetaData(seurat_obj.tmp, metadata = seurat_obj$celltype, col.name = "celltype")
seurat_obj.tmp@reductions$harmonyUMAP <- seurat_obj@reductions$harmonyUMAP
seurat_obj.tmp$age <- seurat_obj$age
seurat_obj.tmp$chromatogram_infor <- seurat_obj$chromatogram_infor
seurat_obj.tmp$prepare_date <- seurat_obj$prepare_date
seurat_obj.tmp$nFeature_SCP <- seurat_obj$nFeature_SCP

seurat_obj.tmp <- subset(seurat_obj.tmp, features = genes_to_keep)

seurat_obj.tmp <- SCTransform(seurat_obj.tmp,assay = "RNA",
                       vars.to.regress = "nFeature_SCP")

seurat_obj.tmp <- RunPCA(seurat_obj.tmp,npcs = 20)

pca_embeddings <- Embeddings(seurat_obj.tmp, reduction = "pca")

pca_embeddings <- pca_embeddings[, c(3:20)]

seurat_obj.tmp[["pca"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PC_", assay = DefaultAssay(seurat_obj.tmp))

seurat_obj.tmp <- RunHarmony(seurat_obj.tmp, group.by.vars = "prepare_date",
                      reduction.use = "pca", dims.use = 1:18)

seurat_obj.tmp <- RunUMAP(seurat_obj.tmp,dims = c(1:6),assay = "SCT",
                   reduction = "harmony",
                   reduction.name = "harmonyUMAP")

seurat_obj.tmp <- FindNeighbors(seurat_obj.tmp,
                     reduction = "harmony",
                     dims = 1:5)

seq <- seq(0.1,1,by = 0.1)
tmp.seu <- seurat_obj.tmp
for(res in seq){
  tmp.seu <- FindClusters(tmp.seu, resolution = res)
}

tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 %in% c("1", "9")] <- "1"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 %in% c("4", "7")] <- "4"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 %in% c("0","5", "6","8")] <- "0"


tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "10"] <- "7"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "4"]  <- "5"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "11"] <- "4"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "2"]  <- "6"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "1"]  <- "2"
tmp.seu$SCT_snn_res.0.7[tmp.seu$SCT_snn_res.0.7 == "0"]  <- "1"


df <- data.frame(
    celltype = tmp.seu@meta.data$celltype,
    celltype_without_phago_genes = tmp.seu@meta.data$`SCT_snn_res.0.7`
)

df <- df %>% group_by(celltype, celltype_without_phago_genes) %>% summarise('Phagocytosis Content' = n())
df$'Cell Type' <- as.character(df$celltype)

df <-ungroup(df)


p1 <- SankeyPlot(df, x = c(".", "Phagocytosis Content"),
    stratum = "celltype_without_phago_genes",  links_fill_by = "Cell Type", xlab = NULL)

p1