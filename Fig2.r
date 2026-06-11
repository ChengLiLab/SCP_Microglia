# Fig2.R
# Purpose: Identify Potential Exogenous Proteins

##----0. load package--------
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(cowplot)
library(Seurat)
library(qs)
source("scp_utils.R")

# -------- **Filter exogenous genes --------------------------------------------------------
# -------- 1. Filter by expression level ------------------------------------------
scp.data <- readRDS('res/SCP_MG_3085_seurat_post_QC.rds')

scp.data <- NormalizeData(scp.data)
scp.data[[]]$group <- paste0(scp.data[[]]$age," ",scp.data[[]]$brain_region)
scp.data[[]]$group <- factor(scp.data[[]]$group,levels = c("2M PFC","14M PFC","24M PFC",
                                                           "2M HIP","14M HIP","24M HIP"))
seu_obj <- scp.data
head(seu_obj[[]])


matrix_list <- readRDS("res/matrix_list_4_datasets.rds")
hammond_all_samples <- qread("Hammond_et-al-2019_Seurat_Converted_v4.qs")
immunity2019 <- subset(hammond_all_samples,
                       features = intersect(rownames(matrix_list[["2019 Immunity"]]), rownames(hammond_all_samples)),
                       cells = intersect(colnames(matrix_list[["2019 Immunity"]]), colnames(hammond_all_samples))) %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData()
immunity2019.rna <- subset(immunity2019, 
                           Age == "P100" & Sex == "Male")

neuron2019 <- read_rds("neuron2019_scRNA.rds")
neuron2019.rna <- subset(neuron2019,
                         features = intersect(rownames(matrix_list[["2019 Neuron"]]), rownames(neuron2019)),
                         cells = intersect(colnames(matrix_list[["2019 Neuron"]]), colnames(neuron2019))) %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 20) %>% 
  RunUMAP(dims = 1:20,
          reduction = "pca",
          reduction.name = "umap")



neuron2019.rna.cortex <- subset(neuron2019.rna, 
                                tissue == "Cortex" & age == "P60")
neuron2019.rna.hippocampus <- subset(neuron2019.rna, 
                                     tissue == "Hippocampus" & age == "P60")


scp.genes <- as.vector(rownames(scp.data))

# Utility function to calculate gene expression percentage
calc_gene_pct <- function(gene_list,
                          count_mat,
                          threshold = 0,
                          as_df = TRUE) {
  # Parameters description
  # gene_list  : Character vector, gene names to check
  # count_mat  : Gene x Cell count matrix (matrix / dgCMatrix)
  # threshold  : Lower limit to determine "expression", default > 0
  # as_df      : TRUE returns a dataframe, FALSE returns a named vector
  
  # Convert genes and matrix rownames to character to avoid factor errors
  genes <- as.character(gene_list)
  rn    <- rownames(count_mat)
  
  # sapply preserves order, returns NA for missing genes
  pct_vec <- sapply(genes, function(g) {
    if (g %in% rn) {
      n_pos <- sum(count_mat[g, ] > threshold)
      100 * n_pos / ncol(count_mat)
    } else {
      NA_real_
    }
  }, USE.NAMES = TRUE)
  
  if (as_df) {
    out <- data.frame(
      gene       = names(pct_vec),
      pct_cells  = pct_vec,
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  } else {
    return(pct_vec)
  }
}

# -------- Calculate expression percentage of SCP genes in scRNA data --------

neuron2019.rna.express.pct.cortex <- calc_gene_pct(scp.genes, 
                                                   neuron2019.rna.cortex@assays$RNA$data, 
                                                   threshold = 0) %>%
  as.data.frame(.)

neuron2019.rna.express.pct.hippocampus <- calc_gene_pct(scp.genes, 
                                                        neuron2019.rna.hippocampus@assays$RNA$data, 
                                                        threshold = 0) %>%
  as.data.frame(.)

immunity2019.rna.express.pct <- calc_gene_pct(scp.genes, 
                                              immunity2019.rna@assays$RNA$data, 
                                              threshold = 0) %>%
  as.data.frame(.)


# General function to calculate the percentile of average gene expression
calc_gene_percentile <- function(gene_list,
                                 expr_mat,
                                 as_df = TRUE,
                                 ties = "average") {
  gene_means <- rowMeans(expr_mat)
  
  gene_percentiles <- rank(gene_means,
                           ties.method = ties) /
    length(gene_means) * 100
  pct_vec <- sapply(gene_list, function(g) {
    if (g %in% names(gene_percentiles)) {
      gene_percentiles[g]
    } else {
      NA_real_
    }
  }, USE.NAMES = FALSE)
  
  if (as_df) {
    out <- data.frame(
      gene        = names(pct_vec),
      percentile  = pct_vec,
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  } else {
    return(pct_vec)
  }
}



# -------- Calculate average expression percentile of SCP genes in scRNA data --------

tmp.neuron2019.rna.express.pct.cortex <- neuron2019.rna.express.pct.cortex %>% filter(pct_cells > 0)
tmp.neuron2019.rna.pct.cortex.matrix <- neuron2019.rna.cortex@assays$RNA$data[tmp.neuron2019.rna.express.pct.cortex$gene, ]
neuron2019.rna.pct.cortex.rank <- calc_gene_percentile(scp.genes, 
                                                       tmp.neuron2019.rna.pct.cortex.matrix, 
                                                       as_df = TRUE)


tmp.neuron2019.rna.express.pct.hippocampus <- neuron2019.rna.express.pct.hippocampus %>% filter(pct_cells > 0)
tmp.neuron2019.rna.pct.hippocampus.matrix <- neuron2019.rna.hippocampus@assays$RNA$data[tmp.neuron2019.rna.express.pct.hippocampus$gene, ]
neuron2019.rna.pct.hippocampus.rank <- calc_gene_percentile(scp.genes, 
                                                            tmp.neuron2019.rna.pct.hippocampus.matrix, 
                                                            as_df = TRUE)

tmp.immunity2019.rna.express.pct <- immunity2019.rna.express.pct %>% filter(pct_cells > 0)
tmp.immunity2019.rna.pct.matrix <- immunity2019.rna@assays$RNA$data[tmp.immunity2019.rna.express.pct$gene, ]
immunity2019.rna.pct.rank <- calc_gene_percentile(scp.genes, 
                                                  tmp.immunity2019.rna.pct.matrix, 
                                                  as_df = TRUE)


# -------- Merge results --------
neuron2019.rna.express.pct.cortex[is.na(neuron2019.rna.express.pct.cortex)] <- 0
neuron2019.rna.pct.cortex.rank[is.na(neuron2019.rna.pct.cortex.rank)] <- 0
neuron2019.rna.express.pct.hippocampus[is.na(neuron2019.rna.express.pct.hippocampus)] <- 0
neuron2019.rna.pct.hippocampus.rank[is.na(neuron2019.rna.pct.hippocampus.rank)] <- 0
immunity2019.rna.express.pct[is.na(immunity2019.rna.express.pct)] <- 0
immunity2019.rna.pct.rank[is.na(immunity2019.rna.pct.rank)] <- 0

immunity2019.rna.summary <- immunity2019.rna.express.pct %>%
  left_join(immunity2019.rna.pct.rank,
            by = "gene")
neuron2019.rna.hippocampus.summary <- neuron2019.rna.express.pct.hippocampus %>%
  left_join(neuron2019.rna.pct.hippocampus.rank,
            by = "gene")
neuron2019.rna.cortex.summary <- neuron2019.rna.express.pct.cortex %>%
  left_join(neuron2019.rna.pct.cortex.rank,
            by = "gene")   

neuron2019.rna.cortex.summary[is.na(neuron2019.rna.cortex.summary)] <- 0
neuron2019.rna.hippocampus.summary[is.na(neuron2019.rna.hippocampus.summary)] <- 0
immunity2019.rna.summary[is.na(immunity2019.rna.summary)] <- 0

neuron2019.rna.cortex.summary <- neuron2019.rna.cortex.summary %>%
  filter(., pct_cells < 5 & percentile < 15)
neuron2019.rna.hippocampus.summary <- neuron2019.rna.hippocampus.summary %>%
  filter(., pct_cells < 5 & percentile < 15)
immunity2019.rna.summary <- immunity2019.rna.summary %>%
  filter(., pct_cells < 5 & percentile < 15)

ex.genes.cortex <- neuron2019.rna.cortex.summary$gene
ex.genes.hippocampus <- neuron2019.rna.hippocampus.summary$gene
ex.genes.immunity2019 <- immunity2019.rna.summary$gene



scp.data.cortex <- subset(scp.data, brain_region == "PFC")

scp.data.cortex.matrix <- scp.data.cortex@assays$SCP$counts
gene2keep <- rownames(scp.data.cortex.matrix[rowSums(scp.data.cortex.matrix) > 0, ])
scp.data.cortex <- subset(
  scp.data.cortex,
  features = gene2keep
)

scp.data.hippocampus <- subset(scp.data, brain_region == "HIP")
scp.data.hippocampus.matrix <- scp.data.hippocampus@assays$SCP$counts
gene2keep <- rownames(scp.data.hippocampus.matrix[rowSums(scp.data.hippocampus.matrix) > 0, ])

scp.data.hippocampus <- subset(
  scp.data.hippocampus,
  features = gene2keep
)

scp.data.cortex <- NormalizeData(scp.data.cortex)

scp.data.hippocampus <- NormalizeData(scp.data.hippocampus)

scp.express.pct.hippocampus.rank <- calc_gene_percentile(scp.genes, 
                                                         scp.data.hippocampus@assays$SCP$data, 
) %>%
  as.data.frame(.)
scp.express.pct.cortex.rank <- calc_gene_percentile(scp.genes, 
                                                    scp.data.cortex@assays$SCP$data, 
) %>%
  as.data.frame(.)

scp.data <- NormalizeData(scp.data)
scp.express.pct.rank <- calc_gene_percentile(scp.genes, 
                                             scp.data@assays$SCP$data, 
) %>%
  as.data.frame(.)

scp.genes.cortex <- scp.express.pct.cortex.rank %>%
  filter(percentile < 70) %>%
  dplyr::select(gene) 

scp.genes.hippocampus <- scp.express.pct.hippocampus.rank %>%
  filter(percentile < 70) %>%
  dplyr::select(gene)

scp.genes.all <- scp.express.pct.rank %>%
  filter(percentile < 70) %>%
  dplyr::select(gene)

scp.genes.exclu <- scp.express.pct.rank %>%
  filter(percentile > 70) %>%
  dplyr::select(gene)

norm_data <- GetAssayData(scp.data, slot = "data")
# Calculate the average expression of each gene across all cells
gene_avg_abundance <- Matrix::rowMeans(norm_data)



## GO Enrichment Analysis of excluded proteins based on mRNA expression
enrich.go <- enrichGO(
  gene = scp.genes.exclu$gene,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE)
enrich.go_simplify <- simplify(enrich.go,cutoff=0.7,by="p.adjust",select_fun=min) 

p1 <- barplot(enrich.go_simplify, showCategory = 30) + 
  ggtitle("Exclude genes GO(CC) Enrichment Analysis (n=1350)") +
  theme(plot.title = element_text(hjust = 0.5))
p1


cortex.ex.genes <- intersect(scp.genes.cortex$gene, ex.genes.cortex)
hippocampus.ex.genes <- intersect(scp.genes.hippocampus$gene, ex.genes.hippocampus)
immunity2019.ex.genes <- intersect(scp.genes.all$gene, ex.genes.immunity2019)


ex.genes_1 <- intersect(immunity2019.ex.genes,unique(union(cortex.ex.genes,hippocampus.ex.genes)))

# -------- 3. Compare with Allen Institute whole brain data -----------------------------------
# sceasy::convertFormat(obj = "public data/nature2025_mouse_adult_male_26w.h5ad",
#                       from = "anndata",
#                       to = "seurat",
#                       outFile = "public data/nature2025_mouse_adult_male_26w.rds")
#rna.brain <- ref
rna.brain <- readRDS("public data/nature2025_mouse_adult_male_26w.rds")
rna.brain <- rna.brain %>% NormalizeData()

SCP::CellDimPlot(rna.brain,pt.size = 1,
                 group.by = "class_label",reduction = "umap")
table(rna.brain[[]]$class_label)

rna.brain@meta.data$celltype <- rna.brain@meta.data$class_label
metadata <- rna.brain[[]]
metadata$celltype <- as.character(metadata$celltype)
metadata[metadata$subclass_label == "Astro-NT NN",][["celltype"]] = 'Astrocyte'
metadata[metadata$subclass_label == "Astro-TE NN",][["celltype"]] = 'Astrocyte'
metadata[metadata$subclass_label == "Ependymal NN",][["celltype"]] = 'Ependymal'
metadata[metadata$subclass_label == "Hypendymal NN",][["celltype"]] = 'Hypendymal'
metadata[metadata$subclass_label == "Tanycyte NN",][["celltype"]] = 'Tanycyte'
metadata[metadata$subclass_label == "Astroependymal NN",][["celltype"]] = 'Astroependymal'
metadata[metadata$class_label == "CNU-HYa GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "CNU-HYa Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "CNU-LGE GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "CNU-MGE GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "CTX-CGE GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "CTX-MGE GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "DG-IMN Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "HY GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "HY Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "HY MM Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "IT-ET Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$subclass_label == "BAM NN",][["celltype"]] = 'BAM'
metadata[metadata$subclass_label == "Microglia NN",][["celltype"]] = 'Microglia'
metadata[metadata$subclass_label == "DC NN",][["celltype"]] = 'DC cells'
metadata[metadata$subclass_label == "Lymphoid NN",][["celltype"]] = 'Lymphoid'
metadata[metadata$class_label == "MB Dopa",][["celltype"]] = 'Dopa Neurons'
metadata[metadata$class_label == "MB GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "MB Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "MB-HB Sero",][["celltype"]] = 'Sero Neurons'
metadata[metadata$class_label == "MY Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "NP-CT-L6b Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "OB-CR Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "OB-IMN GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$subclass_label == "OPC NN",][["celltype"]] = 'OPC'
metadata[metadata$subclass_label == "Oligo NN",][["celltype"]] = 'Oligocyte'
metadata[metadata$class_label == "P GABA",][["celltype"]] = 'GABA Neurons'
metadata[metadata$class_label == "P Glut",][["celltype"]] = 'Glut Neurons'
metadata[metadata$class_label == "Vascular",][["celltype"]] = 'Vascular'
rna.brain@meta.data <- metadata
Idents(rna.brain) <- rna.brain[[]]$celltype

rna.brain[[]]$celltype <- factor(rna.brain[[]]$celltype,
                                 levels = c("Glut Neurons",
                                            "Vascular", 
                                            "Oligocyte",
                                            "Microglia",
                                            "GABA Neurons",
                                            "Dopa Neurons",
                                            "Sero Neurons",
                                            "Astrocyte",
                                            "OPC",
                                            "BAM",
                                            "Lymphoid",
                                            "Ependymal",
                                            "DC cells",
                                            "Tanycyte",
                                            "Hypendymal",
                                            "Astroependymal"))

rna.brain.main <- rna.brain %>% 
  subset(celltype %in% c("Glut Neurons","Vascular","Oligodendrocyte",
                         "GABA Neurons","Astrocyte","OPC","Microglia"))

SCP::CellDimPlot(rna.brain.main,pt.size = 1,theme_use = "theme_blank",
                 group.by = "celltype",reduction = "umap",
                 palcolor =  c("Glut Neurons" = "#fda542",
                               "Vascular" = "#a6cee3",
                               "Oligocyte" = "#63a8a0",
                               "Microglia" = "#98d277",
                               "GABA Neurons" = "#3ba432",
                               "Dopa Neurons" = "#7d54a5",
                               "Sero Neurons" = "#ead27a",
                               "Astrocyte" = "#438ec0",
                               "OPC"="#b9b458",
                               "BAM" = "#fe8214",
                               "Lymphoid"="#fb9684",
                               "Ependymal"="#ec4d4e",
                               "DC cells"="#da4c59",
                               "Tanycyte"="#c3aad2",
                               "Hypendymal"="#b9a499",
                               "Astroependymal"="#b15928")
)

SCP::CellDimPlot(rna.brain,pt.size = 1,theme_use = "theme_blank",
                 group.by = "celltype",reduction = "umap",
                 palcolor =  c("Glut Neurons" = "#fda542",
                               "Vascular" = "#a6cee3",
                               "Oligocyte" = "#63a8a0",
                               "Microglia" = "#98d277",
                               "GABA Neurons" = "#3ba432",
                               "Dopa Neurons" = "#7d54a5",
                               "Sero Neurons" = "#ead27a",
                               "Astrocyte" = "#438ec0",
                               "OPC"="#b9b458",
                               "BAM" = "#fe8214",
                               "Lymphoid"="#fb9684",
                               "Ependymal"="#ec4d4e",
                               "DC cells"="#da4c59",
                               "Tanycyte"="#c3aad2",
                               "Hypendymal"="#b9a499",
                               "Astroependymal"="#b15928")
)

#DEG <- FindAllMarkers(rna.brain, only.pos = T,test.use = 'wilcox')
DEG.vas <- FindMarkers(rna.brain.main, ident.1 = "Vascular",
                       logfc.threshold = 0.2,min.pct = 0.1,
                       only.pos = T,test.use = 'wilcox')
DEG.ast <- FindMarkers(rna.brain.main, ident.1 = "Astrocyte",
                       logfc.threshold = 0.2,min.pct = 0.1,
                       only.pos = T,test.use = 'wilcox')
DEG.olig <- FindMarkers(rna.brain.main, ident.1 = "Oligocyte",
                        logfc.threshold = 0.2,min.pct = 0.1,
                        only.pos = T,test.use = 'wilcox')
DEG.mg <- FindMarkers(rna.brain.main, ident.1 = "Microglia",
                      logfc.threshold = 0.2,min.pct = 0.1,
                      only.pos = T,test.use = 'wilcox')
DEG.gaba <- FindMarkers(rna.brain.main, ident.1 = "GABA Neurons",
                        logfc.threshold = 0.2,min.pct = 0.1,
                        only.pos = T,test.use = 'wilcox')
DEG.glut <- FindMarkers(rna.brain.main, ident.1 = "Glut Neurons",
                        logfc.threshold = 0.2,min.pct = 0.1,
                        only.pos = T,test.use = 'wilcox')
DEG.opc <- FindMarkers(rna.brain.main, ident.1 = "OPC",
                       logfc.threshold = 0.2,min.pct = 0.1,
                       only.pos = T,test.use = 'wilcox')
# DEG.bam <- FindMarkers(rna.brain, ident.1 = "BAM",
#                        only.pos = T,test.use = 'wilcox')
# DEG.lym <- FindMarkers(rna.brain, ident.1 = "Lymphoid",
#                        only.pos = T,test.use = 'wilcox')
# DEG.epen <- FindMarkers(rna.brain, ident.1 = "Ependymal",
#                        only.pos = T,test.use = 'wilcox')
# DEG.dc <- FindMarkers(rna.brain, ident.1 = "DC cells",
#                        only.pos = T,test.use = 'wilcox')
# DEG.tan <- FindMarkers(rna.brain, ident.1 = "Tanycyte",
#                        only.pos = T,test.use = 'wilcox')
# DEG.dopa <- FindMarkers(rna.brain, ident.1 = "Dopa Neurons",
#                        only.pos = T,test.use = 'wilcox')
# DEG.hypen <- FindMarkers(rna.brain, ident.1 = "Hypendymal",
#                        only.pos = T,test.use = 'wilcox')
# DEG.sero <- FindMarkers(rna.brain, ident.1 = "Sero Neurons",
#                        only.pos = T,test.use = 'wilcox')
# DEG.astepen <- FindMarkers(rna.brain, ident.1 = "Astroependymal",
#                        only.pos = T,test.use = 'wilcox')


DEG_list <- list(
  DEG.vas = cbind(cluster = "Vascular", DEG.vas) %>% rownames_to_column("gene"),
  DEG.ast = cbind(cluster = "Astrocyte", DEG.ast)%>% rownames_to_column("gene"),
  DEG.olig = cbind(cluster = "Oligocyte", DEG.olig)%>% rownames_to_column("gene"),
  DEG.mg = cbind(cluster = "Microglia", DEG.mg)%>% rownames_to_column("gene"),
  DEG.gaba = cbind(cluster = "GABA Neurons", DEG.gaba)%>% rownames_to_column("gene"),
  DEG.glut = cbind(cluster = "Glut Neurons", DEG.glut)%>% rownames_to_column("gene"),
  DEG.opc = cbind(cluster = "OPC", DEG.opc)%>% rownames_to_column("gene")
  # DEG.bam = cbind(cluster = "BAM", DEG.bam)%>% rownames_to_column("gene"),
  # DEG.lym = cbind(cluster = "Lymphoid", DEG.lym)%>% rownames_to_column("gene"),
  # DEG.epen = cbind(cluster = "Ependymal", DEG.epen)%>% rownames_to_column("gene"),
  # DEG.dc = cbind(cluster = "DC cells", DEG.dc)%>% rownames_to_column("gene"),
  # DEG.tan = cbind(cluster = "Tanycyte", DEG.tan)%>% rownames_to_column("gene"),
  # DEG.dopa = cbind(cluster = "Dopa Neurons", DEG.dopa)%>% rownames_to_column("gene"),
  # DEG.hypen = cbind(cluster = "Hypendymal", DEG.hypen)%>% rownames_to_column("gene"),
  # DEG.sero = cbind(cluster = "Sero Neurons", DEG.sero)%>% rownames_to_column("gene"),
  # DEG.astepen = cbind(cluster = "Astroependymal", DEG.astepen)%>% rownames_to_column("gene")
)

# Merge by rows
combined.DEG <- do.call(rbind, DEG_list)
DEG_rna <- combined.DEG %>% mutate(pct.diff = pct.1-pct.2)
write.csv(DEG_rna,"25w_brain_rna_main_CellType_DEG_1021.csv")
DEG_rna <- read.csv("25w_brain_rna_main_CellType_DEG_1021.csv")

MG_enrich <- DEG_rna %>% 
  filter(p_val_adj<0.01 & avg_log2FC > 0.5 & cluster == "Microglia")

ex.genes_2_rna <- DEG_rna %>% filter(gene %in% ex.genes_1) %>% 
  filter(p_val_adj<0.01 & avg_log2FC > 0.5 & pct.diff >0.2) %>% 
  filter(!(gene %in% MG_enrich$gene ))
length(unique(ex.genes_2_rna$gene))
ex.genes_2 <- unique(ex.genes_2_rna$gene)
write.csv(ex.genes_2,"data/exo_genes_list.csv")

# Figure 2B
ex.genes_2_rna <- DEG_rna %>% filter(gene %in% ex.genes_2) %>% 
  filter(p_val_adj<0.01 & avg_log2FC > 0.5 & pct.diff >0.2)  
top5 <- ex.genes_2_rna %>% 
  group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice_head(n=5)
ggplot(ex.genes_2_rna, aes(x = factor(cluster), y = pct.diff)) +
  geom_point(aes(color = factor(cluster)), 
             #size = avg_log2FC,
             alpha = 1, position = position_jitter(width = 0.2)) +
  ggrepel::geom_text_repel(data = top5, force = 0.2,
                           max.overlaps=20,
                           aes(label = gene), size = 5, 
                           box.padding = unit(0.1, "lines"), 
                           point.padding = unit(0.1, "lines"), 
                           segment.size = 0.5)+
  scale_color_manual(values = c("Glut Neurons" = "#fda542",
                                "Vascular" = "#a6cee3",
                                "Oligocyte" = "#63a8a0",
                                "Microglia" = "#98d277",
                                "GABA Neurons" = "#3ba432",
                                "Dopa Neurons" = "#7d54a5",
                                "Sero Neurons" = "#ead27a",
                                "Astrocyte" = "#438ec0",
                                "OPC"="#b9b458",
                                "BAM" = "#fe8214",
                                "Lymphoid"="#fb9684",
                                "Ependymal"="#ec4d4e",
                                "DC cells"="#da4c59",
                                "Tanycyte"="#c3aad2",
                                "Hypendymal"="#b9a499",
                                "Astroependymal"="#b15928")) + 
  scale_y_continuous(limits = c(0.2, 1.0)) +
  labs(
    title = "Phagoproteome enriched in brain cells",
    x = "CellType",
    y = "Percentage difference in scRNA-seq",
    color = "CellType"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 13,colour = "black"),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle =45,vjust = 1,hjust = 1),
    legend.position = "none")

#Figure 2C
p <- ggplot(cluster_counts, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("Glut Neurons" = "#fda542",
                               "Vascular" = "#a6cee3",
                               "Oligocyte" = "#63a8a0",
                               "Microglia" = "#98d277",
                               "GABA Neurons" = "#3ba432",
                               "Dopa Neurons" = "#7d54a5",
                               "Sero Neurons" = "#ead27a",
                               "Astrocyte" = "#438ec0",
                               "OPC"="#b9b458",
                               "BAM" = "#fe8214",
                               "Lymphoid"="#fb9684",
                               "Ependymal"="#ec4d4e",
                               "DC cells"="#da4c59",
                               "Tanycyte"="#c3aad2",
                               "Hypendymal"="#b9a499",
                               "Astroependymal"="#b15928")) + 
  coord_polar("y", start = 0) +
  geom_text_repel(aes(label = label), 
                  position = position_stack(vjust = 0.5),force = 1,
                  size = 4, color = "white", fontface = "bold") +
  labs(title = "Distribution of Phagocyted proteins",
       fill = "CellType") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none")
p

#Figure 2D
enrich.go <- enrichGO(
  gene = ex.genes_2,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE)
enrich.go_simplify <- simplify(enrich.go,cutoff=0.7,by="p.adjust",select_fun=min) 

enrich.go_simplify_plot <- enrich.go_simplify %>% as.data.frame() %>% 
  mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>% 
  arrange(desc(`-log10(p.adjust)`)) %>% 
  slice_head(n=15)
enrich.go_simplify_plot$Description <- factor(enrich.go_simplify_plot$Description,
                                              levels = rev(enrich.go_simplify_plot$Description))

p_go <- ggplot(enrich.go_simplify_plot,aes(x = `-log10(p.adjust)`,
                                           y = Description))  +  
  geom_bar(stat="identity",width=0.8 ) + 
  ggthemes::theme_few() +
  #scale_fill_gradient2(low="#FFFFCCFF",mid="#FD8D3CFF",high ="#b04130",midpoint = 12.5) +  
  labs(x="-Log10(p.adjust)",y=" ",title=" ",fill="-log10 (p.adj)") + 
  ggtitle("GO(CC) Enrichment Analysis (n=365)") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=10,color="black"),
        axis.title = element_text(size=16),title = element_text(size=13)) 
p_go

#Figure 2E
p  <-  Seurat::DotPlot(object = scp.data, 
                       features = ex.genes_2, 
                       group.by = "mass") 

avg_expression <- p$data
avg_expression_df <- as.data.frame(avg_expression) 
Top50_genes <- avg_expression_df %>%
  arrange(desc(pct.exp)) %>%
  head(50)
Seurat::DotPlot(scp.data, 
                features = rev(rownames(Top50_genes)),
                group.by = "group") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggtitle("Top 50 engulfed proteins") +
  scale_color_gradientn(colors = sc.hic.orange(100)) +
  #coord_flip()+
  scale_x_discrete(labels = function(x)toupper(x))

#Figure 2F
p <- Seurat::DotPlot(seu_obj,
                     features = c(ex.genes_2),
                     group.by = "group") 
df <- p@data

df2 <- df %>% 
  group_by(id) %>% 
  summarise(
    avg_exp_mean = mean(avg.exp, na.rm = TRUE),
    pct_exp_mean = mean(pct.exp, na.rm = TRUE)
  ) %>% ungroup() %>% 
  mutate(scale_avg_exp_mean = scale(avg_exp_mean)) %>% 
  mutate(features.plot = "Phagoproteome")

seu_obj[[]]$percentage <- PercentageFeatureSet(seu_obj,
                                               features = ex.genes_2)
df2 <- seu_obj[[]] %>% 
  group_by(group) %>% 
  summarise(Percentage = mean(percentage, na.rm = TRUE)) %>% 
  mutate(Percentage = Percentage * 100)

df3 <- df %>% 
  group_by(id) %>% 
  summarise(
    avg_exp_mean = mean(avg.exp, na.rm = TRUE),
    pct_exp_mean = mean(pct.exp, na.rm = TRUE)
  ) %>% ungroup() %>% 
  mutate(scale_avg_exp_mean = scale(avg_exp_mean)) %>% 
  mutate(features.plot = "Phagoproteome") %>% 
  left_join(df2,by = c("id" = "group"))
ggplot(data = df3, mapping = aes_string(x = "features.plot", 
                                        y = "id")) + 
  geom_point(mapping = aes_string(size = "Percentage",
                                  color = "scale_avg_exp_mean[,1]")) + 
  #scale.func(range = c(0, 6), 
  #limits = c(scale.min, scale.max)) + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  guides(size = guide_legend(title = "Average\nPercent Detected")) + 
  labs(x = "Features") + theme_cowplot()+
  ggtitle("") +
  scale_color_gradientn(
    name = "z-score",
    colors = sc.hic.orange(100))+
  theme_bw()+
  #coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(
          hjust = 0.5,  # Horizontal alignment, 0.5 means center
          size = 16,    # Title font size
          face = "bold",  # Title bold
        ),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#Figure S6D
enrich.go <- enrichGO(
  gene = ex.genes_1,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE)
enrich.go_simplify <- simplify(enrich.go,cutoff=0.7,by="p.adjust",select_fun=min) 
enrich.go_simplify_plot <- enrich.go_simplify %>% as.data.frame() %>% 
  mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>% 
  arrange(desc(`-log10(p.adjust)`)) %>% 
  slice_head(n=15)
enrich.go_simplify_plot$Description <- factor(enrich.go_simplify_plot$Description,
                                              levels = rev(enrich.go_simplify_plot$Description))

p_go <- ggplot(enrich.go_simplify_plot,aes(x = `-log10(p.adjust)`,
                                           y = Description))  +  
  geom_bar(stat="identity",width=0.8 ) + 
  ggthemes::theme_few() +
  #scale_fill_gradient2(low="#FFFFCCFF",mid="#FD8D3CFF",high ="#b04130",midpoint = 12.5) +  
  labs(x="-Log10(p.adjust)",y=" ",title=" ",fill="-log10 (p.adj)") + 
  ggtitle("GO(CC) Enrichment Analysis (n=823)") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=10,color="black"),
        axis.title = element_text(size=16),title = element_text(size=13)) 
p_go

# Figure S6E
diff <- setdiff(unique(ex.genes_1), ex.genes_2)
enrich.go <- enrichGO(
  gene = diff,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE)
enrich.go_simplify <- simplify(enrich.go,cutoff=0.7,by="p.adjust",select_fun=min) 
enrich.go_simplify_plot <- enrich.go_simplify %>% as.data.frame() %>% 
  mutate(`-log10(p.adjust)` = -log10(p.adjust)) %>% 
  arrange(desc(`-log10(p.adjust)`)) %>% 
  slice_head(n=15)
enrich.go_simplify_plot$Description <- factor(enrich.go_simplify_plot$Description,
                                              levels = rev(enrich.go_simplify_plot$Description))

p_go <- ggplot(enrich.go_simplify_plot,aes(x = `-log10(p.adjust)`,
                                           y = Description))  +  
  geom_bar(stat="identity",width=0.8 ) + 
  ggthemes::theme_few() +
  #scale_fill_gradient2(low="#FFFFCCFF",mid="#FD8D3CFF",high ="#b04130",midpoint = 12.5) +  
  labs(x="-Log10(p.adjust)",y=" ",title=" ",fill="-log10 (p.adj)") + 
  ggtitle("Excluded proteins in Step 2 and Step 3 (n=458) \n GO(CC Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=10,color="black"),
        axis.title = element_text(size=16),title = element_text(size=13)) 
p_go



