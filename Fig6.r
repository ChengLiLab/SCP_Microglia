# Fig6.R
# Purpose: Analyze age- and brain-region-associated proteomic differences in microglia.

# 0.Load Packages -------------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(ggthemes)
library(patchwork)
library(ggrepel)
library(cowplot)
library(scplotter)
library(msigdbr)
library(UCell)

#----------------- 1.Identification of spatiotemporal subclasses of microglia -----------------------------------
seu <- readRDS("res/SCP_MG_3085_seurat_clustering.rds")
DefaultAssay(seu) <- "SCP"

seu$brain_region <- factor(seu$brain_region, levels = c("PFC", "HIP"))
seu$brain_region_label <- dplyr::recode(as.character(seu$brain_region), "PFC" = "Prefrontal cortex", "HIP" = "Hippocampus")
seu$brain_region_label <- factor(seu$brain_region_label, levels = c("Prefrontal cortex", "Hippocampus"))
age_groups <- list("2M", "14M", "24M")
brain_colors <- c("Prefrontal cortex" = "#9cca62", "Hippocampus" = "#00758b")
colnames(seu@reductions$bbknnUMAP@cell.embeddings) <- c("UMAP 1", "UMAP 2")

umap_plots <- lapply(age_groups, function(age) {
  seu_sub <- subset(seu, age == !!age)
  scplotter::CellDimPlot(seu_sub,
                         group_by = "brain_region_label", 
                         facet_by = "brain_region_label", 
                         reduction = "bbknnUMAP",
                         theme = "theme_blank", 
                         legend.position = "top",
                         highlight_stroke = 0,
                         highlight_size = 0.5,
                         palcolor = brain_colors) +
    labs(x = NULL, y = NULL, color = NULL, fill = NULL) +
    guides(color = guide_legend(nrow = 1,
                                byrow = TRUE,
                                override.aes = list(size = 2.5, alpha = 1)),
           fill = "none") +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          panel.spacing.x = unit(0.8, "cm"),
          panel.spacing.y = unit(0.35, "cm"),
          plot.margin = margin(t = 2, r = 3, b = 12, l = 10),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.25, "cm"),
          legend.spacing.x = unit(0.25, "cm"))
})
combined_plot <- patchwork::wrap_plots(umap_plots,
                                       ncol = 1,
                                       guides = "collect") &
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.25, "cm"))

# cell proportion
summarydata <- seu@meta.data %>% 
  group_by(brain_region, age, plot_cluster) %>% 
  dplyr::summarise(.,count=n())%>%
  group_by(brain_region, age) %>% 
  dplyr::mutate(Percentage = round(count/sum(count)*100,2)) %>%
  ungroup()
cluster_colors <- c("#EE8172", "#7DBFA7", "#AECDE1", "#EE934E", "#A5917F")
p <- ggplot(summarydata, aes(x = age, y = Percentage, fill = plot_cluster)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.9, width = 0.7) +
  facet_wrap(~ brain_region, ncol = 2) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "",
       x = "Age",
       y = "Proportion (%)",
       fill = "Cluster") +
  theme_cowplot(font_size = 15) +
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5))


# 2.Brain-region differential proteins within each age group --------------------------------------------------------------------------------------------
DefaultAssay(seu) <- "SCP"

compute_brain_dep_by_age <- function(seu_obj, 
                                     age_col = "age",
                                     brain_col = "brain_region",
                                     min_cells = 10,
                                     logfc_threshold = 0.5,
                                     pval_threshold = 0.05) {
  
  ages <- unique(seu_obj@meta.data[[age_col]])
  all_results <- list()
  
  for (age in ages) {
    cat("Processing age:", age, "\n")
    
    # Select cells from the current age group
    age_cells <- colnames(seu_obj)[seu_obj@meta.data[[age_col]] == age]
    
    if (length(age_cells) < min_cells) {
      cat(" Too few cells; skipped.\n")
      next
    }
    
    age_metadata <- seu_obj@meta.data[age_cells, ]
    
    seu_subset <- subset(seu_obj, cells = age_cells)
    Idents(seu_subset) <- seu_subset@meta.data[[brain_col]]
    
    markers <- FindMarkers(seu_subset,
                           ident.1 = "PFC",
                           ident.2 = "HIP",
                           min.pct = 0.1,
                           logfc.threshold = logfc_threshold,
                           test.use = "wilcox")
    
    if (nrow(markers) > 0) {
      markers$age <- age
      markers$protein <- rownames(markers)
      all_results[[as.character(age)]] <- markers
      
      cat("  Differential proteins identified:", nrow(markers), "个\n")
    } else {
      cat("  No differential proteins found.\n")
    }
  }
  
  # Combine results
  if (length(all_results) > 0) {
    combined_results <- bind_rows(all_results)
    return(combined_results)
  } else {
    warning("No differential proteins were found.")
    return(NULL)
  }
}

dep_results <- compute_brain_dep_by_age(seu)


# Prepare significant proteins for plotting
markers <- dep_results %>% 
  filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.05) %>%
  group_by(age) %>%
  mutate(type = ifelse(avg_log2FC > 0.5, "PFC-enriched", "HIP-enriched")) 

markers$type <- factor(markers$type, levels = c("PFC-enriched","HIP-enriched"))
markers$age <- factor(markers$age, levels = c("2M","14M","24M"))

p_jitter <- ggplot(markers, aes(x = age, y = avg_log2FC, color = type)) +
  geom_jitter(stroke = 0, size = 1, alpha = 0.8, 
              position = position_jitter(width = 0.3, seed = 42)) +
  scale_color_manual(values = c("PFC-enriched"="#9cca62", 
                                "HIP-enriched"="#00758b"),
                     name = "") +  
  labs(x = "", y = "log2FC (PFC vs. HIP)", title = "") +
  guides(color = guide_legend(override.aes = list(size = 2, shape = 19))) +
  theme_few() +  
  theme(
    panel.background = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.title   = element_text(size = 11),
    legend.position = c(0.86, 0.92),
    legend.key = element_blank(),
    legend.key.height = unit(0.3, "cm"),   
    legend.key.width = unit(0.5, "cm"),      
    legend.background = element_blank(),
    legend.text = element_text(size = 9, margin = margin(l = -1.5)),
    plot.title = element_text(size = 11, hjust = 0.5))


### Add known transcriptomic brain-region heterogeneity genes
# Tan, YL., Yuan, Y. & Tian, L. Microglial regional heterogeneity and its role in the brain. Mol Psychiatry 25, 351–367 (2020). https://doi.org/10.1038/s41380-019-0609-8
genes_fig1 <- c(
  "CD300E", "ITGAX", "CRYBB1", "CD22", "GPR84", 
  "S100A8", "S100A9", "SOCS3", "SLC2A5", "MRC1", 
  "PTPRC", "TGFB1", "ADGRE1", "MAFB", "SALL1", 
  "OLFML3","SMAD7","RNASE4","ITGB5","TYROBP",
  "TYROBP","RGS10",
  "TREM2", "P2RY13", "SLCO2B1", "GPR34", "SMAD3", 
  "CSF1R", "P2RY12", "TNFAR1", "CX3CR1",  
  "PRDX5", "HEXB", "LGNN", "SYNGR1", "BIN1"
)
genes_fig1 <- str_to_title(tolower(genes_fig1))   
genes_table1 <- c("Itgam","Hoxb8","Tim3","P2ry12","H2-Aa", "H2ab1", "H2-Eb1",
                  "H2-D1", "H2-K1")

# Grabert, K., Michoel, T., Karavolos, M. et al. Microglial brain region−dependent diversity and selective regional sensitivities to aging. Nat Neurosci 19, 504–516 (2016). https://doi.org/10.1038/nn.4222
# qPCR-validated genes
genes_qpcr <- c("Camp","H2ab1")
genes_fig4b <- c(
  "Clec4e", "Clec7a", "Fpr1", "Fpr2", "Fcnb", 
  "CD36", "Ccl44", "Itg4a", "Itgal", "Itgb2l", "Itgb7",
  "Ccl5", "Cxcl4", "Cxcl16", "Ccr1", "Cxcr2",
  "Nlcr5", "Zbp1", "Irf7", "Stat1", "Oasl1",
  "H2-Aa", "H2ab1", "H2-Eb1", "Cd74", "Ciita",
  "H2-D1", "H2-K1", "Camp", "Ngp", "Lyz1","Lyz2",
  "Aldoa", "Prdx5", "Cox5b", "Ogdh"
)

genes_fig4g <- c(
  "Aldoa", "Gpi1", "Hk1", "Pfkp", "Pgk1", "Pdha1",
  "Aco2", "Mdh1", "Sucla2", "Ogdh",
  "Ndufa1", "Ndufa12", "Ndufa13", "Ndufb6", "Ndufb7", "Ndufb10",
  "Uqcr11", "Uqcrq", "Uqcrh", "Uqcrb",
  "Atp5a1", "Atp5c", "Atp5d", "Atp5e", 
  "Cat", "Gpx4", "Gpx8", "Prdx2", "Prdx5", "Sod1", "Sod2",
  "Pparg", "Ppargc1a","Cox5a","Cox5b","Cox6a1","Cox7c"
)

genes_fig5 <- c(
  "Cd22", "Cd33", "Siglec5", "Siglece", 
  "Siglecg", "Siglech", 
  "Trem1","Trem2", "Trem3", "Cd300a", "Cd300lb", 
  "Cd300ld", "Cd300lf", "Cd200", "Cd200r1", 
  "Cd200r3", "Cd200r4", "Sirpa", "Sirpb1a", "Cd47"
)

genes_supp_fig5 <- c(
  "Cd86", "Tmem37", "Cd33", "Adora3", "Havcr2", 
  "Trem2", "Tlr12", "Cd79b", "Tnfrsf17", "Csf1r", 
  "Cx3cr1", "Entpd1", "Upk1b", "P2ry13", "Selplg", 
  "Fcrl1", "Tmem8c", "Slc7a7", "Tmem173", "Siglech", 
  "Cd53", "Fcgr3", "Tnfrsf1b", "Fcer1g", "Fcgr2b", 
  "Cd22", "Slamf9", "Fcgr4", "Cd52", "Clec4b1", 
  "Ifitm6", "Cd74", "Clec7a", "Cxcl16","Fcrls"
)

# de Haas AH, Boddeke HW, Biber K. Region-specific expression of immunoregulatory proteins on microglia in the healthy CNS. Glia. 2008 Jun;56(8):888-94. doi: 10.1002/glia.20663.
# FACS-validated
genes_facs <- c("Itgam","Cd40","Ptprc","Cd80","Cd86","Adgre","Trem2","Cxcr3","Ccr9")


brain_genes <- unique(c(genes_fig1,genes_table1,genes_qpcr,genes_fig4b,genes_fig4g,genes_fig5,genes_supp_fig5,genes_facs,"Cd68"))  
brain_proteins <- intersect(brain_genes, rownames(seu))  

markers_plot <- markers %>%
  mutate(
    region = case_when(
      protein %in% brain_proteins ~ "known",
      TRUE ~ "unknown"
    ))

label_df <- markers_plot %>%
  filter(region == "known") %>%
  group_by(age) %>%
  ungroup()

p_jitter <- ggplot(markers_plot, aes(x = age, y = avg_log2FC, color = type)) +
  geom_jitter(stroke = 0, size = 1, alpha = 0.8, 
              position = position_jitter(width = 0.3, seed = 42)) +
  geom_jitter(data = markers_plot %>% filter(protein %in% label_df$protein), aes(fill = type),
              stroke = 0.3, size = 1, shape = 21, color = "black", 
              position = position_jitter(width = 0.3, seed = 42),
              show.legend = FALSE) + 
  geom_text_repel(data = label_df, 
                  aes(label = toupper(protein)), 
                  color = "black", 
                  size = 2.5, 
                  position = position_jitter(width = 0.3, seed = 42),
                  segment.color = "gray50",
                  segment.size = 0.4,
                  segment.alpha = 0.7,
                  min.segment.length = 0.1,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = c("PFC-enriched"="#9cca62", 
                                "HIP-enriched"="#00758b"),
                     name = "") +  
  scale_fill_manual(values = c("PFC-enriched"="#9cca62", 
                               "HIP-enriched"="#00758b"),
                    name = "") + 
  labs(x = "", y = "log2FC (PFC vs. HIP)", title = "Regional-DEPs Across Ages") +
  guides(color = guide_legend(override.aes = list(shape = 15,
                                                  size = 4,
                                                  color = c("#9cca62", "#00758b"))),
         fill = "none") +
  theme_few() +  
  theme(panel.background = element_blank(),
        axis.text.y  = element_text(size = 10),
        axis.title   = element_text(size = 11),
        legend.position = c(0.05, 0.05),
        legend.justification = c(0, 0),
        legend.key = element_blank(),
        legend.key.height = unit(0.3, "cm"),  
        legend.key.width = unit(0.5, "cm"),   
        legend.background = element_blank(),  
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 11, hjust = 0.5))


# Bar plot showing DEP counts
count_data <- markers %>%
  group_by(age, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(age,
                  type = c("PFC-enriched", "HIP-enriched"),
                  fill = list(n = 0)) %>%
  mutate(y = ifelse(type == "PFC-enriched", n, -n),
         type = factor(type, levels = c("PFC-enriched", "HIP-enriched")))

p_bar <- ggplot(count_data, aes(x = age, y = y, fill = type)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  scale_fill_manual(values = c("PFC-enriched" = "#9cca62",
                               "HIP-enriched" = "#00758b"),
                    name = "") +
  scale_y_continuous(labels = abs) +
  labs(title = "", x = "", y = "Protein counts") +
  theme_few() +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = "right",
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))


# Violin plot for selected proteins
Idents(seu) <- "age"
seu@meta.data$brain_region <- factor(seu@meta.data$brain_region, levels = c("PFC","HIP"))

vln_features <- c("Cox6a1", "Ogdh", "Gpr34", "Itgb5")
feature_labels <- setNames(toupper(vln_features), vln_features)
VlnPlot(seu, 
        features = vln_features, 
        stack = T,
        split.by = "brain_region",
        flip = T, 
        add.noise = T,
        cols = c("#9cca62", "#00758b")) +
  facet_grid(feature ~ .,
             scales = "free_y",
             labeller = labeller(feature = as_labeller(feature_labels))) +
  theme_few() +
  theme(axis.title.x = element_blank(),             
        axis.text.x = element_text(angle = 0),       
        axis.text.y = element_text(face = "plain"),  
        strip.text.y.right = element_text(angle = 0, face = "plain"),   
        legend.position = "top",
        legend.title = element_blank())


#------------- 3.Phagoproteome total intensity ---------------------------------
ex_df <- read.csv("data/exo_genes_list.csv", header = F)
ex_genes <- ex_df$V1
gene_sets <- list(Phago = ex_genes)

calculate_geneset_counts <- function(seurat_obj, 
                                     geneset, 
                                     assay = "ALL", 
                                     layer = "counts") {
  mat <- GetAssayData(seurat_obj, 
                      assay = assay, 
                      layer = layer)
  existing_genes <- intersect(geneset, rownames(mat))
  
  if(length(existing_genes) == 0) {
    warning("No genes from the input gene set were found in the assay matrix.")
    return(rep(0, ncol(mat)))
  }
  colSums(mat[existing_genes, , drop = FALSE])
}

for(set_name in names(gene_sets)) {
  seu[[paste0(set_name, "_total")]] <- calculate_geneset_counts(seu, gene_sets[[set_name]])
}
seu@meta.data$log <- log10(seu@meta.data$Phago_total+1)

p_phago_total <- plotthis::ViolinPlot(seu[[]],x = "brain_region",
                                      y = "log", 
                                      sig_label = "p.signif",
                                      sig_labelsize = 4,
                                      step_increase = 0,
                                      legend.position = "none",
                                      xlab = "",
                                      ylab = "log10 (Intensity + 1)",
                                      title = "Phagoproteome Total Intensity",
                                      add_box = T,
                                      pt_alpha = 0.1,
                                      box_color = "black",
                                      box_width = 0.15,
                                      hide_ns = TRUE,
                                      palcolor = c("#9cca62", "#00758b"))+
  ggsignif::geom_signif(comparisons = list(c("PFC", "HIP")),
                        test = "t.test",color = "black",
                        y_position = 6.5,
                        map_signif_level = TRUE)+
  coord_cartesian(ylim = c(-0.5, 7.5), expand = FALSE) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(~ factor(age, levels = c("2M", "14M", "24M")))


#------------------ 4. Calculate UCell pathway scores --------------------------------------------
## Build GO-based gene sets from msigdbr
GOCC_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:CC") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

GOBP_list <- msigdbr(species = "Mus musculus",
                     category = "C5",
                     subcategory = "GO:BP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

phagocytosis_genes <- unique(c(GOBP_list$GOBP_PHAGOCYTOSIS,
                               GOBP_list$GOBP_PHAGOCYTOSIS_RECOGNITION))
lysosome_genes <- unique(c(GOBP_list$GOBP_LYSOSOMAL_TRANSPORT,
                           GOBP_list$GOBP_LYSOSOMAL_PROTEIN_CATABOLIC_PROCESS,
                           GOBP_list$GOBP_LYSOSOMAL_MEMBRANE_ORGANIZATION,
                           GOBP_list$GOBP_LYSOSOMAL_LUMEN_ACIDIFICATION,
                           GOCC_list$GOCC_LYSOSOMAL_LUMEN))

pathway_gene_sets <- list(lysosome = lysosome_genes,
                          phagocytosis = phagocytosis_genes)

## Calculate UCell pathway scores 
seu_obj <- AddModuleScore_UCell(seu,
                                features = pathway_gene_sets,
                                name = "_UCell",
                                assay = "ALL",
                                slot = "counts")

p_phagocytosis <- plotthis::ViolinPlot(seu_obj[[]],x = "brain_region",
                                       y = "phagocytosis_UCell", 
                                       sig_label = "p.signif",
                                       sig_labelsize = 4,
                                       step_increase = 0,
                                       legend.position = "none",
                                       xlab = "",
                                       ylab = "UCell Score",
                                       title = "Phagocytosis Pathway Score",
                                       add_box = T,
                                       pt_alpha=0.1,
                                       box_color = "black",
                                       box_width = 0.15,
                                       hide_ns = TRUE,
                                       palcolor = c("#9cca62", "#00758b"))+
  ggsignif::geom_signif(comparisons = list(c("PFC", "HIP")),
                        test = "t.test",color = "black",
                        y_position = 0.22,
                        map_signif_level = TRUE)+
  coord_cartesian(ylim = c(0, 0.25), expand = FALSE) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(~ factor(age, levels = c("2M", "14M", "24M")))


p_lysosome <- plotthis::ViolinPlot(seu_obj[[]],x = "brain_region",
                                   y = "lysosome_UCell", 
                                   sig_label = "p.signif",
                                   sig_labelsize = 4,
                                   step_increase = 0,
                                   legend.position = "none",
                                   xlab = "",
                                   ylab="Lysosome Pathway Score",
                                   add_box = T,
                                   pt_alpha=0.1,
                                   box_color = "black",
                                   box_width = 0.15,
                                   hide_ns = TRUE,
                                   palcolor = c("#9cca62", "#00758b"))+
  ggsignif::geom_signif(comparisons = list(c("PFC", "HIP")),
                        test = "t.test",color = "black",
                        y_position = 0.23,
                        map_signif_level = TRUE)+
  coord_cartesian(ylim = c(0, 0.27), expand = FALSE) +
  theme(strip.background = element_rect(color = "beige", fill="beige"),)+
  facet_grid(~ factor(age, levels = c("2M", "14M", "24M")))


