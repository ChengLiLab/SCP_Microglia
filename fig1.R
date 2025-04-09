#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: SCP dataset analysis                                            ##
## Data: 2025.03.26                                                                ##
## Author: Bijia Chen                                                              ##
#####################################################################################

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
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(GseaVis)
library(cowplot)
library(eulerr)
library(ggsci)
source("code/scp_utils.R")

# 0.Load Files -------------------------------------------------------------------------------
# SCP
scp_seu <- readRDS("data/SCP_MG_Seurat.rds")
DefaultAssay(scp_seu) <- "SCP"
scp_count <- as.data.frame(GetAssayData(scp_seu, layer = "counts"))
scp_data <- as.data.frame(LogNormalize(data = scp_count))

protein_data <- read.csv("data/SCP_MG_All_Sample_metadata.csv", header = TRUE) 
protein_data <- subset(protein_data,celltype!="CD45neg")

protein_seu <- readRDS("data/SCP_MG_All_Sample_Seurat.rds")
protein_seu <- protein_seu %>%
  subset(celltype != "CD45neg")
cell_seu <- protein_seu %>%
  subset(mass=="1c")
cell_seu <- cell_seu[rowSums(cell_seu@assays$SCP@layers$counts) > 0, ]
cell_count <- as.data.frame(GetAssayData(cell_seu, layer = "counts"))
cell_data <- as.data.frame(LogNormalize(data = cell_count))

cells_seu <- protein_seu %>%
  subset(mass=="20c")
cells_seu <- cells_seu[rowSums(cells_seu@assays$SCP@layers$counts) > 0, ]
cells_count <- as.data.frame(GetAssayData(cells_seu, layer = "counts"))
cells_data <- as.data.frame(LogNormalize(data = cells_count))

# 1.Quantified proteins/peptides count-----------------------------------------------------------------------------------
stats <- protein_data %>%  
  group_by(mass) %>%  
  summarize(n = n(),
            mean = mean(nFeature_SCP),
            sd = sd(nFeature_SCP),
            min = min(nFeature_SCP),
            max = max(nFeature_SCP)) %>%
  mutate(mass=factor(mass,levels=c("Blank","1c","20c")))
protein_data$mass <- factor(protein_data$mass,levels=c("Blank","1c","20c"))

p_num <- ggplot() +  
  geom_errorbar(stats, mapping=aes(x = mass, ymin = mean-sd, ymax = mean+sd),
                width = 0.4, 
                color = c("black"), 
                linewidth = 0.4) + 
  geom_bar(stats, mapping=aes(x = mass, y = mean, fill = mass),
           stat = "identity", 
           position = position_dodge()) +  
  scale_fill_manual(values = c("Blank" = "#b0afb7","1c" = "#da6968", "20c" = "#95a4cb")) +
  geom_jitter(protein_data,mapping=aes(x = mass, y = nFeature_SCP), 
              size = 0.1, 
              color = "black",
              alpha = 0.2,
              width = 0.3) +  
  labs(title = "",  
       x = "",  
       y = "Protein Groups") +  
  theme_few() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=9))

# peptides count
stats_p <- protein_data %>%  
  group_by(mass) %>%  
  summarize(n = n(),
            mean = mean(nFeatures_pep),
            sd = sd(nFeatures_pep)) %>%
  mutate(mass=factor(mass,levels=c("Blank","1c","20c")))

p_num <- ggplot() +  
  geom_errorbar(stats_p, mapping=aes(x = mass, ymin = mean-sd, ymax = mean+sd),
                width = 0.4, 
                color = c("black"), 
                linewidth = 0.4) + 
  geom_bar(stats_p, mapping=aes(x = mass, y = mean, fill = mass),
           stat = "identity",  
           position = position_dodge()) +  
  scale_fill_manual(values = c("Blank" = "#b0afb7","1c" = "#da6968", "20c" = "#95a4cb")) +
  geom_jitter(protein_data,mapping=aes(x = mass, y = nFeatures_pep), 
              size = 0.1, 
              color = "black",
              alpha = 0.2,
              width = 0.3) +  
  labs(title = "",  
       x = "",  
       y = "Peptides") +  
  theme_few() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=9))


# 2.Venn plot of protein groups -----------------------------------------------------------------------------------
overlap <- intersect(rownames(cell_count), rownames(cells_count))
df <- c("1c" = nrow(cell_count)-length(overlap),
        "20c" = nrow(cells_count)-length(overlap), 
        "1c&20c" = length(overlap))
p <- plot(euler(df), 
          quantities = list(type = "counts",cex=1),          
          edges = list(col = "lightgray", lex = 2, lwd=0.6), 
          fills = list(fill = c("#da6968", "#95a4cb"),alpha=0.8),
          main = list(label = "Protein Group Overlap",cex=1.5))

# 3.Downsampling analysis -----------------------------------------------------------------------------------
set.seed(111)   
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

index <- head(which(saturation_results$Identified_Proteins >= 4000), n = 1)  
point <- saturation_results$Cell_Number[index]  

p <- ggplot(saturation_results, aes(x = Cell_Number, y = Identified_Proteins)) +  
  geom_line(color = "gray", linewidth = 0.5) +  
  geom_point(color = "#da6968", size = 0.5) +  
  geom_vline(xintercept = saturation_results$Cell_Number[point], linetype = "dashed", color = "black", alpha = 0.6) +  
  labs(title = "Protein Saturation Curve",  
       x = "Number of Cells",  
       y = "Number of Proteins") +  
  theme_few() 

# 4.Protein completeness -----------------------------------------------------------------------------------
completeness_df <- data.frame(  
  completeness = rowSums(scp_count > 0) / ncol(scp_count) * 100  
) 

p_com <- ggplot(completeness_df, aes(x = completeness)) +  
  geom_histogram(alpha = 0.8,linewidth=0.2,color="black",binwidth = 5,fill = "#da6968") +  
  labs(title = "",  
       x = "Completeness (%)",  
       y = "Protein Groups") +  
  theme_few()  


# 5.Sample correlation -----------------------------------------------------------------------------------
tmp.meta <- scp_seu@meta.data
tmp.meta$age <- factor(tmp.meta$age,levels = c("2M","14M","24M"))
tmp.levels <- levels(tmp.meta$age)

tmp.data <- scp_data
tmp.data <- tmp.data %>%
  mutate(across(everything(), ~ ifelse(. == 0, NA, .)))
tmp.res1 <- lapply(seq_along(tmp.levels), function(ii){
  x <- tmp.levels[ii] 
  tmp.id <- tmp.meta %>%
    filter(mass == "1c") %>%
    filter(age == x) %>%
    rownames() 
  tmp.df <- data.frame(avg_PCC = apply(cor(tmp.data[,tmp.id], use = "pairwise.complete.obs"), 1, mean, na.rm = TRUE),stringsAsFactors = F)
  return(tmp.df)
})
tmp.res1 <- Reduce(f = rbind,tmp.res1) 

tmp.meta <- protein_seu@meta.data
tmp.meta$age <- factor(tmp.meta$age,levels = c("2M","14M","24M"))
tmp.levels <- levels(tmp.meta$age)
tmp.data <- cells_data
tmp.data <- tmp.data %>%
  mutate(across(everything(), ~ ifelse(. == 0, NA, .)))
tmp.res2 <- lapply(seq_along(tmp.levels), function(ii){
  x <- tmp.levels[ii] 
  tmp.id <- tmp.meta %>%
    filter(mass == "20c") %>%
    filter(age == x) %>%
    rownames() 
  tmp.df <- data.frame(avg_PCC = apply(cor(tmp.data[,tmp.id], use = "pairwise.complete.obs"), 1, mean, na.rm = TRUE),stringsAsFactors = F)
  return(tmp.df)
})
tmp.res2 <- Reduce(f = rbind,tmp.res2) 

tmp.pcc <- merge(rbind(tmp.res1,tmp.res2),tmp.meta,by="row.names")

tmp.bar.plot <- tmp.pcc %>%
  group_by(age,mass) %>%
  summarise(median_prot_value = median(avg_PCC, na.rm = TRUE),
            avg_prot_value = mean(avg_PCC, na.rm = TRUE),
            sd_prot_value = sd(avg_PCC, na.rm = TRUE))

p <- ggplot()+
  geom_errorbar(data = tmp.bar.plot,
                aes(x = age,
                    ymin = avg_prot_value-sd_prot_value,
                    ymax = avg_prot_value+sd_prot_value),
                width = 0.4, 
                color = c("black"), 
                linewidth = 0.4)+
  geom_bar(data = tmp.bar.plot,
           aes(age,
               avg_prot_value,
               fill=age),
           stat = "identity",
           width = 0.8)+
  geom_jitter(tmp.pcc,mapping=aes(x = age, y = avg_PCC), 
              size = 0.1, 
              color = "black",
              alpha = 0.2,
              width = 0.3) + 
  scale_fill_manual(values = c("#ceeb9c","#5bb6a9","#5e4fa2"))+
  ylab("Mean Pearson correlation")+
  xlab(NULL)+
  scale_y_continuous(breaks = seq(0,1,by=0.2),limits = c(0,1.1),expand = c(0,0))+
  theme_few() +
  theme(legend.position = "none",
        panel.grid = element_line(colour = "grey85",linewidth = 0.25,linetype = "dashed"),
        axis.line = element_line(linewidth  = 0.5),
        axis.ticks = element_line(size = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1,vjust = 1)) +
  facet_wrap(~ mass)
