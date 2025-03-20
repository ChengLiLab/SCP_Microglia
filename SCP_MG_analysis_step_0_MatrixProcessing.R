
##----0.load package--------
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
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(MSnID)
library(ggsci)
source("scp_utils.R")

##-------1. SCP data loading ----------

###-------load long form data--------
####------pdata----------

### For large datasets , long form data is the only choice
runLoadData <- function(tmp.path){
  tmp.all <- read_tsv(file = tmp.path)
  return(tmp.all)
}

runProcess.pdata <- function(tmp.all){
  tmp.pdata <- tmp.all %>%
    group_by(R.FileName) %>% 
    dplyr::mutate(nFeatures_pep = n()) %>% 
    ungroup() %>% 
    dplyr::select(R.FileName) %>%
    dplyr::rename(id = R.FileName) %>%
    unique() %>%
    dplyr::mutate(loading_date = gsub(pattern = "_.*",
                                      replacement = "",id)) %>%
    dplyr::mutate(prepare_date = gsub(pattern = ".*min_",
                                      replacement = "",id)) %>%
    dplyr::mutate(prepare_date = gsub(pattern = "_.*",
                                      replacement = "",prepare_date)) %>%
    dplyr::mutate(chromatogram_infor=gsub(pattern="[0-9]{8}_",
                                          replacement = "",id)) %>%
    dplyr::mutate(chromatogram_infor = gsub(pattern="_.*",
                                            replacement = "",chromatogram_infor)) %>%
    dplyr::mutate(LC_gradient = gsub(pattern=".*Io[0-9]{2}_",replacement = "",id)) %>%
    dplyr::mutate(LC_gradient = gsub(pattern="_.*",
                                   replacement = "",LC_gradient)) %>%
    dplyr::mutate(method = "SCP") %>%
    dplyr::mutate(celltype = gsub(".*HIP_|.*PFC_|_1c.*$|_20c.*$|_Blank.*$",replacement="",id)) %>%
    dplyr::mutate(age = gsub(pattern=".*_[0-9]{4}_",replacement="",id)) %>%
    dplyr::mutate(age = gsub(pattern="_.*",replacement="",age)) %>%
    dplyr::mutate(gender = gsub(pattern=".*M_",replacement="",id)) %>%
    dplyr::mutate(gender = gsub(pattern="_.*",replacement="",gender)) %>%
    dplyr::mutate(brain_region=gsub(pattern=".*male_",replacement="",id)) %>%
    dplyr::mutate(brain_region=gsub(pattern="_.*",replacement="",brain_region)) %>%
    dplyr::mutate(sample_type="SCP") %>%
    dplyr::mutate(mass = gsub(pattern=".*MG_|.*CD45neg_",
                              replacement = "",id)) %>%
    dplyr::mutate(mass = gsub(pattern="_.*",
                              replacement = "",mass)) %>%
    dplyr::mutate(samle_id = gsub(pattern=".*(Blank|1c|20c)_",
                                  replacement = "",id))
  return(tmp.pdata)
}

tmp.dir <- "../20241222_analysis/SCP"
tmp.files <- list.files(path = tmp.dir,full.names = T)
tmp.files <- grep(pattern = ".txt",invert = T,x = tmp.files,value = T)

tmp.path.list <- lapply(seq_along(tmp.files),function(x){
  ii <- tmp.files[x]
  jj <- grep(pattern = "\\(Normal\\).tsv",list.files(ii),value = T)
  if (is_empty(jj)) {
    tmp.path <- jj
    cat(paste0(ii,' is empty'),sep = "\n")
  }
  else{
    tmp.path <- paste0(ii,"/",jj)
  }
  return(tmp.path)
})

tmp.path.list <- unlist(tmp.path.list)

tmp.part.list <- lapply(seq_along(tmp.path.list), 
                        function(xx){
                          tmp.path <- tmp.path.list[xx]
                          tmp.part <- runLoadData(tmp.path)
                          return(tmp.part)
                        })
tmp.all <- Reduce(rbind,tmp.part.list)
tmp.all <- tmp.all[!is.na(tmp.all$R.FileName),]
#tmp.all <- readRDS("tmp.all_20241224.rds")

tmp.pdata <- runProcess.pdata(tmp.all)
tmp.pdata <- tmp.pdata %>%
  unique() %>%
  dplyr::mutate(mass = case_when(
     str_detect(id, "CD45neg_1c") ~ "CD45neg_1c",
     str_detect(id, "CD45neg_20c") ~ "CD45neg_20c",
    str_detect(id, "Debris") ~ "Blank",
    TRUE ~ mass,
  ))
table(tmp.pdata$mass)
saveRDS(tmp.pdata,"tmp.pdata_20241224.rds")

###-------explore -----------
tmp.pdata$age <- factor(tmp.pdata$age,
                        levels = c("2M","14M","24M"))

table(tmp.pdata$mass)

ttt <- tmp.pdata %>%
  dplyr::filter(mass == "1c")
table(ttt$brain_region,
      ttt$prepare_date,
      ttt$age,
      ttt$celltype)

table(ttt$age)
table(ttt$chromatogram_infor,ttt$age)

ttt <- tmp.pdata %>%
  dplyr::filter(age=="2M") 
table(ttt$brain_region,ttt$prepare_date)
grep("MG",ttt$brain_region)
table(ttt$brain_region)

###-----merge export-------
tmp.path <- "../20241222_analysis/Produced.Date.txt"
tmp.data <- read_delim(file = tmp.path)
x <- tmp.data$File.names
x <- gsub(pattern = "\\.raw",replacement = "",x)
x
tmp.data$File.names <- x
colnames(tmp.data)[[5]] <- "id"

setdiff(tmp.pdata$id,tmp.data$id)

head(tmp.pdata)
tmp.pdata <- left_join(tmp.pdata,
                       tmp.data,by="id")
head(tmp.pdata)

write.table(tmp.pdata,
            file = "res/txt/sample_infor_20241207.txt",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)

####-------- fdata------------
runProcess.fdata <- function(tmp.all){
  tmp.fdata <- tmp.all %>%
    dplyr::select(PG.ProteinGroups,
                  PG.Genes,
                  PG.UniProtIds,
                  PG.ProteinNames,
                  PG.ProteinDescriptions) %>%
    unique() %>%
    dplyr::mutate(is_contaminate = str_detect(PG.ProteinGroups, "^CON_")) %>%
    dplyr::filter(!is.na(PG.Genes)) %>%
    dplyr::filter(is_contaminate==FALSE) %>%
    dplyr::select(-is_contaminate)
  return(tmp.fdata)
}

tmp.fdata <- runProcess.fdata(tmp.all) %>% unique()
saveRDS(tmp.fdata,"tmp.fdata_20241223.rds")
####--------extract LFQ --------------
### 
head(tmp.fdata)

tmp.fdata[duplicated(tmp.fdata$PG.ProteinGroups),]
tmp.fdata[is.na(tmp.fdata$PG.ProteinGroups),]

runProcess.matrix <- function(tmp.all,tmp.fdata){
  tmp.mat <- tmp.all %>%
    dplyr::select(R.FileName,PG.ProteinGroups,PG.Quantity) %>%
    unique()
  order <- tmp.mat$R.FileName %>% unique()
  tmp.mat <- merge(tmp.mat,
                   tmp.fdata,
                   by="PG.ProteinGroups")
  tmp.mat.out <- tmp.mat %>%
    dplyr::select(R.FileName,PG.Genes,PG.Quantity) %>%
    pivot_wider(names_from = R.FileName, 
                values_from = PG.Quantity) %>% 
    column_to_rownames("PG.Genes")
  tmp.mat.out <- subset(tmp.mat.out,select=order)
  tmp.mat.out <- myRemoveNA(tmp.mat.out)
  return(tmp.mat.out)
}

tmp.mat.out <- tmp.all %>% 
  runProcess.matrix(tmp.fdata = tmp.fdata)  
rownames_to_column("genes")

####--------extract IBAQ---------------

runProcess.ibaq <- function(tmp.all,tmp.fdata){
  tmp.mat <- tmp.all %>%
    dplyr::select(R.FileName,PG.ProteinGroups,PG.IBAQ) %>%
    unique()
  tmp.mat <- merge(tmp.mat,
                   tmp.fdata,
                   by="PG.ProteinGroups")
  tmp.mat.ibaq <- tmp.mat 
  return(tmp.mat.ibaq)
}

tmp.ibaq.out <- runProcess.ibaq(tmp.all = tmp.all,
                                tmp.fdata = tmp.fdata)

saveRDS(tmp.ibaq.out,file = "res/r/SCP_MG_3440_ibaq_zhr_20241207.rds")


####-------- combine to Seurat------

tmp.fdata <- tmp.fdata %>%
  unique() %>%
  column_to_rownames("PG.Genes")

low_quality <- c("20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_A22",
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_B21",
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_B23",    
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_B24",
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_D24",
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_I18",
                 "20241023_Io02_10min_1015_2M_male_PFC_MG_1c_P2_P20",
                 "20241005_Io01_10min_0911_14M_male_PFC_MG_1c_P1_H22",
                 "20241001_Io01_10min_0815_2M_male_PFC_MG_1c_P1_P12",
                 "20241104_Io03_10min_0913_14M_male_PFC_MG_1c_P2_I15",
                 "20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_F11",
                 "20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_M8",
                 "20241115_Io03_10min_1105_24M_male_HIP_MG_1c_P2_F9") # samples with incomplete LC_gradient
tmp.mat.out <- tmp.mat.out[,-which(colnames(tmp.mat.out)%in%low_quality)]

doublet <- c('20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_N1', 
             '20240927_Io01_10min_0911_14M_male_HIP_MG_1c_P1_K5', 
             '20241005_Io01_10min_0911_14M_male_PFC_MG_1c_P1_F7', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_A10', 
             '20240927_Io01_10min_0911_14M_male_HIP_MG_1c_P1_E8', 
             '20241104_Io03_10min_0913_14M_male_PFC_MG_1c_P2_I21', 
             '20241020_Io02_10min_1015_2M_male_PFC_MG_1c_P1_F22', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_P10', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_I4', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_L12', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_L13', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_M13', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_N19', 
             '20241027_Io02_10min_0822_24M_male_PFC_MG_1c_P1_C10', 
             '20240918_Io01_10min_0815_2M_male_HIP_MG_1c_P1_C22')
tmp.mat.out <- tmp.mat.out[,-which(colnames(tmp.mat.out)%in%doublet)]

tmp.mat.out <- tmp.mat.out %>%
  column_to_rownames("genes")

tmp.meta <- tmp.pdata %>%
  column_to_rownames("id") 

tmp.meta <- tmp.meta[-which(rownames(tmp.meta)%in%low_quality),]
tmp.meta <- tmp.meta[-which(rownames(tmp.meta)%in%doublet),]



head(setdiff(rownames(tmp.meta),colnames(tmp.mat.out)))
nrow(tmp.meta)
ncol(tmp.mat.out)
head(rownames(tmp.meta))
head(colnames(tmp.mat.out))
tail(colnames(tmp.mat.out))
tail(rownames(tmp.meta))

setdiff(rownames(tmp.meta),colnames(tmp.mat.out))
setdiff(colnames(tmp.mat.out),rownames(tmp.meta))

###remove histone, antibodies, Amylase, trypsin and Albumin
tmp.mat.out <- tmp.mat.out[rowSums(tmp.mat.out)!=0,]

pg_to_rm <- grep("Igkv|Ighv|H2ac|H3-|Hist1|H2bc|H1-|H4c|Hist2|H2az|Try10|Alb|Amy",
                 rownames(tmp.mat.out))
tmp.mat.out <- tmp.mat.out[-pg_to_rm,]
saveRDS(tmp.mat.out,"tmp.mat.out_6626_rmig_His_try_alb_20241224.rds")

#---------Matrix proccessing
tmp.mat.out <- tmp.mat.out %>% 
  rownames_to_column("genes")
df <- as.data.frame(tmp.mat.out)
split_groups <- strsplit(df$genes, ";")
has_intersection <- function(current_group, other_groups) {
  any(map_lgl(other_groups, ~ {
    if (!identical(current_group, .x)) {
      length(intersect(current_group, .x)) > 0
    } else {
      FALSE
    }
  }))
}
df <- df %>%
  mutate(
    label = map_lgl(split_groups, ~ has_intersection(.x, split_groups))
  ) %>% 
  mutate(label = if_else(label, "duplicate", "unique"))
df.filt <- df %>%
  separate_rows(genes, sep = ";") %>%  
  group_by(genes) %>%
  mutate(type=if_else(n() > 1, 1, 0)) %>% 
  ungroup() %>% 
  mutate(label_type=paste0(label,"_",type)) %>% 
  subset(label_type != "duplicate_0") %>% 
  dplyr::select(-c("label_type","label","type"))
tmp.mat.out.max <- df.filt %>%
  group_by(genes) %>%
  summarise(across(where(is.numeric), max), .groups = 'drop') %>% 
  filter(genes != "") %>% 
  column_to_rownames("genes")
nrow(tmp.mat.out)
nrow(tmp.mat.out.max)
head(tmp.mat.out.max)
saveRDS(tmp.mat.out.max,"20241228_pocessed_exp_mat_rmig_his_try_amy_krt.rds")

# create seurat object
seu <- CreateSeuratObject(counts = tmp.mat.out.max,
                          meta.data = tmp.meta,
                          assay = "SCP")
#head(seu[[]])

write.csv(seu@meta.data,"res/r/20241223_metadata.csv")
# To add cell level information, 
# add to the Seurat object. If adding feature-level metadata, add to the Assay object (e.g. object[["RNA"]])

seu[["SCP"]] <- AddMetaData(object = seu[["SCP"]],
                            metadata = tmp.fdata)

### Warning this slot contain 20c data etc,
ncol(seu)
nrow(seu)

