# 0_MatrixProcessing_QC.R
# Purpose: preprocess label-free SCP long-format search reports, build a gene-level quantity matrix, create a Seurat object, and perform cell-level QC.

##----0.load package--------
library(tidyverse)
library(cowplot)
library(Seurat)
library(ggsci)

##-------1. Matrix Processing ----------
tmp.dir <- "data/SCP"
tmp.files <- list.files(path = tmp.dir,full.names = T)

tmp.path.list <- lapply(seq_along(tmp.files),function(x){
  ii <- tmp.files[x]
  jj <- grep(pattern = "\\SCP_Report_Normal \\(Normal\\).tsv$",list.files(ii),value = T)
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
tmp.path.list

runLoadData <- function(tmp.path){
  tmp.all <- read_tsv(file = tmp.path)
  return(tmp.all)
}

tmp.part.list <- lapply(seq_along(tmp.path.list), 
                        function(xx){
                          tmp.path <- tmp.path.list[xx]
                          tmp.part <- runLoadData(tmp.path) %>% 
                            dplyr::select(-any_of(c("R.ModifiedSequencesIdentified",
                                                    "R.TotalIonCurrent (MS1)",
                                                    "EG.TotalQuantity (Settings)")))
                          return(tmp.part)
                        })

common_cols <- Reduce(intersect, lapply(tmp.part.list, colnames))
tmp.all <- Reduce(rbind, lapply(tmp.part.list, function(x) x[, common_cols]))
tmp.all <- tmp.all[!is.na(tmp.all$R.FileName),]
tmp.all$R.FileName <- ifelse(!grepl("HIP|PFC", tmp.all$R.FileName),  
                             gsub("(_male)(?=.*_MG_)", "\\1_NA", tmp.all$R.FileName, perl = TRUE),
                             tmp.all$R.FileName)  


runProcess.pdata <- function(tmp.all){
  tmp.pdata <- tmp.all %>%
    dplyr::select(R.FileName) %>%
    dplyr::rename(id = R.FileName) %>%
    group_by(id) %>% 
    summarise() %>% 
    ungroup() %>% 
    dplyr::mutate(loading_date = gsub(pattern = "_.*",
                                      replacement = "",id)) %>%
    dplyr::mutate(prepare_date = gsub(pattern = ".*min_",
                                      replacement = "",id)) %>%
    dplyr::mutate(prepare_date = gsub(pattern = "_.*",
                                      replacement = "",prepare_date)) %>%
    dplyr::mutate(column_id=gsub(pattern="[0-9]{8}_",
                                 replacement = "",id)) %>%
    dplyr::mutate(column_id = gsub(pattern="_.*",
                                   replacement = "",column_id)) %>%
    dplyr::mutate(LC_gradient = gsub(pattern=".*Io[0-9]{2}_",replacement = "",id)) %>%
    dplyr::mutate(LC_gradient = gsub(pattern="_.*",
                                     replacement = "",LC_gradient)) %>%
    dplyr::mutate(method = "SCP") %>%
    dplyr::mutate(celltype = gsub(".*NA_|.*HIP_|.*PFC_|_1c.*$|_20c.*$|_Blank.*$|_Debris.*$",replacement="",id)) %>%
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
    dplyr::mutate(samle_id = gsub(pattern=".*(Blank|1c|20c|Debris)_",
                                  replacement = "",id))
  return(tmp.pdata)
}

tmp.pdata <- runProcess.pdata(tmp.all)
table(tmp.pdata$celltype)
tmp.pdata <- tmp.pdata %>% subset(celltype != "CD45neg")
table(tmp.pdata$mass)
# 1c 
# 3317 


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


####--------intensity matrix --------------
myRemoveNA <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}
runProcess.matrix <- function(tmp.all,tmp.fdata){
  tmp.mat <- tmp.all %>%
    dplyr::select(R.FileName,
                  PG.Quantity,PG.Genes) %>%
    group_by(PG.Genes, R.FileName) %>%    
    summarise(
      PG.Quantity = sum(PG.Quantity, na.rm = TRUE),  
      .groups = "drop"                               
    ) %>% 
    unique()
  order <- tmp.mat$R.FileName %>% unique()
  tmp.mat.out <- tmp.mat %>%
    filter(
      !is.na(PG.Genes),          
      PG.Genes != "",            
      !grepl("^\\s*$", PG.Genes) 
    ) %>% 
    pivot_wider(names_from = R.FileName, 
                values_from = PG.Quantity) %>% 
    column_to_rownames("PG.Genes")
  tmp.mat.out <- subset(tmp.mat.out,select=order)
  tmp.mat.out <- myRemoveNA(tmp.mat.out)
  return(tmp.mat.out)
}
tmp.mat.out <- tmp.all %>% 
  runProcess.matrix(tmp.fdata = tmp.fdata) 

tmp.meta <- tmp.pdata %>%
  mutate(id = factor(id, levels = colnames(tmp.mat.out))) %>%
  arrange(id) %>%                                  
  mutate(id = as.character(id)) %>%            
  column_to_rownames("id") 

####--------filtering ---------------------------
# samples with incomplete LC_gradient
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
                 "20241115_Io03_10min_1105_24M_male_HIP_MG_1c_P2_F9",
                 "20241014_Io01_10min_0911_14M_male_HIP_MG_Blank_P2_A4",
                 "20241014_Io01_10min_0911_14M_male_HIP_MG_Blank_P2_A6",
                 "20240918_Io01_10min_0815_2M_male_HIP_MG_Blank_P1_A8")
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
tmp.meta <- tmp.meta[-which(rownames(tmp.meta) %in% c(low_quality, doublet)), ]

tmp.mat.out <- tmp.mat.out[, which(colnames(tmp.mat.out) %in% rownames(tmp.meta))]

# remove histone, antibodies, Amylase, trypsin and Albumin
tmp.mat.out <- tmp.mat.out[rowSums(tmp.mat.out)>0,]
pg_to_rm <- grep("Igkv|Ighv|Try10|Alb|Krt|Amy",
                 rownames(tmp.mat.out))
tmp.mat.out <- tmp.mat.out[-pg_to_rm,]

valid_genes <- grep("^[A-Za-z0-9;]+$", rownames(tmp.mat.out), value = TRUE)
tmp.mat.out <- tmp.mat.out[ which(rownames(tmp.mat.out) %in% valid_genes),]


# multiple protein groups 
df <- tmp.mat.out %>% 
  rownames_to_column("Gene") %>% 
  data.frame()
df$Gene <- gsub(";$", "", df$Gene)

multi_gene_rows <- df %>%
  filter(grepl(";", Gene)) %>%
  separate_rows(Gene, sep = ";") %>%
  unique()
multi_gene_list <- multi_gene_rows %>%
  dplyr::select("Gene")
single_gene_rows <- df %>%
  filter(!grepl(";", Gene))
single_gene_list <- single_gene_rows %>%
  dplyr::select("Gene")

common_elements <- intersect(multi_gene_list$Gene, single_gene_list$Gene)
common_elements
sum_rows <- multi_gene_rows[which(multi_gene_rows$Gene %in% common_elements), ]

tmp.mat.out.sum <- rbind(single_gene_rows, multi_gene_rows[which(multi_gene_rows$Gene %in% common_elements), ]) %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  column_to_rownames("Gene")
colnames(tmp.mat.out.sum) <- sub("^X", "", colnames(tmp.mat.out.sum))
identical(colnames(tmp.mat.out.sum), colnames(tmp.mat.out))

rownames(df) <- NULL
tmp.mat.out_filtered <- df %>%
  filter(grepl(";", Gene)) %>%
  as.data.frame() %>% 
  column_to_rownames("Gene") %>%
  anti_join(sum_rows[,-1])
tmp.mat.out_filtered$Gene <- trimws(sapply(strsplit(rownames(tmp.mat.out_filtered), ";"), `[`, 1))
tmp.mat.out_filtered <- tmp.mat.out_filtered %>%
  group_by(Gene) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  column_to_rownames("Gene")
colnames(tmp.mat.out_filtered) <- sub("^X", "", colnames(tmp.mat.out_filtered))
identical(colnames(tmp.mat.out.sum), colnames(tmp.mat.out_filtered))
tmp.mat.out.sum <- rbind(tmp.mat.out.sum,tmp.mat.out_filtered)

#----------Seurat Object
head(setdiff(rownames(tmp.meta),colnames(tmp.mat.out.sum)))

seu <- CreateSeuratObject(counts = tmp.mat.out.sum,
                          meta.data = tmp.meta,
                          assay = "SCP")
ncol(seu)
# 3291
nrow(seu)
# 4349


##-------2. QC ----------------------------------------------------------------
###------define function----------
numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  ### this makes the line move above the line
  return(sum(myVector >=(xPts*slope+b)))
}


calculate_cutoff <- function(inputVector = abs(tmp.data.plot$nFeature_RNA),tmp_ylab, drawPlot=TRUE,...){
  
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),
                        myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  
  
  ## second slope add by wlt to get soft cutoff
  slope_final <- (y_cutoff - min(inputVector))/xPt
  xPt_final <- floor(optimize(numPts_below_line,
                              lower=1,
                              upper=xPt,
                              myVector= inputVector,
                              slope=slope_final)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff_final <- inputVector[xPt_final] #The y-value at this x point. This is our cutoff.
  
  
  
  if(drawPlot){  #if TRUE, draw the plot
    
    tmp.data.plot <- data.frame(rank = 1:length(inputVector),
                                PCC = inputVector,
                                stringsAsFactors = F)
    
    b <- y_cutoff-(slope* xPt)
    b_final <- y_cutoff_final-(slope_final * xPt_final)
    
    p <- ggplot(tmp.data.plot,aes(rank,PCC))+
      geom_line(size = 1)+
      geom_point(x = 1,y = inputVector[1],color = "skyblue",size = 3)+
      geom_point(x = length(inputVector),y = inputVector[length(inputVector)],color = "skyblue",size = 3 )+
      geom_path(data = data.frame(x = c(1,length(inputVector)),
                                  y = c(inputVector[1],inputVector[length(inputVector)]),
                                  stringsAsFactors = F),
                mapping = aes(x,y),
                size = 1,color = "skyblue")+
      geom_path(data = data.frame(x = c(1,xPt),
                                  y = c(inputVector[1],inputVector[xPt]),
                                  stringsAsFactors = F),
                mapping = aes(x,y),
                size = 1,color = "skyblue")+
      geom_abline(intercept = b,slope = slope,color = "skyblue",size = 1)+
      geom_abline(intercept = b_final,slope = slope_final,color = "skyblue",size = 1)+
      geom_point(x = xPt,y = y_cutoff,color = "skyblue",size = 3)+
      geom_point(x = xPt_final,y = y_cutoff_final,color = "skyblue",size = 3)+
      geom_text(x = 1+100,
                y = inputVector[1]+0.03,
                label = "A (start)",
                size = 8)+
      geom_text(x = length(inputVector)-100,
                y = inputVector[length(inputVector)]-300,
                label = "B (end)",
                size = 8)+
      geom_text(x = xPt,
                y = y_cutoff-0.03,
                label = "C",
                size = 8)+
      geom_text(x = xPt_final+500,
                y = y_cutoff_final-30,
                label = paste0("D (cutoff),","(","x=",xPt_final,",","y=",round(y_cutoff_final,digits = 4),")"),
                size = 8)+
      ggtitle(label = paste0("cutoff is ",signif(y_cutoff_final,digits = 3)))+
      ylab(tmp_ylab)+
      coord_cartesian(clip = "off")+
      theme_cowplot(font_size = 20)+
      theme(plot.title = element_text(hjust = 0.5))
    
    
  }
  
  return(list(absolute=y_cutoff_final,
              overMedian=y_cutoff_final/median(inputVector),
              overMean=y_cutoff_final/mean(inputVector),
              plot = p))
}


###-------get QC cell cutoff --------------
tmp.seu <- seu

tmp.data.plot <- tmp.seu[[]] %>%
  dplyr::arrange(-nFeature_SCP) %>%
  dplyr::mutate(rank = row_number()) %>%
  dplyr::select(nFeature_SCP,rank)

tmp.res.1 <- calculate_cutoff(inputVector = abs(tmp.data.plot$nFeature_SCP),
                              tmp_ylab = "Protein number")

tmp.cutoff <- tmp.res.1$absolute
p_cutoff <- tmp.res.1$plot+
  ggtitle("SCP")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))


###-------analysis QC results-------
# before cell QC
p <- VlnPlot(object = tmp.seu,
             features = "nFeature_SCP",
             cols = ggsci::pal_npg()(3)[2],
             group.by = "mass")+
  ggtitle("before cell QC")+
  ylab("Protein number")+
  xlab(NULL)+
  scale_x_discrete(label = c("SCP"))+
  theme_cowplot(font_size = 22)+
  NoLegend()+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))

tmp.levels <- c("2M","14M","24M")
tmp.seu$age <- factor(tmp.seu$age,levels = tmp.levels)
tmp.data.plot <- tmp.seu[[]] %>%
  group_by(brain_region,age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()
table(tmp.seu@meta.data$age,tmp.seu@meta.data$brain_region)
p <- ggplot(tmp.data.plot,
            aes(brain_region,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("before cell QC")+
  ylab("proportion (%)")+
  xlab(NULL)+
  scale_fill_manual(values = c("#CEEB9C","#5BB6A9","#5E4FA2"))+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))

# after cell QC
tmp.seu <- subset(tmp.seu,nFeature_SCP >= tmp.cutoff)
ncol(tmp.seu)
tmp.seu <- tmp.seu[rowSums(tmp.seu@assays$SCP@layers$counts) > 0, ] 
dim(GetAssayData(tmp.seu,layer ="counts"))

p <- VlnPlot(object = tmp.seu,
             features = "nFeature_SCP",
             cols = ggsci::pal_npg()(3)[2],
             group.by = "mass")+
  ylab("Protein number")+
  xlab(NULL)+
  ggtitle("after cell QC")+
  scale_x_discrete(label = c("SCP"))+
  theme_cowplot(font_size = 22)+
  NoLegend()+
  scale_y_continuous(limits = c(0,2800),
                     breaks = seq(0,2500,by=500))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))

tmp.data.plot <- tmp.seu[[]] %>%
  group_by(brain_region,age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()
p <- ggplot(tmp.data.plot,aes(brain_region,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("after cell QC")+
  ylab("proportion (%)")+
  xlab(NULL)+
  scale_fill_manual(values = c("#CEEB9C","#5BB6A9","#5E4FA2"))+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))

saveRDS(tmp.seu,
        file = "res/SCP_MG_3085_seurat_post_QC.rds")


