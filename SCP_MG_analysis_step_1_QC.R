#####################################################################################
## Project: SCP_Microglia                                                          ##
## Script Purpose: Quality Control                                                 ##
## Data: 2025.03.20                                                                ##
## Author: Wang Longteng, Zhang Haoran                                             ##
#####################################################################################
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

##-------1. SCP data Quality Control ----------
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
  # slope <- (max(inputVector)-inputVector)/length(inputVector)
  # tmp.idx <- which(abs(inputVector - (max(inputVector)-0.1)) < 0.0001)
  # inputVector[length(inputVector) - tmp.idx]
  # inputVector[1]
  # slope <- max(inputVector)
  
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
      #geom_hline(yintercept = y_cutoff_final,color = "skyblue",size = 1,linetype = "dashed")+
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
      # scale_x_continuous(limits = c(0,700),expand = c(0,0))+
      # scale_y_continuous(limits = c(0,1),expand = c(0,0))+
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

###-------overview data-------------
tmp.path <- "res/r/SCP_MG_3411_seurat_rmhis_ig_try10_alb_amy_krt_zhr_20241228.rds"
seu <- readRDS(file = tmp.path)



seu$mass <- plyr::mapvalues(seu$mass,
                            from = c("Debris"),
                            to = "Blank")
seu$mass <- factor(seu$mass,levels = c("Blank","1c","20c"))
head(seu[[]])

table(seu$mass)


p <- VlnPlot_scCustom(seurat_object = seu,
                      features = 'nFeature_SCP',
                      group.by = "mass")+
  scale_fill_npg()+
  xlab(NULL)+
  ggtitle(NULL)+
  ylab("Protein number")+
  NoLegend()+
  scale_x_discrete(label = c("Blank","SCP","20c"))
p
myggsavePro(p = p,
            prefix = "res/fig/SCP_MG_SCP_stat_rmhis_ig_try10_alb_Amy_krt",
            width = 3,
            height = 4)



tmp.cells <- WhichCells(object = seu,
                        expression = (mass == "1c" & 
                                        age == "2M"))
tmp.seu <- seu[,tmp.cells]

table(tmp.seu$brain_region,tmp.seu$prepare_date)

tmp.cells <- WhichCells(object = seu,
                        expression = (mass == "1c" & 
                                        age == "14M"))
tmp.seu <- seu[,tmp.cells]
table(tmp.seu$brain_region,tmp.seu$prepare_date)


tmp.cells <- WhichCells(object = seu,
                        expression = (mass == "1c" & 
                                        age == "24M"))
tmp.seu <- seu[,tmp.cells]
table(tmp.seu$brain_region,tmp.seu$prepare_date)

###-------get QC cell cutoff --------------
tmp.cells <- WhichCells(object = seu,
                        expression = (mass == "1c" & 
                                        celltype == "MG"))
tmp.seu <- seu[,tmp.cells]

#ncol(seu)
ncol(tmp.seu)

tmp.data.plot <- tmp.seu[[]] %>%
  dplyr::arrange(-nFeature_SCP) %>%
  dplyr::mutate(rank = row_number()) %>%
  dplyr::select(nFeature_SCP,rank)

tmp.res.1 <- calculate_cutoff(inputVector = abs(tmp.data.plot$nFeature_SCP),
                              tmp_ylab = "Protein number")

tmp.cutoff <- tmp.res.1$absolute
p <- tmp.res.1$plot+
  ggtitle("SCP")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
p
#p
myggsavePro(p = p,
            prefix = "res/fig/20241228SCP_MG_get_SCP_428cellcutof_QC_allcell",
            width = 8,
            height = 8,
            dpi = 350)


###-------analysis QC results-------
### Vlnplot
tmp.cells <- WhichCells(object = seu,
                        expression = (mass == "1c"& 
                          celltype == "MG"))
tmp.seu <- seu[,tmp.cells]
ncol(tmp.seu)
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
p
myggsavePro(p = p,
            prefix = paste0("res/fig/SCP_MG_",
                            "SCP_","before_QC_",
                            "rmhis_ig_try10_alb_Amy_krt",
                            "protein_number_vlnplot"),
            width = 6,
            height = 6,
            dpi = 350)
tmp.levels <- c("2M","14M","24M")
tmp.seu$age <- factor(tmp.seu$age,levels = tmp.levels)
tmp.data.plot <- tmp.seu[[]] %>%
  group_by(brain_region,age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

p <- ggplot(tmp.data.plot,
            aes(brain_region,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("before cell QC")+
  ylab("proportion (100%)")+
  xlab(NULL)+
  #scale_x_discrete(label = c("Hippocampus","Prefrontal cortex"))+
  scale_fill_manual(values = col.spectral(6)[4:6])+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))
# theme(plot.title = element_text(hjust = 0.5),
#       axis.text.x = element_text(angle = 45,
#                                  hjust = 1,
#                                  vjust = 1))
p
myggsavePro(p = p,
            prefix = paste0("res/fig/SCP_MG_",
                            "SCP_","before_QC_",
                            "protein_number_stackplot"),
            width = 6,
            height = 6,
            dpi = 350)

tmp.seu <- subset(tmp.seu,nFeature_SCP >= tmp.cutoff)
#tmp.seu <- subset(tmp.seu,nFeature_SCP >= 767)
ncol(tmp.seu)

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
p
myggsavePro(p = p,
            prefix = paste0("res/fig/SCP_MG_",
                            "SCP_",
                            "after_QC_428cutoff","protein_number_vlnplot"),
            width = 6,
            height = 6,
            dpi = 350)


tmp.data.plot <- tmp.seu[[]] %>%
  group_by(brain_region,age) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup()

p <- ggplot(tmp.data.plot,aes(brain_region,pct,fill=age))+
  geom_bar(stat = "identity",color = "black")+
  ggtitle("after cell QC")+
  ylab("proportion (100%)")+
  xlab(NULL)+
  #scale_x_discrete(label = c("Hippocampus","Prefrontal cortex"))+
  scale_fill_manual(values = col.spectral(6)[4:6])+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1,by=0.25),
                     labels = seq(0,100,by=25))+
  theme_cowplot(font_size = 22)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))
p
myggsavePro(p = p,
            prefix = paste0("res/fig/SCP_MG_",
                            "SCP_","after_QC_428cutoff",
                            "protein_number_stackplot"),
            width = 6,
            height = 6,
            dpi = 350)
ncol(tmp.seu)

saveRDS(tmp.seu,
        file = "res/r/SCP_MG_3186_seurat_post_QC_428cutoff_zhr_20241228.rds")