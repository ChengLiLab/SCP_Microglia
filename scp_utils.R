require(ComplexHeatmap)
require(RColorBrewer)
require(ggplot2)
###-------IO-------
myFileName <- function(prefix,suffix){
  res <- paste0(prefix,"_",format(Sys.time(),"%Y%m%d"),suffix)
  return(res)
}

myggsave <- function(p,prefix,suffix,width=8,height=8,...){
  if(suffix==".pdf"){
    ggsave(p,filename = myFileName(prefix = prefix,suffix = suffix),width = width,height = height,useDingbats=F,...)
  }
  else{
    ggsave(p,filename = myFileName(prefix = prefix,suffix = suffix),width = width,height = height,...)
  }
}
myggsavePro <- function(p,prefix,width=8,height=8,dpi=350){
  myggsave(p = p,
           prefix = prefix,
           suffix = ".jpg",
           width = width,
           height = height,
           dpi = dpi)
  myggsave(p = p,
           prefix = prefix,
           suffix = ".pdf",
           width = width,
           height = height)
}

MyHeatmapSave <- function(ph,prefix,width,height,legend.position="right"){
    pdf(file = myFileName(prefix = prefix,suffix = ".pdf"),
        width = width,height = height)
    draw(ph,
         heatmap_legend_side = legend.position)
    dev.off()
    
    jpeg(file = myFileName(prefix = prefix,suffix = ".jpg"),
         width = width,height = height,units = "in",res = 350)
    draw(ph,heatmap_legend_side = legend.position)
    dev.off()
  
}



###------ Utils----------
get_pheatmap_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

myRemoveNA <- function(df){
  df[is.na(df)] <- 0
  return(df)
}

myRemoveNA <- function(df){
  df[is.na(df)] <- 0
  return(df)
}
median_center <- function(x){
  res <- x - apply(x, 2, median)
  return(res)
}
myMatrixFilter <- function(data.exp,filter=0.95){
  data.exp <- apply(data.exp,2,FUN = function(x) ifelse(x >= quantile(x,filter),quantile(x,filter),x))
  return(data.exp)
}

MyHeatmapSave <- function(ph,prefix,width,height){
  pdf(file = myFileName(prefix = prefix,suffix = ".pdf"),
      width = width,height = height)
  draw(ph)
  dev.off()
  
  jpeg(file = myFileName(prefix = prefix,suffix = ".jpg"),
       width = width,height = height,units = "in",res = 350)
  draw(ph)
  dev.off()
}






###------ define colors ----------------

rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
skyblueYelow <- colorRampPalette(c("skyblue","black","yellow"))
skybluered <- colorRampPalette(c("skyblue","black","orange"))
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
blues <- colorRampPalette(colors = brewer.pal(9,"Blues"))
ylord <- colorRampPalette(colors = brewer.pal(9,"YlOrRd"))
hic.red <- colorRampPalette(c("white","red"))
# "#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913" "#D94801" "#A63603"
# "#7F2704"
hic.orage <- colorRampPalette(c("#FFF5EB","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
sc.hic.orange <- colorRampPalette(c("grey85","#EAD397","#FDB31A","#E42A2A","#A31D1D"))

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)

onlywarm <-colorRampPalette(colors = c("#f7fcf0", "#E5D8D1", "#F2CBB7",
                                       "#F7B89C", "#F6A081", "#EE8568",
                                       "#E0654F", "#CC4039", "#B40426"))
coolwarm <- colorRampPalette(colors = c("#3B4CC0", "#4F6AD9", "#6585EC", "#7B9FF9", "#93B5FF", 
                                        "#AAC7FD", "#C0D4F5", "#D4DBE6", "#E5D8D1", "#F2CBB7",
                                        "#F7B89C", "#F6A081", "#EE8568",
                                        "#E0654F", "#CC4039", "#B40426"))
ramp <- colorRampPalette(c("white","pink","red","black"))

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
hic.pca.red <- colorRampPalette(c("blue","gray1","red"))
hic.pca.redwhite <- colorRampPalette(c("#1d1856","navyblue","white","red4","#861617"))
hic.pca.orange <- colorRampPalette(c("#2f2583","black","#f9b232"))
hic.pca.skyblue <- colorRampPalette(c("skyblue","black","orange"))

okabe_ito <- colorRampPalette(colors = c("#e59f01","#56b4e8","#009f73","#f0e442","#0072b1","#d55e00","#cc79a7","#999999","#000000"))
OrBl_div <- colorRampPalette(colors = c("#9f3d22","#be4d21","#db6525","#ef8531","#f1ac73","#d8d4c9","#a1bccf","#6fa3cb","#5689b6","#4171a1","#2b5b8b"))
OrBl_v2 <- colorRampPalette(colors = rev(c("#940025","#e50c2f","#F15F30" ,"#F7962E", "#FCEE2B",
                                       "#88CEEF","#248AF3","#14B3FF","#3361A5","#004377")))

light_molande <- colorRampPalette(colors = c("#413C58","#A3C4BC","#BFD7B5","#E7EFC5","#F2DDA4"))


coldwarm2.0 <- colorRampPalette(c("#3C1518","#69140E","#A44200","#D58936","#F2F3AE",
                                  "#D4E4BC","#96ACB7","#36558F","#40376E","#48233C"))

###gundam
strike_freedom_gundam_color <- colorRampPalette(colors = c("brown3","#eff3ff","navy","#373b35","#b1b0c2","#f2be58"))
npg.color <- colorRampPalette(colors = ggsci:::ggsci_db$npg$nrc) 

my_stepped3.color <- colorRampPalette(colors = pals::stepped3(20))
my_stepped.color <- colorRampPalette(colors = pals::stepped(24))
#mycolor.bar(hic.pca.red(100),min = -1,max = 1)
divergentcolor <- function (n) {
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#F29403", "#F781BF", "#BC9DCC", 
                  "#A65628", "#54B0E4", "#222F75", "#1B9E77", 
                  "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A", 
                  "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", 
                  "#B3DE69", "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

####some more divergent color
####from monet
monet.sunset <- colorRampPalette(colors = c("#272924","#323c48","#67879c",
                                            "#7d9390","#d7695a","#ba9f84"))
monet.sunumbrella.women <- colorRampPalette(colors = c("#aed1e4","#bcdde2","#e8e7d5",
                                                       "#f1ede9","#e6cec2","#e49d40",
                                                       "#ddcb3e","#a39524","#89af7a",
                                                       "#744054","#aea3a9"))
monet.cliff <- colorRampPalette(colors = c("#242a26","#3c5653","#bacbc3","#3d4b6e",
                                           "#6f799d","#b3bedc","#f0f0f2","#898788",
                                           "#ddd0bd"))
####show pictures
pic.green <- colorRampPalette(colors = c("#344A23","#726342","#383812","#879979"))
pic.blue <- colorRampPalette(colors = c("#0B2C26","#004054","#95B7B6","#DFE5E1","#C26441"))
pic.green2 <- colorRampPalette(colors = c("#61776B","#A6AA8B","#B6C7D7","#EFEBEC","#E79A61"))
pic.beauty <- colorRampPalette(colors = c("#2F5365","#8A8E8F","#435D2E","#09200B","#CDBBA4","#D04A49"))
pic.beauty2 <- colorRampPalette(colors = c("#4B6032","#436A4D","#DFB28D","#BF4E31","#87988E"))
pic.beauty3 <- colorRampPalette(colors = c("#DB9371","#FCE6D8","#EB9B63","#DBD6C9","#666E20"))
pic.greenblue <- colorRampPalette(colors = c("#86AFA9","#E9A419","#FC8550","#F7C9BA","#8C5E51"))
pic.boy <- colorRampPalette(colors = c("#7C7F62","#956A3D","#625A58","#B1B1B1","#D8D4D5"))
pic.fellows <- colorRampPalette(colors = c("#193B14","#85A063","#8A4925","#C5947E","#ACA9A2"))
pic.girl <- colorRampPalette(colors = c("#1B3E64","#BE1705","#D77363","#84865D","#E2B88D"))
pic.girl2 <- colorRampPalette(colors = c("#BF593F","#DAAE7D","#AAC5DA","#DAC1BA","#141A22"))
molandi.color <- c("#e4d4c5","#c6b1ac","#764e56",
                   "#f9d9c0","#d19477","#93675a",
                   "#f0e9cd","#b9a783","#796656",
                   "#cdc1b1","#a2967e","#656356",
                   "#d8e7e4","#9eb2b1","#5a6873",
                   "#ccd8b0","#7e8563","#50463d")
molandi.color <- colorRampPalette(colors = molandi.color)
#scales::show_col(pic.girl2(10))
pinkblack<- colorRampPalette(colors = c("#333333","#666A86","#95B8D1","#E8DDB5","#EDAFB8"))

beaches <- colorRampPalette(colors = c("#87D2DB" ,"#5BB1CB", "#4F66AF" ,"#F15F30" ,"#F7962E", "#FCEE2B" ))

molandi.color3 <- colorRampPalette(colors = c("#FCD0A1","#B1B695","#A690A4","#5E4B56","#AFD2E9"))

sun_gradual <- colorRampPalette(colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"))


rdbu_dp <- colorRampPalette(colors = c("#004377","#007ebe","#00a7d1","#81d1e6","#d0ebf4",
                                       "#f9f9f9","#ffe1cf","#ffb18f","#ff6e58","#e50c2f","#940025"))

# tmp.color <- colorRampPalette(colors = c("#612C6D","#D46889","#B2B2B2",
#                                          "#869CC3","#895C48"))
# 
# tmp.color <- c("#73808E","#86BDD8","#fae79b","#F46F43",
#                "#edb386","#344F99","#895C48","#DCD7C1")

# tmp.color <- c("#51c4c2","#0d8a8c","#4583b3","#f78e26","#f172ad","#f7afb9",
#                "#c63596","#be86ba","#8b66b8","#4068b2","#512a93","#223271")

science.color <- colorRampPalette(c("#b03d26","#005f81",
                                    "#9ccfe6","#e0897e",
                                    "#a5a7ab","black"))



YlGnBu_gradient <- colorRampPalette(rev(hcl.colors(100,palette = "YlGnBu")))

YlGnBu_gradientv2 <- colorRampPalette(colors = brewer.pal(n = 9,name = "YlGnBu"))


BlYl_gradient <- colorRampPalette(colors = rev(hcl.colors(100,palette = "Blue-Yellow")))

# tmp.color <- c('#000000', '#380000', '#560000', '#760100', '#980300', '#bb0600', '#df0d00', '#f93500', '#fe6800', '#ff9100', '#ffb402', '#ffd407', '#fff324')
# myfire.color <- colorRampPalette(tmp.color)

tmp.color <- c('#000000', '#380000', '#560000', '#760100', '#980300', '#bb0600', '#df0d00', '#f93500', '#fe6800', '#ff9100', '#ffb402', '#ffd407','#fff324','#ffee9d','#fdf7d2')
myfire.color <- colorRampPalette(tmp.color[c(1,2,7,11,12,13,14)])

###-------plot function----------------
pheatmap_fixed <- function (mat, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                         name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, 
                            border_color_colum_ann = "black",
                            border_color_row_ann = NA,
                            border_color = ifelse(nrow(mat) < 100 & ncol(mat) < 100, 
                                                  "grey60", NA), cellwidth = NA, cellheight = NA, 
                            scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
                            clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
                            clustering_method = "complete", clustering_callback = NA, 
                            cutree_rows = NA, cutree_cols = NA, treeheight_row = ifelse(class(cluster_rows) == 
                                                                                          "hclust" || cluster_rows, 50, 0), treeheight_col = ifelse(class(cluster_cols) == 
                                                                                                                                                      "hclust" || cluster_cols, 50, 0), legend = TRUE, 
                            legend_breaks = NA, legend_labels = NA, annotation_row = NA, 
                            annotation_col = NA, annotation = NA, annotation_colors = NA, 
                            annotation_legend = TRUE, annotation_names_row = TRUE, annotation_names_col = TRUE, 
                            drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, 
                            main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, 
                            angle_col = c("270", "0", "45", "90", 
                                          "315"), display_numbers = FALSE, number_format = "%.2f", 
                            number_color = "grey30", fontsize_number = 0.8 * fontsize, 
                            gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL, 
                            filename = NA, width = NA, height = NA, silent = FALSE, na_col = "#DDDDDD", 
                            name = NULL, heatmap_legend_param = list(), ..., run_draw = FALSE) 
{
  if (is.data.frame(mat)) {
    ComplexHeatmap:::warning_wrap("The input is a data frame, convert it to the matrix.")
    mat = as.matrix(mat)
  }
  if (!identical(kmeans_k, NA)) {
    ComplexHeatmap:::warning_wrap("argument `kmeans_k` is not suggested to use in pheatmap -> Heatmap translation because it changes the input matrix. You might check `row_km` and `column_km` arguments in Heatmap().")
    km = kmeans(mat, centers = kmeans_k)
    mat = km$centers
    rownames(mat) = paste0("Cluster: ", seq_along(km$size), 
                           ", Size: ", km$size)
  }
  if ("row" %in% scale) {
    if (any(is.na(mat))) {
      mat = (mat - rowMeans(mat, na.rm = TRUE))/rowSds(mat, 
                                                       na.rm = TRUE)
    }
    else {
      mat = t(scale(t(mat)))
    }
  }
  else if ("column" %in% scale) {
    if (any(is.na(mat))) {
      mat = t((t(mat) - colMeans(mat, na.rm = TRUE))/colSds(mat, 
                                                            na.rm = TRUE))
    }
    else {
      mat = scale(mat)
    }
  }
  ht_param = list(matrix = mat)
  if (!identical(scale, "none") && !identical(breaks, 
                                              NA)) {
    ComplexHeatmap:::warning_wrap("It not suggested to both set `scale` and `breaks`. It makes the function confused.")
  }
  if (is.function(color)) {
    ht_param$col = color
    if (!identical(breaks, NA)) {
      ComplexHeatmap:::warning_wrap("`breaks` is ignored when `color` is set as a color mapping function.")
    }
  }
  else {
    if (identical(breaks, NA)) {
      n_col = length(color)
      if (identical(scale, "row") || identical(scale, 
                                               "column")) {
        lim = max(abs(mat), na.rm = TRUE)
        ht_param$col = circlize::colorRamp2(seq(-lim, lim, length = n_col), 
                                            color)
      }
      else {
        ht_param$col = circlize::colorRamp2(seq(min(mat, na.rm = TRUE), 
                                                max(mat, na.rm = TRUE), length = n_col), color)
      }
    }
    else {
      if (length(breaks) == length(color) + 1) {
        ht_param$col = local({
          breaks = breaks
          color = color
          fun = function(x) {
            n = length(color)
            df = data.frame(start = c(-Inf, breaks[seq_len(n)], 
                                      breaks[n + 1]), end = c(breaks[1], breaks[1 + 
                                                                                  seq_len(n)], Inf))
            ind = numeric(length(x))
            for (i in seq_along(x)) {
              ind[i] = which(df$start <= x[i] & df$end > 
                               x[i])
            }
            ind = ind - 1
            ind[ind < 1] = 1
            ind[ind > n] = n
            color[ind]
          }
          attr(fun, "breaks") = breaks
          fun
        })
      }
      else if (length(breaks) == length(color)) {
        ht_param$col = circlize::colorRamp2(breaks, color)
      }
      else {
        n_col = length(color)
        ht_param$col = circlize::colorRamp2(seq(min(breaks), max(breaks), 
                                                length = n_col), color)
        ComplexHeatmap:::warning_wrap("`breaks` does not have the same length as `color`. The colors are interpolated from the minimal to the maximal of `breaks`.")
      }
    }
  }
  if (!identical(filename, NA)) {
    ComplexHeatmap:::warning_wrap("argument `filename` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(width, NA)) {
    ComplexHeatmap:::warning_wrap("argument `width` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(height, NA)) {
    ComplexHeatmap:::warning_wrap("argument `height` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(silent, FALSE)) {
    ComplexHeatmap:::warning_wrap("argument `silent` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  ht_param$rect_gp = gpar(col = border_color)
  if (nrow(mat) > 1000 || ncol(mat) > 1000) {
    if (!is.na(border_color)) {
      ComplexHeatmap:::warning_wrap("border color is set for the matrix with large numbers of rows or columns. You might only be able to see the border colors in the plot. Set `border_color = NA` to get rid of it.")
    }
  }
  if (!identical(cellwidth, NA)) {
    ht_param$width = ncol(mat) * unit(cellwidth, "pt")
  }
  if (!identical(cellheight, NA)) {
    ht_param$height = nrow(mat) * unit(cellheight, "pt")
  }
  if (identical(clustering_distance_rows, "correlation")) 
    clustering_distance_rows = "pearson"
  if (identical(clustering_distance_cols, "correlation")) 
    clustering_distance_cols = "pearson"
  ht_param$cluster_rows = cluster_rows
  ht_param$cluster_columns = cluster_cols
  ht_param$clustering_distance_rows = clustering_distance_rows
  ht_param$clustering_distance_columns = clustering_distance_cols
  ht_param$clustering_method_rows = clustering_method
  ht_param$clustering_method_columns = clustering_method
  if (!is.na(cutree_rows)) {
    if (inherits(cluster_rows, c("logical", "hclust", 
                                 "dendrogram"))) {
      ht_param$row_split = cutree_rows
      ht_param$row_gap = unit(4, "bigpts")
      ht_param["row_title"] = list(NULL)
    }
  }
  if (!is.na(cutree_cols)) {
    if (inherits(cluster_cols, c("logical", "hclust", 
                                 "dendrogram"))) {
      ht_param$column_split = cutree_cols
      ht_param$column_gap = unit(4, "bigpts")
      ht_param["column_title"] = list(NULL)
    }
  }
  ht_param$row_dend_width = unit(treeheight_row, "pt")
  ht_param$column_dend_height = unit(treeheight_col, "pt")
  ht_param$show_heatmap_legend = legend
  if (identical(scale, "row") || identical(scale, "column")) {
    if (identical(legend_breaks, NA)) {
      lim = quantile(abs(mat), 0.975)
      le = pretty(c(-lim, lim), n = 3)
      if (length(le) == 7 && le[1] == -3) {
        le = c(-3, -1.5, 0, 1.5, 3)
      }
      else if (!0 %in% le) {
        le = c(le[1], le[1]/2, 0, le[length(le)]/2, le[length(le)])
      }
      legend_breaks = le
    }
  }
  if (!identical(legend_breaks, NA)) {
    heatmap_legend_param$at = legend_breaks
  }
  if (!identical(legend_labels, NA)) {
    heatmap_legend_param$labels = legend_labels
  }
  ht_param$heatmap_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                       labels_gp = gpar(fontsize = fontsize*0.8))
  if (identical(annotation_colors, NA)) {
    annotation_colors = list()
  }
  if (!identical(annotation_col, NA)) {
    acn = rownames(annotation_col)
    mcn = colnames(mat)
    if (!is.null(acn)) {
      if (acn[1] %in% mcn) {
        if (length(union(acn, mcn)) == length(mcn)) {
          if (!identical(acn, mcn)) {
            ComplexHeatmap:::warning_wrap("Column annotation has different order from matrix columns. Adjust the column annotation based on column names of the matrix.")
          }
          annotation_col = annotation_col[mcn, , drop = FALSE]
        }
      }
    }
    for (nm in colnames(annotation_col)) {
      if (nm %in% names(annotation_colors)) {
        if (is.null(names(annotation_colors[[nm]])) && 
            is.numeric(annotation_col[, nm])) {
          foo_x = annotation_col[, nm]
          foo_n_col = length(annotation_colors[[nm]])
          annotation_colors[[nm]] = circlize::colorRamp2(seq(min(foo_x), 
                                                             max(foo_x), length = foo_n_col), annotation_colors[[nm]])
        }
      }
    }
    ht_param$top_annotation = HeatmapAnnotation(df = annotation_col[, 
                                                                    ncol(annotation_col):1, drop = FALSE], col = annotation_colors, 
                                                show_legend = annotation_legend, show_annotation_name = annotation_names_col, 
                                                annotation_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                                                               labels_gp = gpar(fontsize = fontsize*0.8)),
                                                gp = gpar(col = border_color_colum_ann), annotation_name_gp = gpar(fontsize = fontsize, 
                                                                                                                   fontface = "bold"), simple_anno_size = unit(10, 
                                                                                                                                                               "bigpts"), gap = unit(2, "bigpts"))
  }
  if (!identical(annotation_row, NA)) {
    arn = rownames(annotation_row)
    mrn = rownames(mat)
    if (!is.null(arn)) {
      if (arn[1] %in% mrn) {
        if (length(union(arn, mrn)) == length(mrn)) {
          if (!identical(arn, mrn)) {
            ComplexHeatmap:::warning_wrap("Row annotation has different order from matrix rows. Adjust the row annotation based on row names of the matrix.")
          }
          annotation_row = annotation_row[mrn, , drop = FALSE]
        }
      }
    }
    for (nm in colnames(annotation_row)) {
      if (nm %in% names(annotation_colors)) {
        if (is.null(names(annotation_colors[[nm]])) && 
            is.numeric(annotation_row[, nm])) {
          foo_x = annotation_row[, nm]
          foo_n_col = length(annotation_colors[[nm]])
          annotation_colors[[nm]] = circlize::colorRamp2(seq(min(foo_x), 
                                                             max(foo_x), length = foo_n_col), annotation_colors[[nm]])
        }
      }
    }
    ht_param$left_annotation = rowAnnotation(df = annotation_row[, 
                                                                 ncol(annotation_row):1, drop = FALSE], col = annotation_colors, 
                                             show_legend = annotation_legend, show_annotation_name = annotation_names_row, 
                                             annotation_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                                                            labels_gp = gpar(fontsize = fontsize*0.8)),
                                             gp = gpar(col = border_color_row_ann), annotation_name_gp = gpar(fontsize = fontsize, 
                                                                                                              fontface = "bold"), simple_anno_size = unit(10, 
                                                                                                                                                          "bigpts"), gap = unit(2, "bigpts"))
  }
  if (!identical(annotation, NA)) {
    ComplexHeatmap:::warning_wrap("argument `annotation` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (identical(drop_levels, FALSE)) {
    ComplexHeatmap:::warning_wrap("argument `drop_levels` is enfored to be TRUE, skip it.")
  }
  ht_param$show_row_names = show_rownames
  ht_param$show_column_names = show_colnames
  ht_param$row_names_gp = gpar(fontsize = fontsize_row)
  ht_param$column_names_gp = gpar(fontsize = fontsize_col)
  angle_col = match.arg(angle_col)[1]
  angle_col = switch(angle_col, `0` = 0, `45` = 45, 
                     `90` = 90, `270` = 90, `315` = -45)
  ht_param$column_names_rot = angle_col
  if (angle_col == 0) {
    ht_param$column_names_centered = TRUE
  }
  if (is.logical(display_numbers)) {
    if (display_numbers) {
      ht_param$layer_fun = local({
        number_format = number_format
        number_color = number_color
        fontsize_number = fontsize_number
        mat = mat
        function(j, i, x, y, w, h, fill) {
          grid.text(sprintf(number_format, pindex(mat, 
                                                  i, j)), x = x, y = y, gp = gpar(col = number_color, 
                                                                                  fontsize = fontsize_number))
        }
      })
    }
  }
  else if (is.matrix(display_numbers)) {
    if (!identical(dim(display_numbers), dim(mat))) {
      stop_wrap("dimension of `display_numbers` should be the same as the input matrix.")
    }
    ht_param$layer_fun = local({
      number_color = number_color
      fontsize_number = fontsize_number
      mat = display_numbers
      function(j, i, x, y, w, h, fill) {
        grid.text(pindex(mat, i, j), x = x, y = y, gp = gpar(col = number_color, 
                                                             fontsize = fontsize_number))
      }
    })
  }
  if (!is.null(labels_row)) {
    ht_param$row_labels = labels_row
  }
  if (!is.null(labels_col)) {
    ht_param$column_labels = labels_col
  }
  if (!is.null(gaps_row)) {
    if (inherits(cluster_rows, c("hclust", "dendrogram"))) {
      stop_wrap("`gaps_row` should not be set when `cluster_rows` is set as a clustering object.")
    }
    if (identical(cluster_rows, TRUE)) {
      stop_wrap("`gaps_row` should not be set when `cluster_rows` is set to TRUE.")
    }
    slices = diff(c(0, gaps_row, nrow(mat)))
    ht_param$row_split = rep(seq_along(slices), times = slices)
    ht_param$row_gap = unit(4, "bigpts")
    ht_param["row_title"] = list(NULL)
  }
  if (!is.null(gaps_col)) {
    if (inherits(cluster_cols, c("hclust", "dendrogram"))) {
      stop_wrap("`gaps_col` should not be set when `cluster_cols` is set as a clustering object.")
    }
    if (identical(cluster_cols, TRUE)) {
      stop_wrap("`gaps_col` should not be set when `cluster_cols` is set to TRUE.")
    }
    slices = diff(c(0, gaps_col, ncol(mat)))
    ht_param$column_split = rep(seq_along(slices), times = slices)
    ht_param$column_gap = unit(4, "bigpts")
    ht_param["column_title"] = list(NULL)
  }
  if (!identical(clustering_callback, NA)) {
    if (!identical(ht_param$cluster_rows, FALSE)) {
      row_hclust = hclust(get_dist(mat, ht_param$clustering_distance_rows), 
                          ht_param$clustering_method_rows)
      row_hclust = clustering_callback(row_hclust, ...)
      ht_param$cluster_rows = row_hclust
    }
    if (!identical(ht_param$cluster_columns, FALSE)) {
      column_hclust = hclust(get_dist(t(mat), ht_param$clustering_distance_columns), 
                             ht_param$clustering_method_columns)
      column_hclust = clustering_callback(column_hclust, 
                                          ...)
      ht_param$cluster_columns = column_hclust
    }
  }
  ht_param$name = name
  ht_param$row_dend_reorder = FALSE
  ht_param$column_dend_reorder = FALSE
  if (!identical(main, NA)) {
    ht_param$column_title = main
    ht_param$column_title_gp = gpar(fontface = "bold", 
                                    fontsize = 1.3 * fontsize)
  }
  ht_param = c(ht_param, list(...))
  ht = do.call(Heatmap, ht_param)
  attr(ht, "translate_from") = "pheatmap"
  if (run_draw) {
    draw(ht)
  }
  else {
    ht
  }
  #return(ht_param)
}

draw_key_text

draw_number_circle <- function (data, params, size) {
  grobTree(pointsGrob(x = 0.5, 
                      y = 0.5, 
                      size = unit(1.5, "char"), 
                      pch = 16, 
                      gp = gpar(col = alpha(data$colour %||% "grey50", 
                                                      data$alpha), 
                                fill = alpha(data$fill %||% "grey50", 
                                             data$alpha), lwd = (data$linewidth %||% 0.5) * .pt, 
                                          lty = data$linetype %||% 1)), 
           textGrob(label = data$label, 
                    x = rep(0.5, 3), 
                    y = rep(0.5, 3), 
                    gp = gpar(col = "black")))
}



### my DimPlot


myGetSeuLayout <- function(object,embed_use = "umap"){
  reduc <- data.frame(Seurat::Embeddings(object, 
                                         reduction = embed_use))
  meta <- object@meta.data
  res_layout <- cbind(reduc, meta)
  return(res_layout)
}

mySeuDimPlot <- function(object,
                         tmp_group,
                         embed_use = "umap",
                         embed_label =  "umap",
                         theme_size = 18,
                         tmp.alpha = 1){
  res_layout <- myGetSeuLayout(object = object,
                               embed_use = embed_use)
  tmp_nlevel <- length(unique(res_layout[[tmp_group]]))
  tmp_number <- 1:tmp_nlevel
  
  umap_1 <- paste0(embed_label,"_1")
  umap_2 <- paste0(embed_label,"_2")
  p <- ggplot(data = res_layout,
              aes(!!sym(umap_1),!!sym(umap_2), 
                  color = !!sym(tmp_group)))+
    geom_point(key_glyph = draw_number_circle,alpha = tmp.alpha)+
    guides(color = guide_legend(override.aes = list(label = tmp_number)))+
    theme_bw(base_size = theme_size,
             base_line_size = 1,
             base_rect_size = 2)+
    theme(panel.grid = element_blank())
  return(p)
}




