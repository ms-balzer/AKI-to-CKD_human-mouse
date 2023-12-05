library(ggplot2)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
set.seed(123)
rm(list=ls())
setwd('/.../aki_to_ckd')




#============= compute jaccard similarity index ===================================
#setup jaccard function
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#example use case (iterate over all combinations of interest)
a <- Abedini_CKD_GOI
b <- Kirita_14days_GOI
jaccard(a,b)



#============= visualize heatmap ===================================
jaccard_mt <- read.csv('jaccard_index.csv', row.names=1)
jaccard_mt[jaccard_mt==1] <- -Inf
jaccard_mt[jaccard_mt==-Inf] <- 1

my_palette <- colorRampPalette(c("white", "red"))(n = 74)
my_palette <- c(my_palette, "#404040")
paletteLength = length(my_palette)
my_breaks <- c(seq(min(jaccard_mt), max_col, length.out=74),1)
species_colors <- c("#85C1E9", "#85C1E9", "#85C1E9", "#85C1E9", 
                    "#85C1E9", "#85C1E9", "#85C1E9",  #7 human
                    "#F8C471", "#F8C471", "#F8C471", "#F8C471",
                    "#F8C471", "#F8C471", "#F8C471", "#F8C471",
                    "#F8C471") #9 mice
names(species_colors) <- c("human", "human", "human", "human",
                           "human", "human", "human",  #7 human
                           "mouse", "mouse", "mouse", "mouse",
                           "mouse", "mouse", "mouse", "mouse",
                           "mouse") #9 mice
group_colors <- c("#336601",   "darkred",      "#336601",   "red" ,     "#336601" ,  "red"    ,  "darkred"     , "#336601" , "yellow", "red"  ,   
                  "darkred",      "#336601"  , "darkred"  ,    "#336601" ,  "red"  ,     "darkred" )
names(group_colors) <- sapply(X = strsplit(colnames(jaccard_mt), split = "_"), FUN = "[", 2)
column_ha_species <- HeatmapAnnotation(
  Species = names(species_colors),
  col = list(Species=species_colors),
  simple_anno_size = unit(2, "mm"),
  annotation_name_gp = gpar(fontsize = 11),
  annotation_name_side = "left",
  show_annotation_name=FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Species = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 11, 
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 11))))
column_ha_group <- HeatmapAnnotation(
  group = names(group_colors),
  col = list(group=group_colors),
  simple_anno_size = unit(2, "mm"),
  annotation_name_gp = gpar(fontsize = 11),
  annotation_name_side = "left",
  show_annotation_name=FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    group = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 11, 
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 11))))
my_heatmap <- pheatmap(as.matrix(jaccard_mt),
                       breaks=my_breaks,
                       cluster_rows = T,
                       cluster_cols = T,
                       scale="none",
                       heatmap_legend_param = list(title = "Jaccard index",
                                                   legend_height = unit(5, "cm")),
                       color=my_palette,
                       treeheight_row=10, treeheight_col=10, border_color="white",
                       show_rownames=T,
                       top_annotation= column_ha_species,
                       bottom_annotation= column_ha_group,
                       cellwidth=10, cellheight=10,
                       na_col="white")
pdf(file = paste0("../Results/DEGs_jaccard_pheatmap.pdf"), width=7.5, height=5.5)
my_heatmap
dev.off()
