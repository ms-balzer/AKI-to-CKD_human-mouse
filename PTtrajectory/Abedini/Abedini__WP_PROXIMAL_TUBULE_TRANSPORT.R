rm(list=ls())
library(slingshot)
library(Seurat)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)
library(RColorBrewer)
library(biomaRt)

library(dplyr)
library(data.table)
library(cowplot)
library(ggplot2)
library(monocle3)
library(htmlwidgets)
set.seed(123)

#load cds object
cds <- readRDS('Abedini__cds.PT.Innjured_PT1_2.def.8.10.2022.rds')


#---rename NA values in Groups.Def----
cds@colData@listData$Clustergroup <- plyr::mapvalues(
  x = cds@colData@listData$orig.ident,
  from = c("HK2558","HK2596","HK2663","HK2739","HK2770","HK2774","HK2833","HK2844","HK2862","HK2867","HK2868","HK2886","HK2891","HK2893","HK2895","HK2898","HK2899"),
  to = c("CKD", "CKD", "Control", "CKD", "CKD", "CKD", "Control", "CKD", "CKD", "Control", "CKD", "CKD", "CKD", "Control", "Control", "Control", "Control"))


cds@colData@listData$Clustergroup <- factor(cds@colData@listData$Clustergroup,
                                            levels = c("Control", "CKD"))



#================= Re-Cluster the cells =================
res=3e-3
cds = cluster_cells(cds, resolution=res)
#save clusters
colData(cds)$clustersNEW <- cds@clusters@listData$UMAP$clusters

cds@clusters@listData$UMAP$clusters <- plyr::mapvalues(
  x = cds@clusters@listData$UMAP$clusters,
  from = c("1", "2", "3", "4",
           "5", "6", "7","8",
           "9", "10", "11", "12",
           "13"),
  to = c("3","1","3","4", 
         "5","2","4","1",
         "7","7", "6", "2",
         "6"))

cds@clusters@listData$UMAP$clusters <- factor(cds@clusters@listData$UMAP$clusters,
                                              levels = c("1", "2", "3", "4", "5", "6", "7"))

colData(cds)$clustersNEW <- cds@clusters@listData$UMAP$clusters

#by clustersNEW
col_palette_short <- c("#4575b4", 
                                "#91bfdb",
                                "#e0f3f8",
                                "#ffffbf",
                                "#fee090",
                                "#fc8d59",
                                "#d73027")

pdf(file = paste0("plots_3_Dimensions_UMAP_by_7clustersNEW.pdf"), width=8, height=7)
plot_cells(cds, label_branch_points=F, label_cell_groups=F,label_roots =F,
           cell_size=1) + scale_color_manual(values=col_palette_short)
dev.off()



#================= Convert edited CDS to Seurat, preprocess and get diffusion map =================
#============= convert CDS to Seurat
count.mat <- assay(cds)
meta.df <- as.data.frame(colData(cds))

my.seurat <- CreateSeuratObject(counts = count.mat,
                                project = "my.project",
                                assay = "RNA",
                                meta.data = meta.df)



#======================== PATHWAY SCORING ===============================
#======================== get pathway of interest genes ===============================
pathwayofinterest <- "WP_PROXIMAL_TUBULE_TRANSPORT"

#===== Load POI and get GOI

#msig__H= msigdbr::msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
msig__C2 = msigdbr::msigdbr(species = "Homo sapiens", category = "C2") #curated gene sets
#msig__C5 = msigdbr::msigdbr(species = "Homo sapiens", category = "C5") #ontologygene sets

#msig_CP_BioCarta = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:BIOCARTA")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:BIOCARTA")])
#msig_CP_KEGG = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:KEGG")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:KEGG")])
#msig_CP_PID = data.frame(msig__C2$gs_name[which(msig__C2$gs_subcat == "CP:PID")], msig__C2$gene_symbol[which(msig__C2$gs_subcat == "CP:PID")])
#msig_CP_Reactome = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:REACTOME")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:REACTOME")])
msig_CP_WikiPathways = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")])
#msig_GO_BP = data.frame(msig__C5$gs_name[msig__C5$gs_subcat %in% c("GO:BP")], msig__C5$gene_symbol[msig__C5$gs_subcat %in% c("GO:BP")])
#msig_hallmark = data.frame(msig__H$gs_name, msig__H$gene_symbol)
colnames(msig_CP_WikiPathways) <- c("pathway", "gene")

POI = data.frame(msig_CP_WikiPathways$pathway[msig_CP_WikiPathways$pathway %in% c(pathwayofinterest)], msig_CP_WikiPathways$gene[msig_CP_WikiPathways$pathway %in% c(pathwayofinterest)])
colnames(POI) <- c("pathway", "gene")
GOI <- unique(POI$gene)





#======================== score this gene set  ===============================
my.seurat <- NormalizeData(my.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
my.seurat <- FindVariableFeatures(my.seurat, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(my.seurat)
my.seurat <- ScaleData(my.seurat, features = all.genes)
my.seurat <- RunPCA(my.seurat, features = VariableFeatures(object = my.seurat))
my.seurat <- RunUMAP(my.seurat, dims = 1:15)

#then, copy over monocle3 UMAP into Seurat
my.seurat@reductions$umap@cell.embeddings <- cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]]


dat <- my.seurat@assays$RNA@data #export normalized gene expression (not-scaled) gene by cell matrix
df <- as.data.frame(dat) #make df


#sort columns by Idents cell barcodes
df_sorting <- data.frame(barcode=names(my.seurat$Idents),
                         Idents=my.seurat$Idents) #export df with barcodes
df_sorting <- df_sorting[order(df_sorting$Idents),] #sort by Idents factor levels
df <- df[,df_sorting$barcode]


#score
reads_single_phase = df
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (GOI) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
GO_score = apply(reads_single_phase_restricted,2,mean))



#======================== visualize ===============================
#===== Violinplot

GO_score <- GO_score[order(factor(names(GO_score), 
                                  levels=colnames(my.seurat)))]
my.seurat$GO_score <- GO_score



#---for clustersNEW----

#Seurat style
#Idents(my.seurat) <- my.seurat$clustersNEW
#pdf(paste0(pathwayofinterest,'_score2_by_clustersNEW_VlnPlot.pdf'), width = 4, height=4)
#VlnPlot(my.seurat, features="GO_score", pt.size=0) +NoLegend()
#dev.off()



#plot with vioplot
library(vioplot)

makeTransparent = function(..., alpha=0.5) {
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent,alpha=alpha)
  return(newColor)
}


col_palette_short <- c("#4575b4",
                                "#91bfdb",
                                "#e0f3f8",
                                "#ffffbf",
                                "#fee090",
                                "#fc8d59",
                                "#d73027")

col_palette_short_trans = makeTransparent(col_palette_short, alpha = 0.3)


cluster1_barcodes <- WhichCells(my.seurat, idents="1")
cluster2_barcodes <- WhichCells(my.seurat, idents="2")
cluster3_barcodes <- WhichCells(my.seurat, idents="3")
cluster4_barcodes <- WhichCells(my.seurat, idents="4")
cluster5_barcodes <- WhichCells(my.seurat, idents="5")
cluster6_barcodes <- WhichCells(my.seurat, idents="6")
cluster7_barcodes <- WhichCells(my.seurat, idents="7")

dummy_order <- c(cluster1_barcodes, cluster2_barcodes, cluster3_barcodes, 
                 cluster4_barcodes, cluster5_barcodes, cluster6_barcodes, 
                 cluster7_barcodes
)

GO_score_vln <- GO_score[order(factor(names(GO_score), 
                                      levels=dummy_order))]



pdf(paste0(pathwayofinterest,'_score2_clustersNEW_',res,'.pdf'), width = 9, height=6)
facts = c(rep(1, table(my.seurat$clustersNEW)[1]),
          rep(2, table(my.seurat$clustersNEW)[2]),
          rep(3, table(my.seurat$clustersNEW)[3]),
          rep(4, table(my.seurat$clustersNEW)[4]), 
          rep(5, table(my.seurat$clustersNEW)[5]), 
          rep(6, table(my.seurat$clustersNEW)[6]),
          rep(7, table(my.seurat$clustersNEW)[7])
)
plot(GO_score_vln, 
  pch = NA,main = pathwayofinterest, 
  xaxt = "n", 
  xlab = "", 
  ylab = "Score",
  xlim = c(0.5,7.5),
  # ylim = c(0,0.9)
)
axis(side=1, #side 1 below, 2 left, 3 above, 4 right 
     las=1, #rotate perpendicular to axis
     at=1:length(table(my.seurat$clustersNEW)), names(table(my.seurat$clustersNEW)))
vioplot(GO_score_vln ~ facts, 
        plotCentre="line", col = col_palette_short_trans, border = col_palette_short, 
        lwd = 3, add = T
)
dev.off()


#---for Control + CKD (DKD + HKD)----

#Seurat style
#Idents(my.seurat) <- my.seurat$Clustergroup
#pdf(paste0(pathwayofinterest,'_score2_by_Clustergroup_VlnPlot.pdf'), width = 4, height=4)
#VlnPlot(my.seurat, features="GO_score", pt.size=0) +NoLegend()
#dev.off()

col_palette_short <- c("#91bfdb", "#b2182b")
col_palette_short_trans = makeTransparent(col_palette_short, alpha = 0.3)


# for Clustergroup
cluster1_barcodes <- WhichCells(my.seurat, idents="Control")
cluster2_barcodes <- WhichCells(my.seurat, idents="CKD")

dummy_order <- c(cluster1_barcodes, cluster2_barcodes)


GO_score_vln <- GO_score[order(factor(names(GO_score), 
                                      levels=dummy_order))]


pdf(paste0(pathwayofinterest,'_score2_Clustergroup_',res,'.pdf'), width = 3.5, height=4)
facts = c(rep(1, table(my.seurat$Clustergroup)[1]),
          rep(2, table(my.seurat$Clustergroup)[2]))
plot(GO_score_vln, 
  pch = NA,
  main = pathwayofinterest, 
  xaxt = "n", 
  xlab = "", 
  ylab = "Score",
  xlim = c(0.5,2.5))
axis(side=1, #side 1 below, 2 left, 3 above, 4 right 
     las=2, #rotate perpendicular to axis
     at=1:length(table(my.seurat$Clustergroup)), names(table(my.seurat$Clustergroup)))
vioplot(GO_score_vln ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, 
        lwd = 3, add = T
)
dev.off()





#======================== visualize GO score in UMAP space  ===============================
#plot pathway enrichment in umap space
rd <- my.seurat@reductions$umap@cell.embeddings
colnames(rd) <- c("UMAP_1", "UMAP_2")
colnames(my.seurat@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")


library(viridis)
nc <- 1
GOscore_df <- data.frame(my.seurat$GO_score)
GOscore_mt <- matrix(my.seurat$GO_score)
rownames(GOscore_mt) <- rownames(GOscore_df)
colnames(GOscore_mt) <- paste0(pathwayofinterest," Score")
nms <- paste0(pathwayofinterest," Score")
nr <- ceiling(length(nms)/nc)
pal <- rev(inferno(100, end = 1))
pdf(paste0(pathwayofinterest,'_GOScore_inferno.pdf'), width=nc*4.5, height=5)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(GOscore_mt[,i], breaks = 100)]
  plot(rd, col = colors, pch = 16, cex = 0.5, main = i)
}
dev.off()



#Seurat featureplot with GO_score
#library(viridis)
#pdf(file = paste0(pathwayofinterest,"_FeaturePlot_umap.pdf"), width=5, height=5)
#FeaturePlot(my.seurat,
#           features = "GO_score", 
#          combine = FALSE,
#         cols=rev(inferno(100)),
#        pt.size=1)
#dev.off()





