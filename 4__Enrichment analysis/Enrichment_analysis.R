library(devtools)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(RColorBrewer)
library(UpSetR)
require(DOSE)
library(dplyr)
library(Seurat)
library(data.table)
library(cowplot)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
require(remotes)
library(europepmc)
library(viridis)
set.seed(123)
rm(list=ls())
setwd('/.../4__Enrichment\ analysis/')



#---Load MSigDB get gene sets information and combine them to "msig_combined_gene_sets"----
#http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
msig__H= msigdbr::msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
msig__C2 = msigdbr::msigdbr(species = "Homo sapiens", category = "C2") #curated gene sets
msig__C5 = msigdbr::msigdbr(species = "Homo sapiens", category = "C5") #ontology gene sets

msig_CP_BioCarta = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:BIOCARTA")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:BIOCARTA")])
msig_CP_KEGG = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:KEGG")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:KEGG")])
msig_CP_PID = data.frame(msig__C2$gs_name[which(msig__C2$gs_subcat == "CP:PID")], msig__C2$gene_symbol[which(msig__C2$gs_subcat == "CP:PID")])
msig_CP_Reactome = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:REACTOME")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:REACTOME")])
msig_CP_WikiPathways = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")])
msig_GO_BP = data.frame(msig__C5$gs_name[msig__C5$gs_subcat %in% c("GO:BP")], msig__C5$gene_symbol[msig__C5$gs_subcat %in% c("GO:BP")])
msig_hallmark = data.frame(msig__H$gs_name, msig__H$gene_symbol)

#changing column names
colnames(msig_CP_BioCarta) <- c('pathway', 'genes')
colnames(msig_CP_KEGG) <- c("pathway", "genes")
colnames(msig_CP_PID) <- c("pathway", "genes")
colnames(msig_CP_Reactome) <- c("pathway", "genes")
colnames(msig_CP_WikiPathways) <- c("pathway", "genes")
colnames(msig_GO_BP) <- c("pathway", "genes")
colnames(msig_hallmark) <- c("pathway", "genes")

#combining columns
msig_combined_gene_sets <- rbind(msig_CP_BioCarta, msig_CP_KEGG, msig_CP_PID, msig_CP_Reactome, msig_CP_WikiPathways, msig_GO_BP, msig_hallmark)



#---Load DEGs----
human_hinze_2groups = read.csv("Hinze__human.csv", header=T)
human_lake_3groups = read.csv("Lake__human.csv", header=T)
human_abedini_2groups = read.csv("Abedini__human.csv")
mouse_balzer_7groups = read.csv("Balzer__mouse__human_genenames.csv", header=T)
mouse_doke_2groups = read.csv("Doke__mouse__human_genenames.csv", header=T)
mouse_kirita_6groups = read.csv("Kirita__mouse__human_genenames.csv", header=T)

mouse_balzer_7groups$gene <- NULL
mouse_doke_2groups$gene <- NULL
mouse_kirita_6groups$gene <- NULL
names(mouse_balzer_7groups)[names(mouse_balzer_7groups) == "HGNC.symbol"] <- "gene"
names(mouse_doke_2groups)[names(mouse_doke_2groups) == "HGNC.symbol"] <- "gene"
names(mouse_kirita_6groups)[names(mouse_kirita_6groups) == "HGNC.symbol"] <- "gene"
names(human_abedini_2groups)[names(human_abedini_2groups) == "X"] <- "gene"
names(mouse_balzer_7groups)[names(mouse_balzer_7groups) == "avg_logFC"] <- "avg_log2FC"



#---Prepare lists for clusterprofiler----
# mouse Balzer 7 groups
mouse_balzer7_Control <- subset(mouse_balzer_7groups, cluster=="Control")
mouse_balzer7_Control <- subset(mouse_balzer7_Control, avg_log2FC > 0)
mouse_balzer7_Control <- subset(mouse_balzer7_Control, p_val_adj < 0.05)
mouse_balzer7_Control <- mouse_balzer7_Control[order(mouse_balzer7_Control$avg_log2FC, decreasing=T),]
mouse_balzer7_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_Control
mouse_balzer7_Control <- mouse_balzer7_Control$gene

mouse_balzer7_AKI_1d <- subset(mouse_balzer_7groups, cluster=="IRI_short_1")
mouse_balzer7_AKI_1d <- subset(mouse_balzer7_AKI_1d, avg_log2FC > 0)
mouse_balzer7_AKI_1d <- subset(mouse_balzer7_AKI_1d, p_val_adj < 0.05)
mouse_balzer7_AKI_1d <- mouse_balzer7_AKI_1d[order(mouse_balzer7_AKI_1d$avg_log2FC, decreasing=T),]
mouse_balzer7_AKI_1d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_AKI_1d
mouse_balzer7_AKI_1d <- mouse_balzer7_AKI_1d$gene

mouse_balzer7_AKI_3d <- subset(mouse_balzer_7groups, cluster=="IRI_short_3")
mouse_balzer7_AKI_3d <- subset(mouse_balzer7_AKI_3d, avg_log2FC > 0)
mouse_balzer7_AKI_3d <- subset(mouse_balzer7_AKI_3d, p_val_adj < 0.05)
mouse_balzer7_AKI_3d <- mouse_balzer7_AKI_3d[order(mouse_balzer7_AKI_3d$avg_log2FC, decreasing=T),]
mouse_balzer7_AKI_3d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_AKI_3d
mouse_balzer7_AKI_3d <- mouse_balzer7_AKI_3d$gene

mouse_balzer7_AKI_14d <- subset(mouse_balzer_7groups, cluster=="IRI_short_14")
mouse_balzer7_AKI_14d <- subset(mouse_balzer7_AKI_14d, avg_log2FC > 0)
mouse_balzer7_AKI_14d <- subset(mouse_balzer7_AKI_14d, p_val_adj < 0.05)
mouse_balzer7_AKI_14d <- mouse_balzer7_AKI_14d[order(mouse_balzer7_AKI_14d$avg_log2FC, decreasing=T),]
mouse_balzer7_AKI_14d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_AKI_14d
mouse_balzer7_AKI_14d <- mouse_balzer7_AKI_14d$gene

mouse_balzer7_CKD_1d <- subset(mouse_balzer_7groups, cluster=="IRI_long_1")
mouse_balzer7_CKD_1d <- subset(mouse_balzer7_CKD_1d, avg_log2FC > 0)
mouse_balzer7_CKD_1d <- subset(mouse_balzer7_CKD_1d, p_val_adj < 0.05)
mouse_balzer7_CKD_1d <- mouse_balzer7_CKD_1d[order(mouse_balzer7_CKD_1d$avg_log2FC, decreasing=T),]
mouse_balzer7_CKD_1d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_CKD_1d
mouse_balzer7_CKD_1d <- mouse_balzer7_CKD_1d$gene

mouse_balzer7_CKD_3d <- subset(mouse_balzer_7groups, cluster=="IRI_long_3")
mouse_balzer7_CKD_3d <- subset(mouse_balzer7_CKD_3d, avg_log2FC > 0)
mouse_balzer7_CKD_3d <- subset(mouse_balzer7_CKD_3d, p_val_adj < 0.05)
mouse_balzer7_CKD_3d <- mouse_balzer7_CKD_3d[order(mouse_balzer7_CKD_3d$avg_log2FC, decreasing=T),]
mouse_balzer7_CKD_3d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_CKD_3d
mouse_balzer7_CKD_3d <- mouse_balzer7_CKD_3d$gene

mouse_balzer7_CKD_14d <- subset(mouse_balzer_7groups, cluster=="IRI_long_14")
mouse_balzer7_CKD_14d <- subset(mouse_balzer7_CKD_14d, avg_log2FC > 0)
mouse_balzer7_CKD_14d <- subset(mouse_balzer7_CKD_14d, p_val_adj < 0.05)
mouse_balzer7_CKD_14d <- mouse_balzer7_CKD_14d[order(mouse_balzer7_CKD_14d$avg_log2FC, decreasing=T),]
mouse_balzer7_CKD_14d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_balzer7_CKD_14d
mouse_balzer7_CKD_14d <- mouse_balzer7_CKD_14d$gene

mouse_balzer_7groups <- list(mouse_balzer7_Control=mouse_balzer7_Control,
                             mouse_balzer7_AKI_1d=mouse_balzer7_AKI_1d,
                             mouse_balzer7_AKI_3d=mouse_balzer7_AKI_3d,
                             mouse_balzer7_AKI_14d=mouse_balzer7_AKI_14d,
                             mouse_balzer7_CKD_1d=mouse_balzer7_CKD_1d,
                             mouse_balzer7_CKD_3d=mouse_balzer7_CKD_3d,
                             mouse_balzer7_CKD_14d=mouse_balzer7_CKD_14d)

# mouse Doke
mouse_doke_Control <- subset(mouse_doke_2groups, cluster=="Control")
mouse_doke_Control <- subset(mouse_doke_Control, avg_log2FC > 0)
mouse_doke_Control <- subset(mouse_doke_Control, p_val_adj < 0.05)
mouse_doke_Control <- mouse_doke_Control[order(mouse_doke_Control$avg_log2FC, decreasing=T),]
mouse_doke_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_doke_Control
mouse_doke_Control <- mouse_doke_Control$gene

mouse_doke_UUO <- subset(mouse_doke_2groups, cluster=="UUO")
mouse_doke_UUO <- subset(mouse_doke_UUO, avg_log2FC > 0)
mouse_doke_UUO <- subset(mouse_doke_UUO, p_val_adj < 0.05)
mouse_doke_UUO <- mouse_doke_UUO[order(mouse_doke_UUO$avg_log2FC, decreasing=T),]
mouse_doke_UUO %>% top_n(n = 100, wt = avg_log2FC) -> mouse_doke_UUO
mouse_doke_UUO <- mouse_doke_UUO$gene

mouse_doke_2groups <- list(mouse_doke_Control=mouse_doke_Control,
                           mouse_doke_UUO=mouse_doke_UUO)

# mouse Kirita
mouse_kirita_Control <- subset(mouse_kirita_6groups, cluster=="Control")
mouse_kirita_Control <- subset(mouse_kirita_Control, avg_log2FC > 0)
mouse_kirita_Control <- subset(mouse_kirita_Control, p_val_adj < 0.05)
mouse_kirita_Control <- mouse_kirita_Control[order(mouse_kirita_Control$avg_log2FC, decreasing=T),]
mouse_kirita_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_kirita_Control
mouse_kirita_Control <- mouse_kirita_Control$gene

mouse_kirita_4hours <- subset(mouse_kirita_6groups, cluster=="4hours")
mouse_kirita_4hours <- subset(mouse_kirita_4hours, avg_log2FC > 0)
mouse_kirita_4hours <- subset(mouse_kirita_4hours, p_val_adj < 0.05)
mouse_kirita_4hours <- mouse_kirita_4hours[order(mouse_kirita_4hours$avg_log2FC, decreasing=T),]
mouse_kirita_4hours %>% top_n(n = 100, wt = avg_log2FC) -> mouse_kirita_4hours
mouse_kirita_4hours <- mouse_kirita_4hours$gene

mouse_kirita_12hours <- subset(mouse_kirita_6groups, cluster=="12hours")
mouse_kirita_12hours <- subset(mouse_kirita_12hours, avg_log2FC > 0)
mouse_kirita_12hours <- subset(mouse_kirita_12hours, p_val_adj < 0.05)
mouse_kirita_12hours <- mouse_kirita_12hours[order(mouse_kirita_12hours$avg_log2FC, decreasing=T),]
mouse_kirita_12hours %>% top_n(n = 100, wt = avg_log2FC) -> mouse_kirita_12hours
mouse_kirita_12hours <- mouse_kirita_12hours$gene

mouse_kirita_6weeks <- subset(mouse_kirita_6groups, cluster=="6weeks")
mouse_kirita_6weeks <- subset(mouse_kirita_6weeks, avg_log2FC > 0)
mouse_kirita_6weeks <- mouse_kirita_6weeks[order(mouse_kirita_6weeks$avg_log2FC, decreasing=T),]
mouse_kirita_6weeks %>% top_n(n = 100, wt = avg_log2FC) -> mouse_kirita_6weeks
mouse_kirita_6weeks <- mouse_kirita_6weeks$gene

mouse_kirita_6groups <- list(mouse_kirita_Control=mouse_kirita_Control,
                             mouse_kirita_4hours=mouse_kirita_4hours,
                             mouse_kirita_12hours=mouse_kirita_12hours,
                             mouse_kirita_6weeks=mouse_kirita_6weeks)

# hinze 2 groups
human_hinze_Control <- subset(human_hinze_2groups, cluster=="Control")
human_hinze_Control <- subset(human_hinze_Control, avg_log2FC > 0)
human_hinze_Control <- subset(human_hinze_Control, p_val_adj < 0.05)
human_hinze_Control <- human_hinze_Control[order(human_hinze_Control$avg_log2FC, decreasing=T),]
human_hinze_Control %>% top_n(n = 100, wt = avg_log2FC) -> human_hinze_Control
human_hinze_Control <- human_hinze_Control$gene

human_hinze_AKI <- subset(human_hinze_2groups, cluster=="AKI")
human_hinze_AKI <- subset(human_hinze_AKI, avg_log2FC > 0)
human_hinze_AKI <- subset(human_hinze_AKI, p_val_adj < 0.05)
human_hinze_AKI <- human_hinze_AKI[order(human_hinze_AKI$avg_log2FC, decreasing=T),]
human_hinze_AKI %>% top_n(n = 100, wt = avg_log2FC) -> human_hinze_AKI
human_hinze_AKI <- human_hinze_AKI$gene

human_hinze_2groups <- list(human_hinze_Control=human_hinze_Control,
                            human_hinze_AKI=human_hinze_AKI)

# lake 3 groups
human_lake_LD <- subset(human_lake_3groups, cluster=="LD")
human_lake_LD <- subset(human_lake_LD, avg_log2FC > 0)
human_lake_LD <- subset(human_lake_LD, p_val_adj < 0.05)
human_lake_LD <- human_lake_LD[order(human_lake_LD$avg_log2FC, decreasing=T),]
human_lake_LD %>% top_n(n = 100, wt = avg_log2FC) -> human_lake_LD
human_lake_LD <- human_lake_LD$gene

human_lake_AKI <- subset(human_lake_3groups, cluster=="AKI")
human_lake_AKI <- subset(human_lake_AKI, avg_log2FC > 0)
human_lake_AKI <- subset(human_lake_AKI, p_val_adj < 0.05)
human_lake_AKI <- human_lake_AKI[order(human_lake_AKI$avg_log2FC, decreasing=T),]
human_lake_AKI %>% top_n(n = 100, wt = avg_log2FC) -> human_lake_AKI
human_lake_AKI <- human_lake_AKI$gene

human_lake_CKD <- subset(human_lake_3groups, cluster=="HCKD")
human_lake_CKD <- subset(human_lake_CKD, avg_log2FC > 0)
human_lake_CKD <- subset(human_lake_CKD, p_val_adj < 0.05)
human_lake_CKD <- human_lake_CKD[order(human_lake_CKD$avg_log2FC, decreasing=T),]
human_lake_CKD %>% top_n(n = 100, wt = avg_log2FC) -> human_lake_CKD
human_lake_CKD <- human_lake_CKD$gene

human_lake_3groups <- list(human_lake_LD=human_lake_LD,
                           human_lake_AKI=human_lake_AKI,
                           human_lake_CKD=human_lake_CKD)

#abedini 2 groups
human_abedini_Control<-subset(human_abedini_2groups, avg_log2FC>0)
human_abedini_Control<-subset(human_abedini_Control, p_val<0.05)
human_abedini_Control <- human_abedini_Control[order(human_abedini_Control$avg_log2FC, decreasing=T),]
human_abedini_Control %>% top_n(n = 100, wt = avg_log2FC) -> human_abedini_Control
human_abedini_Control <- human_abedini_Control$gene

human_abedini_CKD<-subset(human_abedini_2groups, avg_log2FC<0)
human_abedini_CKD<-subset(human_abedini_CKD, p_val<0.05)
human_abedini_CKD <- human_abedini_CKD[order(human_abedini_CKD$avg_log2FC, decreasing=F),]
human_abedini_CKD %>% top_n(n = -100, wt = avg_log2FC) -> human_abedini_CKD
human_abedini_CKD <- human_abedini_CKD$gene

human_abedini_2groups <- list(human_abedini_Control=human_abedini_Control,
                              human_abedini_CKD=human_abedini_CKD)



#---Prepare dataset list----
mice_x_human_reduced <- list(Abedini_Health=human_abedini_Control,
                             Balzer_Health=mouse_balzer7_Control,
                             Doke_Health=mouse_doke_Control,
                             Hinze_Health=human_hinze_Control,
                             Kirita_Health=mouse_kirita_Control,
                             Lake_Health=human_lake_LD,
                             Balzer_Recovery=mouse_balzer7_AKI_14d,
                             Balzer_AKI=mouse_balzer7_CKD_1d,
                             Hinze_AKI=human_hinze_AKI,
                             Kirita_AKI=mouse_kirita_12hours,
                             Lake_AKI=human_lake_AKI,
                             Abedini_CKD=human_abedini_CKD,
                             Balzer_CKD=mouse_balzer7_CKD_14d,
                             Doke_CKD=mouse_doke_UUO,                     
                             Lake_CKD=human_lake_CKD,
                             Kirita_CKD=mouse_kirita_6weeks)

#UpSetPlot for top100DEG lists
m = make_comb_mat(mice_x_human_reduced, mode = "distinct")
m <- m[comb_degree(m) <= 2]
tm = t(m)
tm[reverse(comb_name(tm)), reverse(set_name(tm))]
pdf('../Results/Top100DEGs_reduceddatasets_UpSetPlot_transposed_v4.pdf', width=6, height=10)
UpSet(tm, 
      set_order = names(mice_x_human_reduced),
      comb_order = order(-comb_size(tm)),
      pt_size = unit(5, "mm"), 
      lwd = 3,
      comb_col = c("black", "red")[comb_degree(tm)],
      top_annotation = HeatmapAnnotation(show_legend=F,
                                         "Set size" = anno_barplot(set_size(tm), 
                                                                   border = FALSE, 
                                                                   gp = gpar(fill = "black"), 
                                                                   width = unit(2, "cm")
                                         ),
                                         Phenotype = c(rep("Health", 6), rep("Recovery", 1), rep("AKI", 4), rep("CKD", 5)),
                                         col=list(Phenotype = c("Health"="#336600", "Recovery"="yellow","AKI"="red","CKD"="darkred")),
                                         annotation_name_side = "left", 
                                         annotation_name_rot = 0
      ),
      right_annotation = rowAnnotation(show_legend=F,
                                       Degree = as.character(comb_degree(tm)),
                                       col=list(Degree = c("1"="black", "2"="red")),
                                       "Intersection\nsize" = anno_barplot(comb_size(tm), 
                                                                           border = FALSE, 
                                                                           gp = gpar(fill = "black"), 
                                                                           height = unit(2, "cm"),
                                                                           add_numbers = T,
                                                                           side = "top"))
      
)
dev.off()



#---compare with Clusterprofiler----
cp_mice_x_human_reduced = clusterProfiler::compareCluster(geneCluster = mice_x_human_reduced, 
                                                          fun = "enricher", 
                                                          TERM2GENE = msig_combined_gene_sets)
cp_mice_x_human_reduced
#.. number of enriched terms found for each gene cluster:
#..   Abedini_Health: 223 
#..   Balzer_Health: 158 
#..   Doke_Health: 179 
#..   Hinze_Health: 33 
#..   Kirita_Health: 16 
#..   Lake_Health: 173 
#..   Balzer_Recovery: 213 
#..   Balzer_AKI: 514 
#..   Hinze_AKI: 298 
#..   Kirita_AKI: 146 
#..   Lake_AKI: 374 
#..   Abedini_CKD: 202 
#..   Balzer_CKD: 953 
#..   Doke_CKD: 39 
#..   Lake_CKD: 122 
#..   Kirita_CKD: 11 



#---Prepare dataset list----
mice_x_human <- list(Kirita_Health=mouse_kirita_Control,
                     Hinze_Health=human_hinze_Control,
                     Lake_Health=human_lake_LD,
                     Balzer_Health=mouse_balzer7_Control,
                     Doke_Health=mouse_doke_Control,
                     Abedini_Health=human_abedini_Control,
                     Balzer_Recovery=mouse_balzer7_AKI_14d,
                     Kirita_AKI_4hr=mouse_kirita_4hours,
                     Kirita_AKI_12hr=mouse_kirita_12hours,
                     Balzer_AKI_1d=mouse_balzer7_AKI_1d,
                     Balzer_AKI_3d=mouse_balzer7_AKI_3d,
                     Hinze_AKI=human_hinze_AKI,
                     Lake_AKI=human_lake_AKI,
                     Abedini_CKD=human_abedini_CKD,
                     Balzer_CKD=mouse_balzer7_CKD_14d,
                     Doke_CKD=mouse_doke_UUO,
                     Lake_CKD=human_lake_CKD)

#---compare with Clusterprofiler----
cp_mice_x_human = clusterProfiler::compareCluster(geneCluster = mice_x_human, 
                                                          fun = "enricher", 
                                                          TERM2GENE = msig_combined_gene_sets)
fwrite(cp_mice_x_human@compareClusterResult, '../Results/Enrichment_analysis.csv', row.names = T)
cp_mice_x_human
#.. number of enriched terms found for each gene cluster:
#..   Kirita_Health: 16 
#..   Hinze_Health: 33 
#..   Lake_Health: 173 
#..   Balzer_Health: 158 
#..   Doke_Health: 179 
#..   Abedini_Health: 223 
#..   Balzer_Recovery: 213 
#..   Kirita_AKI_4hr: 333 
#..   Kirita_AKI_12hr: 146 
#..   Balzer_AKI_1d: 248 
#..   Balzer_AKI_3d: 354 
#..   Hinze_AKI: 298 
#..   Lake_AKI: 374 
#..   Abedini_CKD: 202 
#..   Balzer_CKD: 953 
#..   Doke_CKD: 39 
#..   Lake_CKD: 122

#filter results by q-value
qvalue_cutoff = 0.05
cp_mice_x_human@compareClusterResult = cp_mice_x_human@compareClusterResult[which(cp_mice_x_human@compareClusterResult$qvalue < qvalue_cutoff),]
cp_mice_x_human_mat = xtabs(qvalue ~ ID + Cluster, data = cp_mice_x_human@compareClusterResult)
cp_mice_x_human_mat = -log10(cp_mice_x_human_mat)
cp_mice_x_human_mat[is.infinite(cp_mice_x_human_mat)] = 0

#top5 per cluster
test_mice_x_human <- data.frame(cp_mice_x_human_mat)
subset_mouse_kirita_Control = subset(test_mice_x_human, Cluster == "Kirita_Health")
subset_mouse_kirita_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_kirita_Control

subset_mouse_kirita_4hours = subset(test_mice_x_human, Cluster == "Kirita_AKI_4hr")
subset_mouse_kirita_4hours %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_kirita_4hours

subset_mouse_kirita_12hours = subset(test_mice_x_human, Cluster == "Kirita_AKI_12hr")
subset_mouse_kirita_12hours %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_kirita_12hours

subset_mouse_doke_Control = subset(test_mice_x_human, Cluster == "Doke_Health")
subset_mouse_doke_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_doke_Control

subset_mouse_doke_UUO = subset(test_mice_x_human, Cluster == "Doke_CKD")
subset_mouse_doke_UUO %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_doke_UUO

subset_mouse_balzer_Control = subset(test_mice_x_human, Cluster == "Balzer_Health")
subset_mouse_balzer_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_balzer_Control

subset_mouse_balzer_AKI_14d = subset(test_mice_x_human, Cluster == "Balzer_Recovery")
subset_mouse_balzer_AKI_14d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_balzer_AKI_14d

subset_mouse_balzer_AKI_1d = subset(test_mice_x_human, Cluster == "Balzer_AKI_1d")
subset_mouse_balzer_AKI_1d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_balzer_AKI_1d

subset_mouse_balzer_AKI_3d = subset(test_mice_x_human, Cluster == "Balzer_AKI_3d")
subset_mouse_balzer_AKI_3d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_balzer_AKI_3d

subset_mouse_balzer_CKD_14d = subset(test_mice_x_human, Cluster == "Balzer_CKD")
subset_mouse_balzer_CKD_14d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_balzer_CKD_14d

subset_human_lake_Control = subset(test_mice_x_human, Cluster == "Lake_Health")
subset_human_lake_Control %>% top_n(n = 5, wt = Freq) -> topSOI_human_lake_LD

subset_human_lake_AKI = subset(test_mice_x_human, Cluster == "Lake_AKI")
subset_human_lake_AKI %>% top_n(n = 5, wt = Freq) -> topSOI_human_lake_AKI

subset_human_lake_CKD = subset(test_mice_x_human, Cluster == "Lake_CKD")
subset_human_lake_CKD %>% top_n(n = 5, wt = Freq) -> topSOI_human_lake_CKD

subset_human_hinze_Control = subset(test_mice_x_human, Cluster == "Hinze_Health")
subset_human_hinze_Control %>% top_n(n = 5, wt = Freq) -> topSOI_human_hinze_Control

subset_human_hinze_AKI = subset(test_mice_x_human, Cluster == "Hinze_AKI")
subset_human_hinze_AKI %>%top_n(n = 5, wt = Freq) -> topSOI_human_hinze_AKI

subset_human_abedini_Control = subset(test_mice_x_human, cluster= "Abedini_Health")
subset_human_abedini_Control %>% top_n(n=5, wt=Freq) -> topSOI_human_abedini_Control

subset_human_abedini_CKD<-subset(test_mice_x_human, cluster="Abedini_CKD")
subset_human_abedini_CKD %>% top_n(n=5, wt=Freq) -> topSOI_human_abedini_CKD

topSOI_all <- rbind(topSOI_human_abedini_Control, 
                    topSOI_mouse_balzer_Control, 
                    topSOI_mouse_doke_Control, 
                    topSOI_human_hinze_Control, 
                    topSOI_mouse_kirita_Control,
                    topSOI_human_lake_LD,
                    topSOI_mouse_balzer_AKI_14d,
                    topSOI_mouse_balzer_AKI_1d,
                    topSOI_mouse_balzer_AKI_3d,
                    topSOI_human_hinze_AKI,
                    topSOI_mouse_kirita_4hours,
                    topSOI_mouse_kirita_12hours,
                    topSOI_human_lake_AKI, 
                    topSOI_human_abedini_CKD,
                    topSOI_mouse_balzer_CKD_14d, 
                    topSOI_mouse_doke_UUO,
                    topSOI_human_lake_CKD)

#subset cp_mice_x_human_mat to these top5 per cluster
cp_mice_x_human_mat <- cp_mice_x_human_mat[unique(topSOI_all$ID),]

#any q-value lower than 10e-5 is clamped to produce a more readable heatmap
cp_mice_x_human_mat[cp_mice_x_human_mat > 5] = 5



#---Heatmap of top5 enriched pathways----
cp_mice_x_human_mat2 <- cp_mice_x_human_mat[rownames(cp_mice_x_human_mat)[order(rownames(cp_mice_x_human_mat))],]
colnames(cp_mice_x_human_mat2)
#"Kirita_Health"   "Hinze_Health"    "Lake_Health"     "Balzer_Health"   "Doke_Health"     "Abedini_Health"  
#"Balzer_Recovery"
#"Kirita_AKI_4hr"  "Kirita_AKI_12hr" "Balzer_AKI_1d"   "Balzer_AKI_3d"   "Hinze_AKI"       "Lake_AKI"        
#"Abedini_CKD"  "Balzer_CKD"      "Doke_CKD"        "Lake_CKD"  
cluster_colors <- c("#336600","#336600", "#336600", "#336600","#336600", "#336600", 
                    "yellow",
                    "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", 
                    "#660000", "#660000", "#660000","#660000")
names(cluster_colors) <- c("Health", "Health", "Health", "Health","Health", "Health", 
                           "Recovery",  
                           "AKI", "AKI", "AKI","AKI", "AKI","AKI",
                           "CKD", "CKD", "CKD","CKD")
species_colors <- c("#F8C471","#85C1E9","#85C1E9", "#F8C471", "#F8C471","#85C1E9", "#F8C471", "#F8C471","#F8C471",
                    "#F8C471","#F8C471","#85C1E9", "#85C1E9", "#85C1E9", "#F8C471","#F8C471","#85C1E9")
names(species_colors) <- c("mouse","human","human", "mouse", "mouse","human", "mouse", "mouse","mouse",
                           "mouse","mouse","human", "human", "human", "mouse","mouse","human") #"#F8C471", #mouse; #"#85C1E9", #human
column_ha_cluster <- HeatmapAnnotation(
  Phenotype = names(cluster_colors),
  col = list(Phenotype = cluster_colors),
  simple_anno_size = unit(5, "mm"),
  annotation_name_gp = gpar(fontsize = 20),
  annotation_name_side = "left",
  show_annotation_name=F,
  show_legend = T,
  annotation_legend_param = list(
    Phenotype = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 20, 
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 20))))

column_ha_species <- HeatmapAnnotation(
  Species = names(species_colors),
  col = list(Species=species_colors),
  simple_anno_size = unit(5, "mm"),
  annotation_name_gp = gpar(fontsize = 20),
  annotation_name_side = "left",
  show_annotation_name=F,
  show_legend = T,
  annotation_legend_param = list(
    Species = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 20, 
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 20))))

# rename rows (pathways) with more readable names defined in "new_rownames.csv":
new_rownames<-read.csv("new_rownames.csv", header=T)
rownames_original<-rownames(cp_mice_x_human_mat2)
new_rownames<-new_rownames[match(rownames_original, new_rownames$ID),]
new_rownames$ID<-NULL   # short alternative to subset()
names(new_rownames)[names(new_rownames) == "new_rownames"] <- "ID"
rownames(cp_mice_x_human_mat2) <- new_rownames$ID

# remove duplicate pathways:
cp_mice_x_human_mat3<-cp_mice_x_human_mat2[-c(1,3,5,6,8,10,14,15,16,17,19,20,21,23,30,31,32,34,37,39,41,42,45,46,47,48,52),]

# plot
heatmap<-ComplexHeatmap::Heatmap(
  cp_mice_x_human_mat3, "-log10(q)",
  na_col="white",
  col= colorRampPalette(c("white", "#241346"))(n=200),
  cluster_rows=T,
  width= unit(30, "cm"), height=unit(30,"cm"),
  row_names_gp=gpar(fontsize=28),
  cluster_columns= F,
  top_annotation = column_ha_cluster,
  bottom_annotation= column_ha_species,
  show_heatmap_legend = T,
  column_names_gp=gpar(fontsize=28),
  column_names_rot = 45,
  use_raster = T,
  heatmap_legend_param = list(title = "-log10(q)",
                              title_gp = gpar(col = "black", fontsize = 20, fontface="bold"),
                              legend_height = unit(16, "cm"),
                              legend_width  = unit(16,"cm"),
                              title_position = "topleft",
                              direction = "vertical",
                              labels_gp = gpar(fontsize= 20, fontface = "bold")))
heatmap=draw(heatmap, heatmap_legend_side = "left", annotation_legend_side = "left")
pdf('../Results/Enrichment_heatmap.pdf', width=26, height=16)
heatmap
dev.off()
