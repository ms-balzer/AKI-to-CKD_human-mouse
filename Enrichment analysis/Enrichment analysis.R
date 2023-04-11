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

set.seed(123)
# setwd()


#---MSigDB get gene sets information and combine them to "msig_combined"----
#http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

# msigdbr::msigdbr_species() #show species
# print(msigdbr::msigdbr_collections(), n=23) #show all rows of msigdbr::msigdbr_collections()

msig__H= msigdbr::msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
msig__C2 = msigdbr::msigdbr(species = "Homo sapiens", category = "C2") #curated gene sets
msig__C5 = msigdbr::msigdbr(species = "Homo sapiens", category = "C5") #ontologygene sets

msig_CP_BioCarta = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:BIOCARTA")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:BIOCARTA")])
msig_CP_KEGG = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:KEGG")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:KEGG")])
msig_CP_PID = data.frame(msig__C2$gs_name[which(msig__C2$gs_subcat == "CP:PID")], msig__C2$gene_symbol[which(msig__C2$gs_subcat == "CP:PID")])
msig_CP_Reactome = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:REACTOME")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:REACTOME")])
msig_CP_WikiPathways = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")])
msig_GO_BP = data.frame(msig__C5$gs_name[msig__C5$gs_subcat %in% c("GO:BP")], msig__C5$gene_symbol[msig__C5$gs_subcat %in% c("GO:BP")])
msig_hallmark = data.frame(msig__H$gs_name, msig__H$gene_symbol)

#changing columnames
colnames(msig_CP_BioCarta) <- c('pathway', 'genes')
colnames(msig_CP_KEGG) <- c("pathway", "genes")
colnames(msig_CP_PID) <- c("pathway", "genes")
colnames(msig_CP_Reactome) <- c("pathway", "genes")
colnames(msig_CP_WikiPathways) <- c("pathway", "genes")
colnames(msig_GO_BP) <- c("pathway", "genes")
colnames(msig_hallmark) <- c("pathway", "genes")

#combining columns
msig_combined <- rbind(msig_CP_BioCarta, msig_CP_KEGG, msig_CP_PID, msig_CP_Reactome, msig_CP_WikiPathways, msig_GO_BP, msig_hallmark)


# reading .csv-files
human_Hinze = read.csv("Hinze__human.csv", header=T)
human_Lake = read.csv("Lake__human.csv", header=T)
human_Abedini = read.csv("Abedini__human.csv")

mouse_Balzer = read.csv("Balzer__mouse__human_genenames.csv", header=T)
mouse_Doke = read.csv("Doke__mouse__human_genenames.csv", header=T)
mouse_Kirita = read.csv("Kirita__mouse__human_genenames.csv", header=T)


#remove column mouse "gene"
#mouse_Balzer_3groups$gene <- NULL   # short alternative to subset()
mouse_Balzer$gene <- NULL   # short alternative to subset()
mouse_Doke$gene <- NULL   # short alternative to subset()
mouse_Kirita$gene <- NULL   # short alternative to subset()

# making sure columnnames are matching
#names(mouse_Balzer_3groups)[names(mouse_Balzer_3groups) == "HGNC.symbol"] <- "gene"
names(mouse_Balzer)[names(mouse_Balzer) == "HGNC.symbol"] <- "gene"
names(mouse_Doke)[names(mouse_Doke) == "HGNC.symbol"] <- "gene"
names(mouse_Kirita)[names(mouse_Kirita) == "HGNC.symbol"] <- "gene"
names(human_Abedini)[names(human_Abedini) == "X"] <- "gene"

#names(mouse_Balzer_3groups)[names(mouse_Balzer_3groups) == "avg_logFC"] <- "avg_log2FC"
names(mouse_Balzer)[names(mouse_Balzer) == "avg_logFC"] <- "avg_log2FC"


#---prepare lists for Clusterprofiler----

# define cluster

# mouse Balzer 7 groups
mouse_Balzer_Control <- subset(mouse_Balzer, cluster=="Control")
mouse_Balzer_Control <- subset(mouse_Balzer_Control, avg_log2FC > 0)
mouse_Balzer_Control <- subset(mouse_Balzer_Control, p_val_adj < 0.05)
mouse_Balzer_Control <- mouse_Balzer_Control[order(mouse_Balzer_Control$avg_log2FC, decreasing=T),]
mouse_Balzer_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_Control
mouse_Balzer_Control <- mouse_Balzer_Control$gene

mouse_Balzer_AKI_1d <- subset(mouse_Balzer, cluster=="IRI_short_1")
mouse_Balzer_AKI_1d <- subset(mouse_Balzer_AKI_1d, avg_log2FC > 0)
mouse_Balzer_AKI_1d <- subset(mouse_Balzer_AKI_1d, p_val_adj < 0.05)
mouse_Balzer_AKI_1d <- mouse_Balzer_AKI_1d[order(mouse_Balzer_AKI_1d$avg_log2FC, decreasing=T),]
mouse_Balzer_AKI_1d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_AKI_1d
mouse_Balzer_AKI_1d <- mouse_Balzer_AKI_1d$gene

mouse_Balzer_AKI_3d <- subset(mouse_Balzer, cluster=="IRI_short_3")
mouse_Balzer_AKI_3d <- subset(mouse_Balzer_AKI_3d, avg_log2FC > 0)
mouse_Balzer_AKI_3d <- subset(mouse_Balzer_AKI_3d, p_val_adj < 0.05)
mouse_Balzer_AKI_3d <- mouse_Balzer_AKI_3d[order(mouse_Balzer_AKI_3d$avg_log2FC, decreasing=T),]
mouse_Balzer_AKI_3d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_AKI_3d
mouse_Balzer_AKI_3d <- mouse_Balzer_AKI_3d$gene

mouse_Balzer_AKI_14d <- subset(mouse_Balzer, cluster=="IRI_short_14")
mouse_Balzer_AKI_14d <- subset(mouse_Balzer_AKI_14d, avg_log2FC > 0)
mouse_Balzer_AKI_14d <- subset(mouse_Balzer_AKI_14d, p_val_adj < 0.05)
mouse_Balzer_AKI_14d <- mouse_Balzer_AKI_14d[order(mouse_Balzer_AKI_14d$avg_log2FC, decreasing=T),]
mouse_Balzer_AKI_14d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_AKI_14d
mouse_Balzer_AKI_14d <- mouse_Balzer_AKI_14d$gene

mouse_Balzer_CKD_1d <- subset(mouse_Balzer, cluster=="IRI_long_1")
mouse_Balzer_CKD_1d <- subset(mouse_Balzer_CKD_1d, avg_log2FC > 0)
mouse_Balzer_CKD_1d <- subset(mouse_Balzer_CKD_1d, p_val_adj < 0.05)
mouse_Balzer_CKD_1d <- mouse_Balzer_CKD_1d[order(mouse_Balzer_CKD_1d$avg_log2FC, decreasing=T),]
mouse_Balzer_CKD_1d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_CKD_1d
mouse_Balzer_CKD_1d <- mouse_Balzer_CKD_1d$gene

mouse_Balzer_CKD_3d <- subset(mouse_Balzer, cluster=="IRI_long_3")
mouse_Balzer_CKD_3d <- subset(mouse_Balzer_CKD_3d, avg_log2FC > 0)
mouse_Balzer_CKD_3d <- subset(mouse_Balzer_CKD_3d, p_val_adj < 0.05)
mouse_Balzer_CKD_3d <- mouse_Balzer_CKD_3d[order(mouse_Balzer_CKD_3d$avg_log2FC, decreasing=T),]
mouse_Balzer_CKD_3d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_CKD_3d
mouse_Balzer_CKD_3d <- mouse_Balzer_CKD_3d$gene

mouse_Balzer_CKD_14d <- subset(mouse_Balzer, cluster=="IRI_long_14")
mouse_Balzer_CKD_14d <- subset(mouse_Balzer_CKD_14d, avg_log2FC > 0)
mouse_Balzer_CKD_14d <- subset(mouse_Balzer_CKD_14d, p_val_adj < 0.05)
mouse_Balzer_CKD_14d <- mouse_Balzer_CKD_14d[order(mouse_Balzer_CKD_14d$avg_log2FC, decreasing=T),]
mouse_Balzer_CKD_14d %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Balzer_CKD_14d
mouse_Balzer_CKD_14d <- mouse_Balzer_CKD_14d$gene


mouse_Balzer <- list(mouse_Balzer_Control=mouse_Balzer_Control,
                     mouse_Balzer_AKI_1d=mouse_Balzer_AKI_1d,
                     mouse_Balzer_AKI_3d=mouse_Balzer_AKI_3d,
                     mouse_Balzer_AKI_14d=mouse_Balzer_AKI_14d,
                     mouse_Balzer_CKD_1d=mouse_Balzer_CKD_1d,
                     mouse_Balzer_CKD_3d=mouse_Balzer_CKD_3d,
                     mouse_Balzer_CKD_14d=mouse_Balzer_CKD_14d)


# mouse Doke
mouse_Doke_Control <- subset(mouse_Doke, cluster=="Control")
mouse_Doke_Control <- subset(mouse_Doke_Control, avg_log2FC > 0)
mouse_Doke_Control <- subset(mouse_Doke_Control, p_val_adj < 0.05)
mouse_Doke_Control <- mouse_Doke_Control[order(mouse_Doke_Control$avg_log2FC, decreasing=T),]
mouse_Doke_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Doke_Control
mouse_Doke_Control <- mouse_Doke_Control$gene

mouse_Doke_UUO <- subset(mouse_Doke, cluster=="UUO")
mouse_Doke_UUO <- subset(mouse_Doke_UUO, avg_log2FC > 0)
mouse_Doke_UUO <- subset(mouse_Doke_UUO, p_val_adj < 0.05)
mouse_Doke_UUO <- mouse_Doke_UUO[order(mouse_Doke_UUO$avg_log2FC, decreasing=T),]
mouse_Doke_UUO %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Doke_UUO
mouse_Doke_UUO <- mouse_Doke_UUO$gene

mouse_Doke <- list(mouse_Doke_Control=mouse_Doke_Control,
                   mouse_Doke_UUO=mouse_Doke_UUO)


# mouse Kirita
mouse_Kirita_Control <- subset(mouse_Kirita, cluster=="Control")
mouse_Kirita_Control <- subset(mouse_Kirita_Control, avg_log2FC > 0)
mouse_Kirita_Control <- subset(mouse_Kirita_Control, p_val_adj < 0.05)
mouse_Kirita_Control <- mouse_Kirita_Control[order(mouse_Kirita_Control$avg_log2FC, decreasing=T),]
mouse_Kirita_Control %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Kirita_Control
mouse_Kirita_Control <- mouse_Kirita_Control$gene

mouse_Kirita_4hours <- subset(mouse_Kirita, cluster=="4hours")
mouse_Kirita_4hours <- subset(mouse_Kirita_4hours, avg_log2FC > 0)
mouse_Kirita_4hours <- subset(mouse_Kirita_4hours, p_val_adj < 0.05)
mouse_Kirita_4hours <- mouse_Kirita_4hours[order(mouse_Kirita_4hours$avg_log2FC, decreasing=T),]
mouse_Kirita_4hours %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Kirita_4hours
mouse_Kirita_4hours <- mouse_Kirita_4hours$gene

mouse_Kirita_12hours <- subset(mouse_Kirita, cluster=="12hours")
mouse_Kirita_12hours <- subset(mouse_Kirita_12hours, avg_log2FC > 0)
mouse_Kirita_12hours <- subset(mouse_Kirita_12hours, p_val_adj < 0.05)
mouse_Kirita_12hours <- mouse_Kirita_12hours[order(mouse_Kirita_12hours$avg_log2FC, decreasing=T),]
mouse_Kirita_12hours %>% top_n(n = 100, wt = avg_log2FC) -> mouse_Kirita_12hours
mouse_Kirita_12hours <- mouse_Kirita_12hours$gene

# for Kirita_2d, _14d and _6weeks were only 12 (1,0,11) enriched terms found, they were therefor excluded
mouse_Kirita <- list(mouse_Kirita_Control=mouse_Kirita_Control,
                     mouse_Kirita_4hours=mouse_Kirita_4hours,
                     mouse_Kirita_12hours=mouse_Kirita_12hours)


# Hinze 2 groups
human_Hinze_Control <- subset(human_Hinze, cluster=="Control")
human_Hinze_Control <- subset(human_Hinze_Control, avg_log2FC > 0)
human_Hinze_Control <- subset(human_Hinze_Control, p_val_adj < 0.05)
human_Hinze_Control <- human_Hinze_Control[order(human_Hinze_Control$avg_log2FC, decreasing=T),]
human_Hinze_Control %>% top_n(n = 100, wt = avg_log2FC) -> human_Hinze_Control
human_Hinze_Control <- human_Hinze_Control$gene

human_Hinze_AKI <- subset(human_Hinze, cluster=="AKI")
human_Hinze_AKI <- subset(human_Hinze_AKI, avg_log2FC > 0)
human_Hinze_AKI <- subset(human_Hinze_AKI, p_val_adj < 0.05)
human_Hinze_AKI <- human_Hinze_AKI[order(human_Hinze_AKI$avg_log2FC, decreasing=T),]
human_Hinze_AKI %>% top_n(n = 100, wt = avg_log2FC) -> human_Hinze_AKI
human_Hinze_AKI <- human_Hinze_AKI$gene

human_Hinze <- list(human_Hinze_Control=human_Hinze_Control,
                    human_Hinze_AKI=human_Hinze_AKI)


# Lake 3 groups
human_Lake_LD <- subset(human_Lake, cluster=="LD")
human_Lake_LD <- subset(human_Lake_LD, avg_log2FC > 0)
human_Lake_LD <- subset(human_Lake_LD, p_val_adj < 0.05)
human_Lake_LD <- human_Lake_LD[order(human_Lake_LD$avg_log2FC, decreasing=T),]
human_Lake_LD %>% top_n(n = 100, wt = avg_log2FC) -> human_Lake_LD
human_Lake_LD <- human_Lake_LD$gene

human_Lake_AKI <- subset(human_Lake, cluster=="AKI")
human_Lake_AKI <- subset(human_Lake_AKI, avg_log2FC > 0)
human_Lake_AKI <- subset(human_Lake_AKI, p_val_adj < 0.05)
human_Lake_AKI <- human_Lake_AKI[order(human_Lake_AKI$avg_log2FC, decreasing=T),]
human_Lake_AKI %>% top_n(n = 100, wt = avg_log2FC) -> human_Lake_AKI
human_Lake_AKI <- human_Lake_AKI$gene

human_Lake_CKD <- subset(human_Lake, cluster=="HCKD")
human_Lake_CKD <- subset(human_Lake_CKD, avg_log2FC > 0)
human_Lake_CKD <- subset(human_Lake_CKD, p_val_adj < 0.05)
human_Lake_CKD <- human_Lake_CKD[order(human_Lake_CKD$avg_log2FC, decreasing=T),]
human_Lake_CKD %>% top_n(n = 100, wt = avg_log2FC) -> human_Lake_CKD
human_Lake_CKD <- human_Lake_CKD$gene

human_Lake <- list(human_Lake_LD=human_Lake_LD,
                   human_Lake_AKI=human_Lake_AKI,
                   human_Lake_CKD=human_Lake_CKD)


#Abedini 2 groups
human_Abedini_Control<-subset(human_Abedini, avg_log2FC>0)
human_Abedini_Control<-subset(human_Abedini_Control, p_val<0.05)
human_Abedini_Control <- human_Abedini_Control[order(human_Abedini_Control$avg_log2FC, decreasing=T),]
human_Abedini_Control %>% top_n(n = 100, wt = avg_log2FC) -> human_Abedini_Control
human_Abedini_Control <- human_Abedini_Control$gene

human_Abedini_CKD<-subset(human_Abedini, avg_log2FC<0)
human_Abedini_CKD<-subset(human_Abedini_CKD, p_val<0.05)
human_Abedini_CKD <- human_Abedini_CKD[order(human_Abedini_CKD$avg_log2FC, decreasing=F),]
human_Abedini_CKD %>% top_n(n = -100, wt = avg_log2FC) -> human_Abedini_CKD
human_Abedini_CKD <- human_Abedini_CKD$gene

human_Abedini <- list(human_Abedini_Control=human_Abedini_Control,
                      human_Abedini_CKD=human_Abedini_CKD)




# prepare list for Clusterprofiler
mice_x_human <- list(m_Kirita_Control=mouse_Kirita_Control,
                     h_Hinze_Control=human_Hinze_Control,
                     h_Lake_LD=human_Lake_LD,
                     m_Balzer_Control=mouse_Balzer_Control,
                     m_Doke_Control=mouse_Doke_Control,
                     h_Abedini_Control=human_Abedini_Control,
                     m_Balzer_AKI_14d=mouse_Balzer_AKI_14d, #7 Control
                     m_Kirita_4hours=mouse_Kirita_4hours,
                     m_Kirita_12hours=mouse_Kirita_12hours,
                     m_Balzer_AKI_1d=mouse_Balzer_AKI_1d,
                     h_Lake_AKI=human_Lake_AKI,
                     h_Hinze_AKI=human_Hinze_AKI,
                     m_Balzer_AKI_3d=mouse_Balzer_AKI_3d, #6 AKI
                     h_Abedini_CKD=human_Abedini_CKD,
                     m_Balzer_CKD_3d=mouse_Balzer_CKD_3d,
                     h_Lake_CKD=human_Lake_CKD,
                     m_Doke_UUO=mouse_Doke_UUO,
                     m_Balzer_CKD_1d=mouse_Balzer_CKD_1d,
                     m_Balzer_CKD_14d=mouse_Balzer_CKD_14d) #6 CKD


#---compare DEGs with UpSetR-----
upset_mice_x_human<-UpSetR::fromList(mice_x_human)

pdf("upsetplot_mice_x_human_DEGs.pdf", height=25, width=20)
UpSetR::upset(upset_mice_x_human, nsets=19, shade.color = 329, color.pal = 589,
              line.size=1.4, point.size=3, text.scale = 3,
              keep.order=T,
              sets=c("m_Balzer_CKD_14d", "m_Balzer_CKD_3d", "m_Balzer_CKD_1d",
                     "m_Doke_UUO", "h_Lake_CKD", "h_Abedini_CKD", # 6 CKD
                     #"m_Kirita_6weeks", "m_Kirita_14d", "m_Kirita_2d",
                     "m_Balzer_AKI_3d", "m_Balzer_AKI_1d",
                     "m_Kirita_4hours", "m_Kirita_12hours", 
                     "h_Lake_AKI", "h_Hinze_AKI", # 6 AKI
                     "m_Balzer_AKI_14d",
                     "h_Lake_LD","h_Hinze_Control", "h_Abedini_Control",
                     "m_Balzer_Control", "m_Kirita_Control", "m_Doke_Control"# 7 Control
              ),
              queries=
                list(list(query=intersects,
                          params=list("m_Kirita_4hours", "m_Kirita_12hours"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Doke_Control", "m_Balzer_AKI_14d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Doke_Control", "m_Balzer_Control"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Abedini_Control", "h_Lake_LD"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_CKD_3d", "m_Balzer_AKI_3d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_AKI_3d", "m_Balzer_AKI_14d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Hinze_Control", "m_Kirita_Control"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Hinze_AKI", "h_Lake_AKI"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_AKI_1d", "m_Balzer_CKD_1d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Hinze_Control", "h_Abedini_Control"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_AKI_1d", "m_Kirita_12hours"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_Control", "m_Doke_UUO"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_AKI_3d", "m_Doke_UUO"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_CKD_1d", "h_Lake_LD"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_CKD_1d", "m_Balzer_CKD_3d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Balzer_CKD_1d", "m_Balzer_CKD_14d"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Kirita_Control", "m_Balzer_Control"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Hinze_AKI", "h_Lake_LD"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Hinze_AKI", "h_Abedini_CKD"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("h_Lake_AKI","h_Abedini_CKD"),
                          color="red", active=T),
                     list(query=intersects,
                          params=list("m_Doke_Control", "h_Lake_CKD"),
                          color="red", active=T)
                ))
dev.off()





#---compare with Clusterprofiler----
#calculate enrichment for all clusters at once
cp_mice_x_human = clusterProfiler::compareCluster(geneCluster = mice_x_human, 
                                                  fun = "enricher", 
                                                  TERM2GENE = msig_combined)


#.. number of enriched terms found for each gene cluster:
#..   m_Kirita_Control: 16 
#..   h_Hinze_Control: 33 
#..   h_Lake_LD: 173 
#..   m_Balzer_Control: 158 
#..   m_Doke_Control: 179 
#..   h_Abedini_Control: 223 
#..   m_Balzer_AKI_14d: 213 
#..   m_Kirita_4hours: 333 
#..   m_Kirita_12hours: 146 
#..   m_Balzer_AKI_1d: 248 
#..   h_Lake_AKI: 374 
#..   h_Hinze_AKI: 298 
#..   m_Balzer_AKI_3d: 354 
#..   h_Abedini_CKD: 202 
#..   m_Balzer_CKD_3d: 470 
#..   h_Lake_CKD: 122 
#..   m_Doke_UUO: 39 
#..   m_Balzer_CKD_1d: 514 
#..   m_Balzer_CKD_14d: 953 



#---heatmap-------------------------------------------------------------------
#filter results by q-value
qvalue_cutoff = 0.05
cp_mice_x_human@compareClusterResult = cp_mice_x_human@compareClusterResult[which(cp_mice_x_human@compareClusterResult$qvalue < qvalue_cutoff),]

#reorganize the qvalue result in matrix form, single cell clusters in columns, pathways in rows
cp_mice_x_human_mat = xtabs(qvalue ~ ID + Cluster, data = cp_mice_x_human@compareClusterResult)

#transform q-value so that low values (highly significant) have higher numbers
cp_mice_x_human_mat = -log10(cp_mice_x_human_mat)
cp_mice_x_human_mat[is.infinite(cp_mice_x_human_mat)] = 0


#top5 per cluster
test_mice_x_human <- data.frame(cp_mice_x_human_mat) #43.491   3

subset_mouse_Kirita_Control = subset(test_mice_x_human, Cluster == "m_Kirita_Control")
subset_mouse_Kirita_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_Control #5   3

subset_mouse_Kirita_4hours = subset(test_mice_x_human, Cluster == "m_Kirita_4hours")
subset_mouse_Kirita_4hours %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_4hours #5   3

subset_mouse_Kirita_12hours = subset(test_mice_x_human, Cluster == "m_Kirita_12hours")
subset_mouse_Kirita_12hours %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_12hours #5   3

#subset_mouse_Kirita_2d = subset(test_mice_x_human, Cluster == "m_Kirita_2d")
#subset_mouse_Kirita_2d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_2d 

#subset_mouse_Kirita_14d = subset(test_mice_x_human, Cluster == "m_Kirita_14d")
#subset_mouse_Kirita_14d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_14d

#subset_mouse_Kirita_6weeks = subset(test_mice_x_human, Cluster == "m_Kirita_6weeks")
#subset_mouse_Kirita_6weeks %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Kirita_6weeks #11   3


subset_mouse_Doke_Control = subset(test_mice_x_human, Cluster == "m_Doke_Control")
subset_mouse_Doke_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Doke_Control #5   3

subset_mouse_Doke_UUO = subset(test_mice_x_human, Cluster == "m_Doke_UUO")
subset_mouse_Doke_UUO %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Doke_UUO #5   3

subset_mouse_Balzer_Control = subset(test_mice_x_human, Cluster == "m_Balzer_Control")
subset_mouse_Balzer_Control %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_Control #5   3

subset_mouse_Balzer_AKI_1d = subset(test_mice_x_human, Cluster == "m_Balzer_AKI_1d")
subset_mouse_Balzer_AKI_1d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_AKI_1d #7   3

subset_mouse_Balzer_AKI_3d = subset(test_mice_x_human, Cluster == "m_Balzer_AKI_3d")
subset_mouse_Balzer_AKI_3d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_AKI_3d #7   3

subset_mouse_Balzer_AKI_14d = subset(test_mice_x_human, Cluster == "m_Balzer_AKI_14d")
subset_mouse_Balzer_AKI_14d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_AKI_14d #7   3

subset_mouse_Balzer_CKD_1d = subset(test_mice_x_human, Cluster == "m_Balzer_CKD_1d")
subset_mouse_Balzer_CKD_1d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_CKD_1d #7   3

subset_mouse_Balzer_CKD_3d = subset(test_mice_x_human, Cluster == "m_Balzer_CKD_3d")
subset_mouse_Balzer_CKD_3d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_CKD_3d #7   3

subset_mouse_Balzer_CKD_14d = subset(test_mice_x_human, Cluster == "m_Balzer_CKD_14d")
subset_mouse_Balzer_CKD_14d %>% top_n(n = 5, wt = Freq) -> topSOI_mouse_Balzer_CKD_14d #7   3

subset_human_Lake_Control = subset(test_mice_x_human, Cluster == "h_Lake_LD")
subset_human_Lake_Control %>% top_n(n = 5, wt = Freq) -> topSOI_human_Lake_LD #5   3

subset_human_Lake_AKI = subset(test_mice_x_human, Cluster == "h_Lake_AKI")
subset_human_Lake_AKI %>% top_n(n = 5, wt = Freq) -> topSOI_human_Lake_AKI #5   3

subset_human_Lake_CKD = subset(test_mice_x_human, Cluster == "h_Lake_CKD")
subset_human_Lake_CKD %>% top_n(n = 5, wt = Freq) -> topSOI_human_Lake_CKD #5   3

subset_human_Hinze_Control = subset(test_mice_x_human, Cluster == "h_Hinze_Control")
subset_human_Hinze_Control %>% top_n(n = 5, wt = Freq) -> topSOI_human_Hinze_Control #5   3

subset_human_Hinze_AKI = subset(test_mice_x_human, Cluster == "h_Hinze_AKI")
subset_human_Hinze_AKI %>%top_n(n = 5, wt = Freq) -> topSOI_human_Hinze_AKI #5   3

subset_human_Abedini_Control = subset(test_mice_x_human, cluster= "human_Abedini_Control")
subset_human_Abedini_Control %>% top_n(n=5, wt=Freq) -> topSOI_human_Abedini_Control #5   3

subset_human_Abedini_CKD<-subset(test_mice_x_human, cluster="human_Abedini_CKD")
subset_human_Abedini_CKD %>% top_n(n=5, wt=Freq) -> topSOI_human_Abedini_CKD #5   3



topSOI_all <- rbind(topSOI_human_Hinze_Control, topSOI_mouse_Kirita_Control,
                    topSOI_mouse_Doke_Control, topSOI_mouse_Balzer_Control, 
                    topSOI_human_Lake_LD, topSOI_mouse_Balzer_AKI_14d,
                    topSOI_mouse_Kirita_4hours, topSOI_mouse_Kirita_12hours,
                    topSOI_mouse_Balzer_AKI_1d, topSOI_mouse_Balzer_AKI_3d,
                    topSOI_human_Lake_AKI, topSOI_human_Hinze_AKI,
                    topSOI_mouse_Balzer_CKD_1d, topSOI_mouse_Balzer_CKD_3d,
                    topSOI_mouse_Balzer_CKD_14d, topSOI_human_Lake_CKD, topSOI_mouse_Doke_UUO,
                    topSOI_human_Abedini_Control, topSOI_human_Abedini_CKD)


#these are the combined top5 enriched terms per cluster, not unique!
# topSOI_all$ID
#87   3


#subset cp_mice_x_human_mat to these top5 per cluster
cp_mice_x_human_mat <- cp_mice_x_human_mat[unique(topSOI_all$ID),]

#any q-value lower than 10e-4 is clamped to produce a more readable heatmap
cp_mice_x_human_mat[cp_mice_x_human_mat > 6] = 6
dim(cp_mice_x_human_mat) #58   19


#---color cols by single-cell clusters----
cp_mice_x_human_mat2 <- cp_mice_x_human_mat[rownames(cp_mice_x_human_mat)[order(rownames(cp_mice_x_human_mat))],]
#dim(cp_mice_x_human_mat2) #58   17


# order of columnss (= Cluster):
# 7 Controls: m_Kirita_Control, h_Hinze_Control, h_Lake_LD, m_Balzer_Control, m_Doke_Control, h_Abedini_Control, m_Balzer_AKI_14d
# 6 AKIs: m_Kirita_4hours, m_Kirita_12hours, m_Balzer_AKI_1d, h_Lake_AKI, h_Hinze_AKI, m_Balzer_AKI_3d
# 5 CKDs: h_Abedini_CKD, m_Balzer_CKD_3d, h_Lake_CKD, m_Doke_UUO, #m_Balzer_CKD_1d, m_Balzer_CKD_14d


cluster_colors <- c("#336600","#336600", "#336600", "#336600",
                             "#336600", "#336600", "yellow",#7 Control
                             "#FF0000", "#FF0000", "#FF0000", 
                             "#FF0000", "#FF0000", "#FF0000", #6 aki
                             "#660000", "#660000", "#660000",
                             "#660000", "#660000", "#660000") #6 ckd
                             

names(cluster_colors) <- c("Control", "Control", "Control", "Control",
                           "Control", "Control", "AKI recovered",  #7
                           "AKI", "AKI", "AKI",
                           "AKI", "AKI", "AKI", #6
                           "CKD / severe AKI", "CKD / severe AKI", "CKD / severe AKI",
                           "CKD / severe AKI", "CKD / severe AKI", "CKD / severe AKI") #6


species_colors <- c("#F8C471","#85C1E9", "#85C1E9", "#F8C471",
                             "#F8C471", "#85C1E9", "#F8C471",#7 Control
                             "#F8C471", "#F8C471", "#F8C471", 
                             "#85C1E9", "#85C1E9", "#F8C471", #6 aki
                             "#85C1E9", "#F8C471", "#85C1E9",
                             "#F8C471", "#F8C471", "#F8C471") #6 ckd
                             

names(species_colors) <- c("mouse", "human", "human", "mouse",
                           "mouse", "human", "mouse",  #7
                           "mouse", "mouse", "mouse",
                           "human", "human", "mouse", #6
                           "human", "mouse", "human",
                           "mouse", "mouse", "mouse") #6



column_ha_cluster <- HeatmapAnnotation(
  Data_set = names(cluster_colors),
  col = list(Data_set = cluster_colors),
  simple_anno_size = unit(5, "mm"),
  annotation_name_gp = gpar(fontsize = 20),
  annotation_name_side = "left",
  show_annotation_name=FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Data_set = list(
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
  show_annotation_name=FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    Species = list(
      title_position= "topleft",
      title_gp = gpar(fontsize = 20, 
                      fontface = "bold"),
      labels_gp = gpar(fontsize = 20))))


# rename rows (pathways) with new names defined in "new_rownames.csv":
new_rownames<-read.csv("new_rownames.csv", header=T)

# charakter vektor of rownames(cp_mice_x_human_mat2)
rownames_original<-rownames(cp_mice_x_human_mat2)

# df[match(target, df$name),]
new_rownames<-new_rownames[match(rownames_original, new_rownames$ID),]

# df2 <- subset(df, select = -c(id, name, chapters))
new_rownames$ID<-NULL   # short alternative to subset()

# rename column
names(new_rownames)[names(new_rownames) == "new_rownames"] <- "ID"

# replacing rownames by new rownames (in charakter vektor)
rownames(cp_mice_x_human_mat2) <- new_rownames$ID



# remove rows (pathways), that are similiar to others:
cp_mice_x_human_mat3<-cp_mice_x_human_mat2[
  -c(1,4:6,8,10,14:15,17:18,20:22,24,32:34,36:37,40,42,44:45,48:50,55,57),]


# rename columns
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Abedini_Control"] <- "Abedini"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Abedini_CKD"] <- "Abedini"

colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Hinze_Control"] <- "Hinze"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Hinze_AKI"] <- "Hinze"

colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Lake_LD"] <- "Lake"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Lake_AKI"] <- "Lake"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "h_Lake_CKD"] <- "Lake"

colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Doke_Control"] <- "Doke"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Doke_UUO"] <- "Doke"

colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Kirita_Control"] <- "Kirita"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Kirita_4hours"] <- "Kirita_4hrs"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Kirita_12hours"] <- "Kirita_12hrs"

colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_Control"] <- "Balzer"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_AKI_1d"] <- "Balzer_1d"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_AKI_3d"] <- "Balzer_3d"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_AKI_14d"] <- "Balzer_14d"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_CKD_1d"] <- "Balzer_1d"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_CKD_3d"] <- "Balzer_3d"
colnames(cp_mice_x_human_mat3)[colnames(cp_mice_x_human_mat3) == "m_Balzer_CKD_14d"] <- "Balzer_14d"


#---HEATMAP----

heatmap<-ComplexHeatmap::Heatmap(
  cp_mice_x_human_mat3, "-log10(q)",
  na_col="white",
  col= colorRampPalette(c("white", "#241346"))(n=200),
  cluster_rows=F,
  row_order = c(13,16,21,14,8,12,28,27,9,12,30,3:4,6,18,17,5,19,22,24,11,2,26,7,20,15,29,23,1,10),
  width= unit(30, "cm"), height=unit(30,"cm"),
  row_names_gp=gpar(fontsize=28),
  cluster_columns= F,
  column_order = c(1:2,5,4,6,3,7, #Control
                   8:13, #aki
                   17,14,15,18,19,16), #ckd
  top_annotation = column_ha_cluster,
  bottom_annotation= column_ha_species,
  show_heatmap_legend = T,
  column_names_gp=gpar(fontsize=28),
  column_names_rot = 45,
  use_raster = TRUE,
  heatmap_legend_param = list(title = "-log10(q)",
                              title_gp = gpar(col = "black", fontsize = 20, fontface="bold"),
                              legend_height = unit(16, "cm"),
                              legend_width  = unit(16,"cm"),
                              title_position = "topleft",
                              direction = "vertical",
                              labels_gp = gpar(fontsize= 20, fontface = "bold")))


heatmap=draw(heatmap, heatmap_legend_side = "left", annotation_legend_side = "left")

pdf('HEATMAP.pdf', width=26, height=16)
heatmap
dev.off()

