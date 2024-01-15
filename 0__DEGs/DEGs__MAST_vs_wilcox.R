set.seed(123)
setwd('/.../0__DEG')

#==========================================================
#======================== LOAD DEGs =======================
#==========================================================

#======================== Abedini =======================
#Abedini_Health MAST
Abedini_MAST_deg <- read.csv('Abedini__human__DEGs_MAST.csv')
Abedini_Health_MAST_COI_deg <- Abedini_MAST_deg[Abedini_MAST_deg$cluster=="Health",]
Abedini_Health_MAST_COI_deg <- Abedini_Health_MAST_COI_deg[Abedini_Health_MAST_COI_deg$p_val_adj <0.05,]
Abedini_Health_MAST_COI_deg <- Abedini_Health_MAST_COI_deg[Abedini_Health_MAST_COI_deg$avg_log2FC >0,]
Abedini_Health_MAST_COI_deg <- Abedini_Health_MAST_COI_deg[order(Abedini_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Abedini_Health_MAST_COI_deg$gene)) #479
Abedini_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Abedini_Health_MAST_COI_deg_top100

#Abedini_Health Wilcoxon
Abedini_wilcox_deg <- read.csv('Abedini__human__DEGs_wilcox.csv')
Abedini_Health_wilcox_COI_deg <- Abedini_wilcox_deg[Abedini_wilcox_deg$cluster=="Health",]
Abedini_Health_wilcox_COI_deg <- Abedini_Health_wilcox_COI_deg[Abedini_Health_wilcox_COI_deg$p_val_adj <0.05,]
Abedini_Health_wilcox_COI_deg <- Abedini_Health_wilcox_COI_deg[Abedini_Health_wilcox_COI_deg$avg_log2FC >0,]
Abedini_Health_wilcox_COI_deg <- Abedini_Health_wilcox_COI_deg[order(Abedini_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Abedini_Health_wilcox_COI_deg$gene)) #479
Abedini_Health_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Abedini_Health_wilcox_COI_deg_top100

#Abedini_Health intersect MAST & Wilcoxon
length(intersect(unique(Abedini_Health_MAST_COI_deg$gene), unique(Abedini_Health_wilcox_COI_deg$gene))) #479

#Abedini_CKD MAST
Abedini_MAST_deg <- read.csv('Abedini__human__DEGs_MAST.csv')
Abedini_CKD_MAST_COI_deg <- Abedini_MAST_deg[Abedini_MAST_deg$cluster=="CKD",]
Abedini_CKD_MAST_COI_deg <- Abedini_CKD_MAST_COI_deg[Abedini_CKD_MAST_COI_deg$p_val_adj <0.05,]
Abedini_CKD_MAST_COI_deg <- Abedini_CKD_MAST_COI_deg[Abedini_CKD_MAST_COI_deg$avg_log2FC >0,]
Abedini_CKD_MAST_COI_deg <- Abedini_CKD_MAST_COI_deg[order(Abedini_CKD_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Abedini_CKD_MAST_COI_deg$gene)) #593
Abedini_CKD_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Abedini_CKD_MAST_COI_deg_top100

#Abedini_CKD Wilcoxon
Abedini_wilcox_deg <- read.csv('Abedini__human__DEGs_wilcox.csv')
Abedini_CKD_wilcox_COI_deg <- Abedini_wilcox_deg[Abedini_wilcox_deg$cluster=="CKD",]
Abedini_CKD_wilcox_COI_deg <- Abedini_CKD_wilcox_COI_deg[Abedini_CKD_wilcox_COI_deg$p_val_adj <0.05,]
Abedini_CKD_wilcox_COI_deg <- Abedini_CKD_wilcox_COI_deg[Abedini_CKD_wilcox_COI_deg$avg_log2FC >0,]
Abedini_CKD_wilcox_COI_deg <- Abedini_CKD_wilcox_COI_deg[order(Abedini_CKD_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Abedini_CKD_wilcox_COI_deg$gene)) #494
Abedini_CKD_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Abedini_CKD_wilcox_COI_deg_top100

#Abedini_CKD intersect MAST & Wilcoxon
length(intersect(unique(Abedini_CKD_MAST_COI_deg$gene), unique(Abedini_CKD_wilcox_COI_deg$gene))) #494



#======================== Balzer =======================
#Balzer_Health MAST humanized
Balzer_MAST_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_MAST_HUMAN_genenames.csv')
Balzer_Health_MAST_COI_deg <- Balzer_MAST_HUMAN_deg[Balzer_MAST_HUMAN_deg$cluster=="Health",]
Balzer_Health_MAST_COI_deg <- Balzer_Health_MAST_COI_deg[Balzer_Health_MAST_COI_deg$p_val_adj <0.05,]
Balzer_Health_MAST_COI_deg <- Balzer_Health_MAST_COI_deg[Balzer_Health_MAST_COI_deg$avg_log2FC >0,]
Balzer_Health_MAST_COI_deg <- Balzer_Health_MAST_COI_deg[order(Balzer_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_Health_MAST_COI_deg$gene)) #1176
Balzer_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_Health_MAST_COI_deg_top100

#Balzer_Health wilcox humanized
Balzer_wilcox_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Balzer_Health_wilcox_COI_deg <- Balzer_wilcox_HUMAN_deg[Balzer_wilcox_HUMAN_deg$cluster=="Health",]
Balzer_Health_wilcox_COI_deg <- Balzer_Health_wilcox_COI_deg[Balzer_Health_wilcox_COI_deg$p_val_adj <0.05,]
Balzer_Health_wilcox_COI_deg <- Balzer_Health_wilcox_COI_deg[Balzer_Health_wilcox_COI_deg$avg_log2FC >0,]
Balzer_Health_wilcox_COI_deg <- Balzer_Health_wilcox_COI_deg[order(Balzer_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_Health_wilcox_COI_deg$gene)) #994
Balzer_Health_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_Health_wilcox_COI_deg_top100

#Balzer_Health intersect MAST & Wilcoxon humanized
length(intersect(unique(Balzer_Health_MAST_COI_deg$gene), unique(Balzer_Health_wilcox_COI_deg$gene))) #994

#Balzer_AKI MAST humanized
Balzer_MAST_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_MAST_HUMAN_genenames.csv')
Balzer_AKI_MAST_COI_deg <- Balzer_MAST_HUMAN_deg[Balzer_MAST_HUMAN_deg$cluster=="AKI",]
Balzer_AKI_MAST_COI_deg <- Balzer_AKI_MAST_COI_deg[Balzer_AKI_MAST_COI_deg$p_val_adj <0.05,]
Balzer_AKI_MAST_COI_deg <- Balzer_AKI_MAST_COI_deg[Balzer_AKI_MAST_COI_deg$avg_log2FC >0,]
Balzer_AKI_MAST_COI_deg <- Balzer_AKI_MAST_COI_deg[order(Balzer_AKI_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_AKI_MAST_COI_deg$gene)) #750
Balzer_AKI_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_AKI_MAST_COI_deg_top100

#Balzer_AKI wilcox humanized
Balzer_wilcox_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Balzer_AKI_wilcox_COI_deg <- Balzer_wilcox_HUMAN_deg[Balzer_wilcox_HUMAN_deg$cluster=="AKI",]
Balzer_AKI_wilcox_COI_deg <- Balzer_AKI_wilcox_COI_deg[Balzer_AKI_wilcox_COI_deg$p_val_adj <0.05,]
Balzer_AKI_wilcox_COI_deg <- Balzer_AKI_wilcox_COI_deg[Balzer_AKI_wilcox_COI_deg$avg_log2FC >0,]
Balzer_AKI_wilcox_COI_deg <- Balzer_AKI_wilcox_COI_deg[order(Balzer_AKI_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_AKI_wilcox_COI_deg$gene)) #750
Balzer_AKI_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_AKI_wilcox_COI_deg_top100

#Balzer_AKI intersect MAST & Wilcoxon humanized
length(intersect(unique(Balzer_AKI_MAST_COI_deg$gene), unique(Balzer_AKI_wilcox_COI_deg$gene))) #750

#Balzer_Recovery MAST humanized
Balzer_MAST_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_MAST_HUMAN_genenames.csv')
Balzer_Recovery_MAST_COI_deg <- Balzer_MAST_HUMAN_deg[Balzer_MAST_HUMAN_deg$cluster=="Recovery",]
Balzer_Recovery_MAST_COI_deg <- Balzer_Recovery_MAST_COI_deg[Balzer_Recovery_MAST_COI_deg$p_val_adj <0.05,]
Balzer_Recovery_MAST_COI_deg <- Balzer_Recovery_MAST_COI_deg[Balzer_Recovery_MAST_COI_deg$avg_log2FC >0,]
Balzer_Recovery_MAST_COI_deg <- Balzer_Recovery_MAST_COI_deg[order(Balzer_Recovery_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_Recovery_MAST_COI_deg$gene)) #229
Balzer_Recovery_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_Recovery_MAST_COI_deg_top100

#Balzer_Recovery wilcox humanized
Balzer_wilcox_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Balzer_Recovery_wilcox_COI_deg <- Balzer_wilcox_HUMAN_deg[Balzer_wilcox_HUMAN_deg$cluster=="Recovery",]
Balzer_Recovery_wilcox_COI_deg <- Balzer_Recovery_wilcox_COI_deg[Balzer_Recovery_wilcox_COI_deg$p_val_adj <0.05,]
Balzer_Recovery_wilcox_COI_deg <- Balzer_Recovery_wilcox_COI_deg[Balzer_Recovery_wilcox_COI_deg$avg_log2FC >0,]
Balzer_Recovery_wilcox_COI_deg <- Balzer_Recovery_wilcox_COI_deg[order(Balzer_Recovery_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_Recovery_wilcox_COI_deg$gene)) #225
Balzer_Recovery_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_Recovery_wilcox_COI_deg_top100

#Balzer_Recovery intersect MAST & Wilcoxon humanized
length(intersect(unique(Balzer_Recovery_MAST_COI_deg$gene), unique(Balzer_Recovery_wilcox_COI_deg$gene))) #225

#Balzer_CKD MAST humanized
Balzer_MAST_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_MAST_HUMAN_genenames.csv')
Balzer_CKD_MAST_COI_deg <- Balzer_MAST_HUMAN_deg[Balzer_MAST_HUMAN_deg$cluster=="CKD",]
Balzer_CKD_MAST_COI_deg <- Balzer_CKD_MAST_COI_deg[Balzer_CKD_MAST_COI_deg$p_val_adj <0.05,]
Balzer_CKD_MAST_COI_deg <- Balzer_CKD_MAST_COI_deg[Balzer_CKD_MAST_COI_deg$avg_log2FC >0,]
Balzer_CKD_MAST_COI_deg <- Balzer_CKD_MAST_COI_deg[order(Balzer_CKD_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_CKD_MAST_COI_deg$gene)) #1321
Balzer_CKD_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_CKD_MAST_COI_deg_top100

#Balzer_CKD wilcox humanized
Balzer_wilcox_HUMAN_deg <- read.csv('Balzer__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Balzer_CKD_wilcox_COI_deg <- Balzer_wilcox_HUMAN_deg[Balzer_wilcox_HUMAN_deg$cluster=="CKD",]
Balzer_CKD_wilcox_COI_deg <- Balzer_CKD_wilcox_COI_deg[Balzer_CKD_wilcox_COI_deg$p_val_adj <0.05,]
Balzer_CKD_wilcox_COI_deg <- Balzer_CKD_wilcox_COI_deg[Balzer_CKD_wilcox_COI_deg$avg_log2FC >0,]
Balzer_CKD_wilcox_COI_deg <- Balzer_CKD_wilcox_COI_deg[order(Balzer_CKD_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Balzer_CKD_wilcox_COI_deg$gene)) #1320
Balzer_CKD_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Balzer_CKD_wilcox_COI_deg_top100

#Balzer_CKD intersect MAST & Wilcoxon humanized
length(intersect(unique(Balzer_CKD_MAST_COI_deg$gene), unique(Balzer_CKD_wilcox_COI_deg$gene))) #1320



#======================== Doke =======================
#Doke_Health MAST humanized
Doke_MAST_HUMAN_deg <- read.csv('Doke__mouse__DEGs_MAST_HUMAN_genenames.csv')
Doke_Health_MAST_COI_deg <- Doke_MAST_HUMAN_deg[Doke_MAST_HUMAN_deg$cluster=="Health",]
Doke_Health_MAST_COI_deg <- Doke_Health_MAST_COI_deg[Doke_Health_MAST_COI_deg$p_val_adj <0.05,]
Doke_Health_MAST_COI_deg <- Doke_Health_MAST_COI_deg[Doke_Health_MAST_COI_deg$avg_log2FC >0,]
Doke_Health_MAST_COI_deg <- Doke_Health_MAST_COI_deg[order(Doke_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Doke_Health_MAST_COI_deg$gene)) #110
Doke_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Doke_Health_MAST_COI_deg_top100

#Doke_Health wilcox humanized
Doke_wilcox_HUMAN_deg <- read.csv('Doke__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Doke_Health_wilcox_COI_deg <- Doke_wilcox_HUMAN_deg[Doke_wilcox_HUMAN_deg$cluster=="Health",]
Doke_Health_wilcox_COI_deg <- Doke_Health_wilcox_COI_deg[Doke_Health_wilcox_COI_deg$p_val_adj <0.05,]
Doke_Health_wilcox_COI_deg <- Doke_Health_wilcox_COI_deg[Doke_Health_wilcox_COI_deg$avg_log2FC >0,]
Doke_Health_wilcox_COI_deg <- Doke_Health_wilcox_COI_deg[order(Doke_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Doke_Health_wilcox_COI_deg$gene)) #110
Doke_Health_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Doke_Health_wilcox_COI_deg_top100

#Doke_Health intersect MAST & Wilcoxon humanized
length(intersect(unique(Doke_Health_MAST_COI_deg$gene), unique(Doke_Health_wilcox_COI_deg$gene))) #110

#Doke_CKD MAST humanized
Doke_MAST_HUMAN_deg <- read.csv('Doke__mouse__DEGs_MAST_HUMAN_genenames.csv')
Doke_CKD_MAST_COI_deg <- Doke_MAST_HUMAN_deg[Doke_MAST_HUMAN_deg$cluster=="CKD",]
Doke_CKD_MAST_COI_deg <- Doke_CKD_MAST_COI_deg[Doke_CKD_MAST_COI_deg$p_val_adj <0.05,]
Doke_CKD_MAST_COI_deg <- Doke_CKD_MAST_COI_deg[Doke_CKD_MAST_COI_deg$avg_log2FC >0,]
Doke_CKD_MAST_COI_deg <- Doke_CKD_MAST_COI_deg[order(Doke_CKD_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Doke_CKD_MAST_COI_deg$gene)) #260
Doke_CKD_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Doke_CKD_MAST_COI_deg_top100

#Doke_CKD wilcox humanized
Doke_wilcox_HUMAN_deg <- read.csv('Doke__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Doke_CKD_wilcox_COI_deg <- Doke_wilcox_HUMAN_deg[Doke_wilcox_HUMAN_deg$cluster=="CKD",]
Doke_CKD_wilcox_COI_deg <- Doke_CKD_wilcox_COI_deg[Doke_CKD_wilcox_COI_deg$p_val_adj <0.05,]
Doke_CKD_wilcox_COI_deg <- Doke_CKD_wilcox_COI_deg[Doke_CKD_wilcox_COI_deg$avg_log2FC >0,]
Doke_CKD_wilcox_COI_deg <- Doke_CKD_wilcox_COI_deg[order(Doke_CKD_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Doke_CKD_wilcox_COI_deg$gene)) #260
Doke_CKD_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Doke_CKD_wilcox_COI_deg_top100

#Doke_CKD intersect MAST & Wilcoxon humanized
length(intersect(unique(Doke_CKD_MAST_COI_deg$gene), unique(Doke_CKD_wilcox_COI_deg$gene))) #260



#======================== Hinze =======================
#Hinze_Health MAST
Hinze_MAST_deg <- read.csv('Hinze__human__DEGs_MAST.csv')
Hinze_Health_MAST_COI_deg <- Hinze_MAST_deg[Hinze_MAST_deg$cluster=="Health",]
Hinze_Health_MAST_COI_deg <- Hinze_Health_MAST_COI_deg[Hinze_Health_MAST_COI_deg$p_val_adj <0.05,]
Hinze_Health_MAST_COI_deg <- Hinze_Health_MAST_COI_deg[Hinze_Health_MAST_COI_deg$avg_log2FC >0,]
Hinze_Health_MAST_COI_deg <- Hinze_Health_MAST_COI_deg[order(Hinze_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Hinze_Health_MAST_COI_deg$gene)) #1208
Hinze_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Hinze_Health_MAST_COI_deg_top100

#Hinze_Health Wilcoxon
Hinze_wilcox_deg <- read.csv('Hinze__human__DEGs_wilcox.csv')
Hinze_Health_wilcox_COI_deg <- Hinze_wilcox_deg[Hinze_wilcox_deg$cluster=="Health",]
Hinze_Health_wilcox_COI_deg <- Hinze_Health_wilcox_COI_deg[Hinze_Health_wilcox_COI_deg$p_val_adj <0.05,]
Hinze_Health_wilcox_COI_deg <- Hinze_Health_wilcox_COI_deg[Hinze_Health_wilcox_COI_deg$avg_log2FC >0,]
Hinze_Health_wilcox_COI_deg <- Hinze_Health_wilcox_COI_deg[order(Hinze_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Hinze_Health_wilcox_COI_deg$gene)) #1204

#Hinze_Health intersect MAST & Wilcoxon
length(intersect(unique(Hinze_Health_MAST_COI_deg$gene), unique(Hinze_Health_wilcox_COI_deg$gene))) #1204

#Hinze_AKI MAST
Hinze_MAST_deg <- read.csv('Hinze__human__DEGs_MAST.csv')
Hinze_AKI_MAST_COI_deg <- Hinze_MAST_deg[Hinze_MAST_deg$cluster=="AKI",]
Hinze_AKI_MAST_COI_deg <- Hinze_AKI_MAST_COI_deg[Hinze_AKI_MAST_COI_deg$p_val_adj <0.05,]
Hinze_AKI_MAST_COI_deg <- Hinze_AKI_MAST_COI_deg[Hinze_AKI_MAST_COI_deg$avg_log2FC >0,]
Hinze_AKI_MAST_COI_deg <- Hinze_AKI_MAST_COI_deg[order(Hinze_AKI_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Hinze_AKI_MAST_COI_deg$gene)) #515
Hinze_AKI_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Hinze_AKI_MAST_COI_deg_top100

#Hinze_AKI Wilcoxon
Hinze_wilcox_deg <- read.csv('Hinze__human__DEGs_wilcox.csv')
Hinze_AKI_wilcox_COI_deg <- Hinze_wilcox_deg[Hinze_wilcox_deg$cluster=="AKI",]
Hinze_AKI_wilcox_COI_deg <- Hinze_AKI_wilcox_COI_deg[Hinze_AKI_wilcox_COI_deg$p_val_adj <0.05,]
Hinze_AKI_wilcox_COI_deg <- Hinze_AKI_wilcox_COI_deg[Hinze_AKI_wilcox_COI_deg$avg_log2FC >0,]
Hinze_AKI_wilcox_COI_deg <- Hinze_AKI_wilcox_COI_deg[order(Hinze_AKI_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Hinze_AKI_wilcox_COI_deg$gene)) #513

#Hinze_AKI_AKI intersect MAST & Wilcoxon
length(intersect(unique(Hinze_AKI_MAST_COI_deg$gene), unique(Hinze_AKI_wilcox_COI_deg$gene))) #513



#======================== Kirita =======================
#Kirita_Health MAST humanized
Kirita_MAST_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_MAST_HUMAN_genenames.csv')
Kirita_Health_MAST_COI_deg <- Kirita_MAST_HUMAN_deg[Kirita_MAST_HUMAN_deg$cluster=="Health",]
Kirita_Health_MAST_COI_deg <- Kirita_Health_MAST_COI_deg[Kirita_Health_MAST_COI_deg$p_val_adj <0.05,]
Kirita_Health_MAST_COI_deg <- Kirita_Health_MAST_COI_deg[Kirita_Health_MAST_COI_deg$avg_log2FC >0,]
Kirita_Health_MAST_COI_deg <- Kirita_Health_MAST_COI_deg[order(Kirita_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_Health_MAST_COI_deg$gene)) #357
Kirita_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_Health_MAST_COI_deg_top100

#Kirita_Health wilcox humanized
Kirita_wilcox_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Kirita_Health_wilcox_COI_deg <- Kirita_wilcox_HUMAN_deg[Kirita_wilcox_HUMAN_deg$cluster=="Health",]
Kirita_Health_wilcox_COI_deg <- Kirita_Health_wilcox_COI_deg[Kirita_Health_wilcox_COI_deg$p_val_adj <0.05,]
Kirita_Health_wilcox_COI_deg <- Kirita_Health_wilcox_COI_deg[Kirita_Health_wilcox_COI_deg$avg_log2FC >0,]
Kirita_Health_wilcox_COI_deg <- Kirita_Health_wilcox_COI_deg[order(Kirita_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_Health_wilcox_COI_deg$gene)) #342
Kirita_Health_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_Health_wilcox_COI_deg_top100

#Kirita_Health intersect MAST & Wilcoxon humanized
length(intersect(unique(Kirita_Health_MAST_COI_deg$gene), unique(Kirita_Health_wilcox_COI_deg$gene))) #342

#Kirita_AKI MAST humanized
Kirita_MAST_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_MAST_HUMAN_genenames.csv')
Kirita_AKI_MAST_COI_deg <- Kirita_MAST_HUMAN_deg[Kirita_MAST_HUMAN_deg$cluster=="AKI",]
Kirita_AKI_MAST_COI_deg <- Kirita_AKI_MAST_COI_deg[Kirita_AKI_MAST_COI_deg$p_val_adj <0.05,]
Kirita_AKI_MAST_COI_deg <- Kirita_AKI_MAST_COI_deg[Kirita_AKI_MAST_COI_deg$avg_log2FC >0,]
Kirita_AKI_MAST_COI_deg <- Kirita_AKI_MAST_COI_deg[order(Kirita_AKI_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_AKI_MAST_COI_deg$gene)) #199
Kirita_AKI_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_AKI_MAST_COI_deg_top100

#Kirita_AKI wilcox humanized
Kirita_wilcox_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Kirita_AKI_wilcox_COI_deg <- Kirita_wilcox_HUMAN_deg[Kirita_wilcox_HUMAN_deg$cluster=="AKI",]
Kirita_AKI_wilcox_COI_deg <- Kirita_AKI_wilcox_COI_deg[Kirita_AKI_wilcox_COI_deg$p_val_adj <0.05,]
Kirita_AKI_wilcox_COI_deg <- Kirita_AKI_wilcox_COI_deg[Kirita_AKI_wilcox_COI_deg$avg_log2FC >0,]
Kirita_AKI_wilcox_COI_deg <- Kirita_AKI_wilcox_COI_deg[order(Kirita_AKI_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_AKI_wilcox_COI_deg$gene)) #199
Kirita_AKI_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_AKI_wilcox_COI_deg_top100

#Kirita_AKI intersect MAST & Wilcoxon humanized
length(intersect(unique(Kirita_AKI_MAST_COI_deg$gene), unique(Kirita_AKI_wilcox_COI_deg$gene))) #197

#Kirita_CKD MAST humanized
Kirita_MAST_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_MAST_HUMAN_genenames.csv')
Kirita_CKD_MAST_COI_deg <- Kirita_MAST_HUMAN_deg[Kirita_MAST_HUMAN_deg$cluster=="CKD",]
Kirita_CKD_MAST_COI_deg <- Kirita_CKD_MAST_COI_deg[Kirita_CKD_MAST_COI_deg$p_val_adj <0.05,]
Kirita_CKD_MAST_COI_deg <- Kirita_CKD_MAST_COI_deg[Kirita_CKD_MAST_COI_deg$avg_log2FC >0,]
Kirita_CKD_MAST_COI_deg <- Kirita_CKD_MAST_COI_deg[order(Kirita_CKD_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_CKD_MAST_COI_deg$gene)) #262
Kirita_CKD_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_CKD_MAST_COI_deg_top100

#Kirita_CKD wilcox humanized
Kirita_wilcox_HUMAN_deg <- read.csv('Kirita__mouse__DEGs_wilcox_HUMAN_genenames.csv')
Kirita_CKD_wilcox_COI_deg <- Kirita_wilcox_HUMAN_deg[Kirita_wilcox_HUMAN_deg$cluster=="CKD",]
Kirita_CKD_wilcox_COI_deg <- Kirita_CKD_wilcox_COI_deg[Kirita_CKD_wilcox_COI_deg$p_val_adj <0.05,]
Kirita_CKD_wilcox_COI_deg <- Kirita_CKD_wilcox_COI_deg[Kirita_CKD_wilcox_COI_deg$avg_log2FC >0,]
Kirita_CKD_wilcox_COI_deg <- Kirita_CKD_wilcox_COI_deg[order(Kirita_CKD_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Kirita_CKD_wilcox_COI_deg$gene)) #262
Kirita_CKD_wilcox_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Kirita_CKD_wilcox_COI_deg_top100

#Kirita_CKD intersect MAST & Wilcoxon humanized
length(intersect(unique(Kirita_CKD_MAST_COI_deg$gene), unique(Kirita_CKD_wilcox_COI_deg$gene))) #262



#======================== Lake =======================
#Lake_Health MAST
Lake_MAST_deg <- read.csv('Lake__human__DEGs_MAST.csv')
Lake_Health_MAST_COI_deg <- Lake_MAST_deg[Lake_MAST_deg$cluster=="Health",]
Lake_Health_MAST_COI_deg <- Lake_Health_MAST_COI_deg[Lake_Health_MAST_COI_deg$p_val_adj <0.05,]
Lake_Health_MAST_COI_deg <- Lake_Health_MAST_COI_deg[Lake_Health_MAST_COI_deg$avg_log2FC >0,]
Lake_Health_MAST_COI_deg <- Lake_Health_MAST_COI_deg[order(Lake_Health_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_Health_MAST_COI_deg$gene)) #402
Lake_Health_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Lake_Health_MAST_COI_deg_top100

#Lake_Health Wilcoxon
Lake_wilcox_deg <- read.csv('Lake__human__DEGs_wilcox.csv')
Lake_Health_wilcox_COI_deg <- Lake_wilcox_deg[Lake_wilcox_deg$cluster=="Health",]
Lake_Health_wilcox_COI_deg <- Lake_Health_wilcox_COI_deg[Lake_Health_wilcox_COI_deg$p_val_adj <0.05,]
Lake_Health_wilcox_COI_deg <- Lake_Health_wilcox_COI_deg[Lake_Health_wilcox_COI_deg$avg_log2FC >0,]
Lake_Health_wilcox_COI_deg <- Lake_Health_wilcox_COI_deg[order(Lake_Health_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_Health_wilcox_COI_deg$gene)) #402

#Lake_Health intersect MAST & Wilcoxon
length(intersect(unique(Lake_Health_MAST_COI_deg$gene), unique(Lake_Health_wilcox_COI_deg$gene))) #402

#Lake_AKI MAST
Lake_MAST_deg <- read.csv('Lake__human__DEGs_MAST.csv')
Lake_AKI_MAST_COI_deg <- Lake_MAST_deg[Lake_MAST_deg$cluster=="AKI",]
Lake_AKI_MAST_COI_deg <- Lake_AKI_MAST_COI_deg[Lake_AKI_MAST_COI_deg$p_val_adj <0.05,]
Lake_AKI_MAST_COI_deg <- Lake_AKI_MAST_COI_deg[Lake_AKI_MAST_COI_deg$avg_log2FC >0,]
Lake_AKI_MAST_COI_deg <- Lake_AKI_MAST_COI_deg[order(Lake_AKI_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_AKI_MAST_COI_deg$gene)) #593
Lake_AKI_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Lake_AKI_MAST_COI_deg_top100

#Lake_AKI Wilcoxon
Lake_wilcox_deg <- read.csv('Lake__human__DEGs_wilcox.csv')
Lake_AKI_wilcox_COI_deg <- Lake_wilcox_deg[Lake_wilcox_deg$cluster=="AKI",]
Lake_AKI_wilcox_COI_deg <- Lake_AKI_wilcox_COI_deg[Lake_AKI_wilcox_COI_deg$p_val_adj <0.05,]
Lake_AKI_wilcox_COI_deg <- Lake_AKI_wilcox_COI_deg[Lake_AKI_wilcox_COI_deg$avg_log2FC >0,]
Lake_AKI_wilcox_COI_deg <- Lake_AKI_wilcox_COI_deg[order(Lake_AKI_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_AKI_wilcox_COI_deg$gene)) #582

#Lake_AKI_AKI intersect MAST & Wilcoxon
length(intersect(unique(Lake_AKI_MAST_COI_deg$gene), unique(Lake_AKI_wilcox_COI_deg$gene))) #582

#Lake_CKD MAST
Lake_MAST_deg <- read.csv('Lake__human__DEGs_MAST.csv')
Lake_CKD_MAST_COI_deg <- Lake_MAST_deg[Lake_MAST_deg$cluster=="CKD",]
Lake_CKD_MAST_COI_deg <- Lake_CKD_MAST_COI_deg[Lake_CKD_MAST_COI_deg$p_val_adj <0.05,]
Lake_CKD_MAST_COI_deg <- Lake_CKD_MAST_COI_deg[Lake_CKD_MAST_COI_deg$avg_log2FC >0,]
Lake_CKD_MAST_COI_deg <- Lake_CKD_MAST_COI_deg[order(Lake_CKD_MAST_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_CKD_MAST_COI_deg$gene)) #583
Lake_CKD_MAST_COI_deg %>% top_n(n = 100, wt = avg_log2FC) -> Lake_CKD_MAST_COI_deg_top100

#Lake_CKD Wilcoxon
Lake_wilcox_deg <- read.csv('Lake__human__DEGs_wilcox.csv')
Lake_CKD_wilcox_COI_deg <- Lake_wilcox_deg[Lake_wilcox_deg$cluster=="CKD",]
Lake_CKD_wilcox_COI_deg <- Lake_CKD_wilcox_COI_deg[Lake_CKD_wilcox_COI_deg$p_val_adj <0.05,]
Lake_CKD_wilcox_COI_deg <- Lake_CKD_wilcox_COI_deg[Lake_CKD_wilcox_COI_deg$avg_log2FC >0,]
Lake_CKD_wilcox_COI_deg <- Lake_CKD_wilcox_COI_deg[order(Lake_CKD_wilcox_COI_deg$avg_log2FC, decreasing = T),]
length(unique(Lake_CKD_wilcox_COI_deg$gene)) #404

#Lake_CKD_CKD intersect MAST & Wilcoxon
length(intersect(unique(Lake_CKD_MAST_COI_deg$gene), unique(Lake_CKD_wilcox_COI_deg$gene))) #404



#==========================================================
#======================== VISUALIZE =======================
#==========================================================
library(ggVennDiagram)
library(ggplot2)
library(viridis)
library(Seurat)

#Abedini
pdf('venn_plot__Abedini_Health.pdf', width=2.5, height=2.5)
p_Abedini_Health <- ggVennDiagram(list(Abedini_Health_MAST = unique(Abedini_Health_MAST_COI_deg$gene), 
                                       Abedini_Health_Wilcox = unique(Abedini_Health_wilcox_COI_deg$gene)), 
                                  label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Abedini_Health
dev.off()

pdf('venn_plot__Abedini_CKD.pdf', width=2.5, height=2.5)
p_Abedini_CKD <- ggVennDiagram(list(Abedini_CKD_MAST = unique(Abedini_CKD_MAST_COI_deg$gene), 
                                    Abedini_CKD_Wilcox = unique(Abedini_CKD_wilcox_COI_deg$gene)), 
                               label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Abedini_CKD
dev.off()

#Balzer
pdf('venn_plot__Balzer_Health.pdf', width=2.5, height=2.5)
p_Balzer_Health <- ggVennDiagram(list(Balzer_Health_MAST = unique(Balzer_Health_MAST_COI_deg$gene), 
                                      Balzer_Health_Wilcox = unique(Balzer_Health_wilcox_COI_deg$gene)), 
                                 label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Balzer_Health
dev.off()

pdf('venn_plot__Balzer_AKI.pdf', width=2.5, height=2.5)
p_Balzer_AKI <- ggVennDiagram(list(Balzer_AKI_MAST = unique(Balzer_AKI_MAST_COI_deg$gene), 
                                   Balzer_AKI_Wilcox = unique(Balzer_AKI_wilcox_COI_deg$gene)), 
                              label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Balzer_AKI
dev.off()

pdf('venn_plot__Balzer_Recovery.pdf', width=2.5, height=2.5)
p_Balzer_Recovery <- ggVennDiagram(list(Balzer_Recovery_MAST = unique(Balzer_Recovery_MAST_COI_deg$gene), 
                                        Balzer_Recovery_Wilcox = unique(Balzer_Recovery_wilcox_COI_deg$gene)), 
                                   label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Balzer_Recovery
dev.off()

pdf('venn_plot__Balzer_CKD.pdf', width=2.5, height=2.5)
p_Balzer_CKD <- ggVennDiagram(list(Balzer_CKD_MAST = unique(Balzer_CKD_MAST_COI_deg$gene), 
                                   Balzer_CKD_Wilcox = unique(Balzer_CKD_wilcox_COI_deg$gene)), 
                              label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Balzer_CKD
dev.off()

#Doke
pdf('venn_plot__Doke_Health.pdf', width=2.5, height=2.5)
p_Doke_Health <- ggVennDiagram(list(Doke_Health_MAST = unique(Doke_Health_MAST_COI_deg$gene), 
                                    Doke_Health_Wilcox = unique(Doke_Health_wilcox_COI_deg$gene)), 
                               label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Doke_Health
dev.off()

pdf('venn_plot__Doke_CKD.pdf', width=2.5, height=2.5)
p_Doke_CKD <- ggVennDiagram(list(Doke_CKD_MAST = unique(Doke_CKD_MAST_COI_deg$gene), 
                                 Doke_CKD_Wilcox = unique(Doke_CKD_wilcox_COI_deg$gene)), 
                            label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Doke_CKD
dev.off()

#Hinze
pdf('venn_plot__Hinze_Health.pdf', width=2.5, height=2.5)
p_Hinze_Health <- ggVennDiagram(list(Hinze_Health_MAST = unique(Hinze_Health_MAST_COI_deg$gene), 
                                     Hinze_Health_Wilcox = unique(Hinze_Health_wilcox_COI_deg$gene)), 
                                label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Hinze_Health
dev.off()

pdf('venn_plot__Hinze_AKI.pdf', width=2.5, height=2.5)
p_Hinze_AKI <- ggVennDiagram(list(Hinze_AKI_MAST = unique(Hinze_AKI_MAST_COI_deg$gene), 
                                  Hinze_AKI_Wilcox = unique(Hinze_AKI_wilcox_COI_deg$gene)), 
                             label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Hinze_AKI
dev.off()

#Kirita
pdf('venn_plot__Kirita_Health.pdf', width=2.5, height=2.5)
p_Kirita_Health <- ggVennDiagram(list(Kirita_Health_MAST = unique(Kirita_Health_MAST_COI_deg$gene), 
                                      Kirita_Health_Wilcox = unique(Kirita_Health_wilcox_COI_deg$gene)), 
                                 label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Kirita_Health
dev.off()

pdf('venn_plot__Kirita_AKI.pdf', width=2.5, height=2.5)
p_Kirita_AKI <- ggVennDiagram(list(Kirita_AKI_MAST = unique(Kirita_AKI_MAST_COI_deg$gene), 
                                   Kirita_AKI_Wilcox = unique(Kirita_AKI_wilcox_COI_deg$gene)), 
                              label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Kirita_AKI
dev.off()

pdf('venn_plot__Kirita_CKD.pdf', width=2.5, height=2.5)
p_Kirita_CKD <- ggVennDiagram(list(Kirita_CKD_MAST = unique(Kirita_CKD_MAST_COI_deg$gene), 
                                   Kirita_CKD_Wilcox = unique(Kirita_CKD_wilcox_COI_deg$gene)), 
                              label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Kirita_CKD
dev.off()

#Lake
pdf('venn_plot__Lake_Health.pdf', width=2.5, height=2.5)
p_Lake_Health <- ggVennDiagram(list(Lake_Health_MAST = unique(Lake_Health_MAST_COI_deg$gene), 
                                    Lake_Health_Wilcox = unique(Lake_Health_wilcox_COI_deg$gene)), 
                               label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Lake_Health
dev.off()

pdf('venn_plot__Lake_AKI.pdf', width=2.5, height=2.5)
p_Lake_AKI <- ggVennDiagram(list(Lake_AKI_MAST = unique(Lake_AKI_MAST_COI_deg$gene), 
                                 Lake_AKI_Wilcox = unique(Lake_AKI_wilcox_COI_deg$gene)), 
                            label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Lake_AKI
dev.off()

pdf('venn_plot__Lake_CKD.pdf', width=2.5, height=2.5)
p_Lake_CKD <- ggVennDiagram(list(Lake_CKD_MAST = unique(Lake_CKD_MAST_COI_deg$gene), 
                                 Lake_CKD_Wilcox = unique(Lake_CKD_wilcox_COI_deg$gene)), 
                            label_alpha = 1, label='count', category.names = c("MAST", "Wilcoxon")) + scale_fill_viridis(option = "inferno", limits=c(0, 1320)) + NoLegend()
p_Lake_CKD
dev.off()
