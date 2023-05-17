#1 This takes the already present IRI PT dataset in diffusion map space.
#2 We then extract mouse orthologue genes from selected enriched gene sets.
#3 We then score this gene set in the IRI dataset.
#4 We project the enrichment score on the dataset in diffusion map space and create violin plots per cluster (exp.time) in this case

rm(list=ls())
library(slingshot)
library(Seurat)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)
library(RColorBrewer)
library(biomaRt)
set.seed(123)





#======================== get pathway of interest genes ===============================
pathwayofinterest <- "WP_NETWORK_MAP_OF_SARSCOV2_SIGNALING_PATHWAY"

#===== Load POI and get GOI
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
#msig__H= msigdbr::msigdbr(species = "Homo sapiens", category = "H") #hallmark gene sets
msig__C2 = msigdbr::msigdbr(species = "Homo sapiens", category = "C2") #curated gene sets
#msig__C5 = msigdbr::msigdbr(species = "Homo sapiens", category = "C5") #ontologygene sets

#msig_CP_BioCarta = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:BIOCARTA")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:BIOCARTA")])
#msig_CP_KEGG = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:KEGG")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:KEGG")])
#msig_CP_PID = data.frame(msig__C2$gs_name[which(msig__C2$gs_subcat == "CP:PID")], msig__C2$gene_symbol[which(msig__C2$gs_subcat == "CP:PID")])
#msig_CP_Reactome = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:REACTOME")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:REACTOME")])
msig_WP_WikiPathways = data.frame(msig__C2$gs_name[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")], msig__C2$gene_symbol[msig__C2$gs_subcat %in% c("CP:WIKIPATHWAYS")])
#msig_GO_BP = data.frame(msig__C5$gs_name[msig__C5$gs_subcat %in% c("GO:BP")], msig__C5$gene_symbol[msig__C5$gs_subcat %in% c("GO:BP")])
#msig_hallmark = data.frame(msig__H$gs_name, msig__H$gene_symbol)
colnames(msig_WP_WikiPathways) <- c("pathway", "gene")

POI = data.frame(msig_WP_WikiPathways$pathway[msig_WP_WikiPathways$pathway %in% c(pathwayofinterest)], msig_WP_WikiPathways$gene[msig_WP_WikiPathways$pathway %in% c(pathwayofinterest)])
colnames(POI) <- c("pathway", "gene")
GOI <- unique(POI$gene)



#===== use bioMart to convert human symbols to mouse symbols
human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
GOI_df = getLDS(attributes = c("hgnc_symbol"), 
                filters = "hgnc_symbol", 
                values = GOI , 
                mart = human, 
                attributesL = c("mgi_symbol"),  #need to change to mouse
                martL = mouse, 
                uniqueRows=T)
GOI_mouse <- GOI_df$MGI.symbol




#======================== score this gene set in the IRI dataset ===============================
#load Seurat object and convert to SingleCellExpxeriment object
PTtrajectory_A <- readRDS("PTtrajectory_A.rds")
dat <- PTtrajectory_A@assays$RNA@data #export normalized gene expression (not-scaled) gene by cell matrix



#subset to GOI_mouse
dat <- dat[which(rownames(dat)%in% GOI_mouse),] #subset on lift-over rat genes
df <- as.data.frame(dat) #make df

#sort columns by exp.time cell barcodes
df_sorting <- data.frame(barcode=names(PTtrajectory_A$exp.time),
                         exp.time=PTtrajectory_A$exp.time) #export df with barcodes
df_sorting <- df_sorting[order(df_sorting$exp.time),] #sort by exp.time factor levels
df <- df[,df_sorting$barcode]

#score
reads_single_phase = df
reads_single_phase = as.matrix(reads_single_phase[rownames(reads_single_phase) %in% (GOI_mouse) ,])
combined_matrix = rbind(reads_single_phase,average=apply(reads_single_phase,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
reads_single_phase_restricted = reads_single_phase[rownames(reads_single_phase) %in% names(cor_vector[cor_vector >= 0.1]),]
GO_score = apply(reads_single_phase_restricted,2,mean))



#======================== visualize ===============================
#===== Violinplot
PTtrajectory_A$GO_score <- GO_score
PTtrajectory_A$exp.time <- plyr::mapvalues(
  x = PTtrajectory_A$exp.time,
  from = c("Control_Control", 
           "ULIRI_short_1", "ULIRI_short_3", "ULIRI_short_14", 
           "ULIRI_long_1", "ULIRI_long_3", "ULIRI_long_14"),
  to = c("Control", 
         "IRI_short_1", "IRI_short_3", "IRI_short_14",
         "IRI_long_1", "IRI_long_3", "IRI_long_14")
)
table(PTtrajectory_A$exp.time)



#plot
library(vioplot)
makeTransparent = function(..., alpha=0.5) {
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}
col_palette_short <- c(
  "gray", #Control
  "lightblue", #IRI_short_1
  "blue", #IRI_short_3
  "darkblue", #IRI_short_14
  "coral", #IRI_long_1
  "red", #IRI_long_3
  "darkred" #IRI_long_14
)
col_palette_short_trans = makeTransparent(col_palette_short, alpha = 0.3)

pdf(paste0(pathwayofinterest,'_score_by_exp.time_AXISTITLE.pdf'), width = 3.5, height=4)
facts = c(rep(1, table(PTtrajectory_A$exp.time)[1]), 
          rep(2, table(PTtrajectory_A$exp.time)[2]), 
          rep(3, table(PTtrajectory_A$exp.time)[3]), 
          rep(4, table(PTtrajectory_A$exp.time)[4]), 
          rep(5, table(PTtrajectory_A$exp.time)[5]), 
          rep(6, table(PTtrajectory_A$exp.time)[6]), 
          rep(7, table(PTtrajectory_A$exp.time)[7]))
par(mfrow = c(1,1),
    mai=c(1.52,0.82,0.42,0.42)) #bottom, left, top and right in inches
plot(#facts + rnorm(length(facts),0,0.08), 
  GO_score, 
  pch = NA, 
  #col = makeTransparent(col_palette_short[facts], 0.8), 
  #cex = 1.5,  
  main = pathwayofinterest, 
  xaxt = "n", 
  xlab = "", 
  ylab = "Score",
  xlim = c(0.5,7.5)
axis(side=1, #side 1 below, 2 left, 3 above, 4 right 
     las=2, #rotate perpendicular to axis
     at=1:length(table(PTtrajectory_A$exp.time)), names(table(PTtrajectory_A$exp.time)))
vioplot(GO_score ~ facts, plotCentre="line", col = col_palette_short_trans, border = col_palette_short, lwd = 3, add = TRUE)
dev.off()





#===== Dimplot in DiffMap space
load('rd2')
load('crv2')

library(viridis)
nc <- 1
GOscore_df <- data.frame(PTtrajectory_A$GO_score)
GOscore_mt <- matrix(PTtrajectory_A$GO_score)
rownames(GOscore_mt) <- rownames(GOscore_df)
colnames(GOscore_mt) <- paste0(pathwayofinterest," Score")
nms <- paste0(pathwayofinterest," Score")
nr <- ceiling(length(nms)/nc)
pal <- rev(inferno(100, end = 1))
pdf(paste0(pathwayofinterest,'_Score_inferno.pdf'), width=nc*4.5, height=5)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(GOscore_mt[,i], breaks = 100)]
  plot(rd2, col = colors, pch = 16, cex = 0.5, main = i)
  lines(crv2, lwd = 3, col = 'black')
}
dev.off()

