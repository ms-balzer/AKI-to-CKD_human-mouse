library(biomaRt)
set.seed(123)
setwd('/.../1__1to1mapping/')



#load mouse DEG lists
mouse_balzer <- read.csv("Balzer__mouse.csv", header=TRUE)
mouse_doke <- read.csv("Doke__mouse.csv", header=TRUE)
mouse_kirita <- read.csv("Kirita__mouse.csv", header=TRUE)



#see different datasets available within a biomaRt:
mart = useMart('ensembl', host = "https://dec2021.archive.ensembl.org")
#listDatasets(mart)
# 79           hsapiens_gene_ensembl                                     Human genes (GRCh38.p13)                        GRCh38.p13
# 106         mmusculus_gene_ensembl                                         Mouse genes (GRCm39)                            GRCm39


human = useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")
mouse = useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl", verbose = TRUE, host = "https://dec2021.archive.ensembl.org")

# create x_genelist:
mouse_balzer_genelist = mouse_balzer$gene
mouse_doke_genelist = mouse_doke$gene
mouse_kirita_genelist = mouse_kirita$gene

#listAttributes(human)
# 63                                       hgnc_symbol                                                                  HGNC symbol
#listAttributes(mouse)
# 60                                        mgi_symbol                                                                   MGI symbol

#---getLDS: repeat code below for other groups as well----
# getLDS(): function that links 2 datasets and retrieves information from these linked BioMart datasets
genesV2_mouse_balzer = getLDS(attributes = c("mgi_symbol"), 
                              filters = "mgi_symbol", 
                              values = mouse_balzer_genelist, 
                              mart = mouse, 
                              attributesL = c("hgnc_symbol"),
                              martL = human, 
                              uniqueRows=T)

genesV2_mouse_doke = getLDS(attributes = c("mgi_symbol"), 
                            filters = "mgi_symbol", 
                            values = mouse_doke_genelist, 
                            mart = mouse, 
                            attributesL = c("hgnc_symbol"),
                            martL = human, 
                            uniqueRows=T)

genesV2_mouse_kirita = getLDS(attributes = c("mgi_symbol"), 
                              filters = "mgi_symbol", 
                              values = mouse_kirita_genelist, 
                              mart = mouse, 
                              attributesL = c("hgnc_symbol"),
                              martL = human, 
                              uniqueRows=T)

# Create a new dataframe which has only the genes which were successfully mapped
mouse_balzer= mouse_balzer[mouse_balzer$gene %in% genesV2_mouse_balzer$MGI.symbol,]
mouse_doke= mouse_doke[mouse_doke$gene %in% genesV2_mouse_doke$MGI.symbol,]
mouse_kirita= mouse_kirita[mouse_kirita$gene %in% genesV2_mouse_kirita$MGI.symbol,]
names(genesV2_mouse_balzer)[1] <- "gene"
names(genesV2_mouse_doke)[1] <- "gene"
names(genesV2_mouse_kirita)[1] <- "gene"
genesV2_mouse_balzer = genesV2_mouse_balzer[!duplicated(genesV2_mouse_balzer[c("gene")]),]
genesV2_mouse_doke = genesV2_mouse_doke[!duplicated(genesV2_mouse_doke[c("gene")]),]
genesV2_mouse_kirita = genesV2_mouse_kirita[!duplicated(genesV2_mouse_kirita[c("gene")]),]

# merge gene list with mouse dataset to create the new dataset with human gene names
mouse_balzer <- merge(mouse_balzer, genesV2_mouse_balzer, by="gene", sort = F)
mouse_doke <- merge(mouse_doke, genesV2_mouse_doke, by="gene", sort = F)
mouse_kirita <- merge(mouse_kirita, genesV2_mouse_kirita, by="gene", sort = F)


write.csv(mouse_balzer, "../0__DEGs/Balzer__mouse__human_genenames.csv", row.names = F)
write.csv(mouse_balzer, "../0__DEGs/Balzer__mouse__human_genenames.csv", row.names = F)
write.csv(mouse_doke, "../0__DEGs/Doke__mouse__human_genenames.csv", row.names = F)
write.csv(mouse_kirita, "../0__DEGs/Kirita__mouse__human_genenames.csv", row.names = F)
