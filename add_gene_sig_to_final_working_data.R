##############################################
# to add gene signatures to final.working.data_v3.with.survival.txt
# Jayesh Kumar Tiwari
# 18 August 2021
###############################################

library(tidyverse)

# creating vectors of the gene signatures 
# (raw.data is test dataframe of any cancer type from TCGA just to make sure that all genes are present or not)
# E.g.- raw.data <- read.table("../firebrowse-data/HNSC/HNSC.uncv2.mRNAseq_RSEM_Z_Score.txt", header = T, sep = '\t')
# 1. 6 gene IFNG signature
ifng.6gene <- c('IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA', 'STAT1', 'IFNG')
# check if all genes are present
intersect(raw.data$symbol, ifng.6gene)  # all present

# 2. activated stroma
activated.stroma <- c('SPARC', 'COL1A2', 'COL3A1', 'POSTN', 'COL5A2', 'COL1A1', 'THBS2', 'FN1', 'COL10A1', 'COL5A1', 'SFRP2', 'CDH11', 'CTHRC1', 'FNDC1', 'SULF1', 'FAP', 'LUM', 'COL11A1', 'ITGA11', 'MMP11', 'INHBA', 'VCAN', 'GREM1', 'COMP')
# check if all genes are present
length(intersect(raw.data$symbol, activated.stroma))  # all present

# 3. stromal enruichment score
stromal.escore <- read.table("deletethis.txt")$V1
stromal.escore[1] <- "DCN"
# check if all genes are present
length(intersect(raw.data$symbol, stromal.escore))  # all present

# 4. macrophages
macrophages <- c('FN1', 'MSR1', 'CD68', 'CCL7', 'PPBP', 'CXCL5')
# check if all genes are present
length(intersect(raw.data$symbol, macrophages))  # all present
# 
# # 5. M1 macrophages
# 
# # 6. m2 macrophages
# macrophages.m2 <- c('CLECSF12', 'IL1Ra', 'Ym1', 'Chi3l3')
# # check if all genes are present
# length(intersect(raw.data$symbol, macrophages.m2))  # all present

# 7. 3 Gene exhausted CD8
cd8.3gene <- c('LAG3', 'CD244', 'EOMES')
# check if all genes are present
length(intersect(raw.data$symbol, cd8.3gene))  # all present

# 8. 9 Gene exhausted CD8
cd8.9gene <- c('CXCL3', 'LAG3', 'CCL5', 'CD244', 'CSF3R', 'CXCL13', 'CYBB', 'KLRK1', 'MSR1')
# check if all genes are present
length(intersect(raw.data$symbol, cd8.9gene))  # all present

# 9. Gajeski 13 gene inflamatory
gajeski <- c('CXCL9', 'CD8A', 'CXCL10', 'CCL2', 'CCL3', 'CCL4', 'GZMK', 'HLA-DMA', 'HLA-DOA', 'HLA-DMB', 'HLA-DOB', 'ICOS', 'IRF1')
# check if all genes are present
length(intersect(raw.data$symbol, gajeski))  # all present

# 10. Cytolytic activity
cytolytic <- c('GZMA', 'PRF1')
# check if all genes are present
length(intersect(raw.data$symbol, cytolytic))  # all present

# 11. Hypoxia
hypoxia <- c("VEGFA","PGAM1","ENO1","LDHA","TPI1","P4HA1","MRPS17","ADM","NDRG1","TUBB6","ALDOA","MIF","SLC2A1","CDKN3","ACOT7")
# check if all genes are present
length(intersect(raw.data$symbol, hypoxia))  # 1 out of 15 Absent
# which not present?
setdiff(hypoxia,intersect(raw.data$symbol, hypoxia))  # VEGF (FIXED ABOVE)

# 12. Granzyme A (GZMA)
# 13. Granzyme B (GZMB)
# 14. Perforin PRF1
length(intersect(raw.data$symbol, c('GZMA','GZMB','PRF1')))   # all present

# integrating all genes (length 222)
all_genes <- c(ifng.6gene, activated.stroma, stromal.escore, macrophages, 
               cd8.3gene, cd8.9gene, gajeski, cytolytic, hypoxia, "GZMA",
               "GZMB","PRF1")

# create an empty dataframe to store the final data
final.data <- data.frame()

# load existing working data
data <- read.table("./txt/final.working.data_v3.with.survival.txt", sep = '\t', header = T)

# make a list of cancers to be run on
cancers <- unique(data$ctype)

for (i in cancers){
  # for testing
  # i = "GBM"
  raw.data <- read.table(paste("../firebrowse-data/",i,"/",i,".uncv2.mRNAseq_RSEM_Z_Score.txt", sep = ""), header = T, sep = '\t')
  
  # cleaning and formatting data
  foo <- data.frame(do.call('rbind', strsplit(as.character(raw.data$gene), '|', fixed = T)))
  raw.data$symbol <- foo$X1
  raw.data$entrez_id <- foo$X2
  raw.data <- raw.data[,-1]
  raw.data <- raw.data %>% select(contains(".01"))
  raw.data$symbol <- foo$X1
  raw.data$entrez_id <- foo$X2
  raw.data <- raw.data %>% select(symbol, entrez_id, everything())
  
  # get only required genes
  required_data <- raw.data %>% filter(symbol %in% all_genes)
  required_data <- data.frame(t(required_data))
  colnames(required_data) <- required_data[1,]
  required_data <- required_data[c(-1,-2),]
  required_data <- mutate_all(required_data, function(x) as.numeric(as.character(x)))
  
  # add gene signature to required_data
  required_data$ifng6gene.sig <- rowMeans(required_data[,which(colnames(required_data) %in% ifng.6gene)])
  required_data$activated.stroma.sig <- rowMeans(required_data[,which(colnames(required_data) %in% activated.stroma)])
  required_data$stromal.escore.sig <- rowMeans(required_data[,which(colnames(required_data) %in% stromal.escore)])
  required_data$macrophages.sig <- rowMeans(required_data[,which(colnames(required_data) %in% macrophages)])
  required_data$cd8.3gene.sig <- rowMeans(required_data[,which(colnames(required_data) %in% cd8.3gene)])
  required_data$cd8.9gene.sig <- rowMeans(required_data[,which(colnames(required_data) %in% cd8.9gene)])
  required_data$gajeski.sig <- rowMeans(required_data[,which(colnames(required_data) %in% gajeski)])
  required_data$cytolytic.sig <- rowMeans(required_data[,which(colnames(required_data) %in% cytolytic)])
  required_data$hypoxia.sig <- rowMeans(required_data[,which(colnames(required_data) %in% hypoxia)])
  required_data$GZMA.sig <- rowMeans(required_data[,which(colnames(required_data) %in% "GZMA"), drop=FALSE])
  required_data$GZMB <- rowMeans(required_data[,which(colnames(required_data) %in% "GZMB"), drop=FALSE])
  required_data$PRF1 <- rowMeans(required_data[,which(colnames(required_data) %in% "PRF1"), drop=FALSE])
  
  # bind all signatures to the required dataframe
  required_data$Sample.ID <- rownames(required_data)
  final.data <- bind_rows(final.data, merge(data %>% filter(ctype == i), required_data, by = "Sample.ID"))
  }

# write to a file
write.table(final.data, file = "final.working.data_v4.txt", row.names = F, sep = '\t')