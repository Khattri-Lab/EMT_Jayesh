########################
# for preparing new master data of individual genes' expression in all cancer types
# The samples have been be classified in 4 quantiles, the 1st and 4th will be EMT-high and EMT-low
# The data of gene expression is merged with the xCell's data
# It is mesen-epi. High score ones are "mesen" and low score ones are "epi"
# Jayesh Kumar Tiwari
# 13 June 2021 (completed on 18 June 2021)
#########################

# loading libraries----------
library(tidyverse)

# master.data <- data.frame() # ran this is version with classifier_output_v2
master.data_v3 <- data.frame()
# load classifier data of all cancers to get their sample id and their category
#for (i in c(BLCA BRCA CESC COADREAD ESCA FPPP GBM GBMLGG HNSC KIPAN KIRC KIRP LGG LIHC LUAD LUSC OV PAAD PCPG PRAD READ SARC SKCM STAD STES TGCT THCA THYM UCEC))
for (i in c("BRCA" ,"CESC", "COADREAD", "ESCA","GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", 'LUAD', "LUSC", "OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "UCEC")){
  # # do it once with BLCA alone to make a master.data_v3() bcoz rbind can't be used on empty dataframe
  # i = "BLCA" # for testing
  cat("Starting ",i, "... \n")
  cldata <- read.table(paste("./classifier_output_v3/",i,"_classification.txt",sep = ""), header = T, sep = '\t')
  mdata <- read.table(paste("../firebrowse-data/gdac.broadinstitute.org_",i,".mRNAseq_Preprocess.Level_3.2016012800.0.0/",i,".uncv2.mRNAseq_RSEM_Z_Score.txt",sep = ""), header = T, sep = '\t')
  
  # cleaning and formatting data of the raw files
  foo <- data.frame(do.call('rbind', strsplit(as.character(mdata$gene), '|', fixed = T)))
  mdata <- mdata %>% select(contains(".01"))
  mdata$symbol <- foo$X1
  mdata$entrez_id <- foo$X2
  # mdata <- mdata[,-1] # this was for old one when there used to be an symbol+gene column
  
  # get the checkpoint, cytokines and t-cell markers from raw data
  checkpoint <- read.table("immunecheckpoint.txt", header = F, sep = '\t')
  inflam <- read.table("immuno-inflamatory.txt", header = F, sep = '\t')
  cd8 <- read.table("cd8.txt", header = F, sep = '\t')
  data <- mdata %>% select(symbol, entrez_id, everything())
  
  # creating vectors of genes
  cd8=intersect(cd8$V1,data[,1])
  checkpoint=intersect(checkpoint$V1,data[,1])
  inflam=intersect(inflam$V1,data[,1])
  
  # preparing dataframes of individual gene sets form master data
  data[c(which(data[,1] %in% cd8)),] -> d_cd8
  rownames(d_cd8) <- d_cd8[,1]
  d_cd8 <- d_cd8[,-c(1,2)]
  d_cd8 <- data.frame(t(d_cd8))
  d_cd8$Sample.ID <- rownames(d_cd8)
  data[c(which(data[,1] %in% checkpoint)),] -> d_checkpoint
  rownames(d_checkpoint) <- d_checkpoint[,1]
  d_checkpoint <- d_checkpoint[,-c(1,2)]
  d_checkpoint <- data.frame(t(d_checkpoint))
  d_checkpoint$Sample.ID <- rownames(d_checkpoint)
  data[c(which(data[,1] %in% inflam)),] -> d_inflam
  rownames(d_inflam) <- d_inflam[,1]
  d_inflam <- d_inflam[,-c(1,2)]
  d_inflam <- data.frame(t(d_inflam))
  d_inflam$Sample.ID <- rownames(d_inflam)
  
  # binding all genes together
  fdata <- as.data.frame(list(d_cd8, d_checkpoint, d_inflam) %>% reduce(full_join, by = "Sample.ID"))
  cat("For",i, " Dimensions of all genes is: ",dim(fdata))
  cat("\n For",i, " Dimensions of category is: ",dim(cldata))
  cat("\n For",i, " Dimensions of cd8, checkpoint and inflam genes is: ", dim(d_cd8), " and ", dim(d_checkpoint), " and ", dim(d_inflam))
  ffdata <- merge(fdata, cldata, by.x = "Sample.ID", by.y = "name")
  ffdata$ctype = i
  cat("\n Dimensions of final data are: ", dim(ffdata))
  # master.data <- rbind(master.data, fdata)  # already ran for classifier_output_v2
  # IN FIRST CANCER we will make the master.data_v3, THEN RUN for loop for FURTHER ITERATIONS
  # master.data_v3 <- ffdata # RUN THIS ONLY ON FIRST CANCER
  master.data_v3 <- rbind(master.data_v3, ffdata) # RUN THIS IN FOR LOOP
  cat("\nDim of master data version 2 are", dim(master.data_v3), "after ",i,"\n")
}
master.data_v3 <- master.data_v3%>% select(Sample.ID, score, class, ctype, everything())
write.table(master.data_v2, file = "master.data_v4.txt", row.names = F, sep = "\t")

# adding xcell data
# working.data <- merge(master.data, xcell.data, by.x = "Sample.ID", by.y = "Sample.ID")
working.data_v3 <- merge(master.data_v3, xcell.data, by.x = "Sample.ID", by.y = "Sample.ID")
working.data_v3$class <- gsub('epi','chymal',working.data_v3$class)
working.data_v3$class <- gsub('mesen','thelial',working.data_v3$class)
working.data_v3$class <- gsub('chymal','mesenchymal',working.data_v3$class)
working.data_v3$class <- gsub('thelial','epithelial',working.data_v3$class)


write.table(working.data_v3, file = "final.working.data_v3.txt", row.names = F, sep = "\t")

# generate stats of the working data ---------------------
epivec <- c()
mesenvec <- c()
intvec <- c()
pepi <- c()
pmesen <- c()
pint <- c()
cancer <- c()
for (i in c("BLCA", "BRCA" ,"CESC", "COADREAD", "ESCA","GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", 'LUAD', "LUSC", "OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "TGCT", "UCEC")){
  epithelial = length(((working.data_v3 %>% filter(ctype == i & class == "epithelial")))$Sample.ID)
  mesenchymal = length(((working.data_v3 %>% filter(ctype == i & class == "mesenchymal")))$Sample.ID)
  intermediate = length(((working.data_v3 %>% filter(ctype == i & class == "intermediate")))$Sample.ID)
  percent_epi = (epithelial/(epithelial+mesenchymal+intermediate))*100
  percent_mesen = (mesenchymal/(epithelial+mesenchymal+intermediate))*100
  percent_intmdt = (intermediate/(epithelial+mesenchymal+intermediate))*100
  epivec <- c(epivec, epithelial)
  mesenvec <- c(mesenvec, mesenchymal)
  intvec <- c(intvec, intermediate)
  pepi <- c(pepi, percent_epi)
  pmesen <- c(pmesen, percent_mesen)
  pint <- c(pint, percent_intmdt)
  cancer <- c(cancer, i)
}
data_stats_working_v6 <- data.frame(cancer.type = cancer, epithelial = epivec, mesenchymal = mesenvec, intermediate = intvec, epi_percent = pepi, mesen_percent = pmesen, intmdt_percent = pint)

write.table(data_stats_working_v6, file = "stats.final.working.data_v3.txt", row.names = F, sep = "\t")

# testing coadread data
classifier_coadread <- master.data_v2 %>% filter(ctype == "COADREAD")
xcelled_coadread <- working.data_v2 %>% filter(ctype == "COADREAD")
raw_classifier_coadread <- read.table("./classifier_output_v3/COADREAD_classification.txt", header = T)
