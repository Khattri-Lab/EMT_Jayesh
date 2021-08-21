#####################################################################
# to plot corrplots and heatmap of p-values (asterisks) of gene signatures
# Jayesh Kumar Tiwari
# 20 August 2021
#########################################################################

library(data.table)
library(corrplot)

# load data
final.data <- read.table("./txt/final.working.data_v4.txt", header = T, sep = '\t')

# filter useful columns
wdata <- final.data %>% select("Sample.ID":"PFI.time", "ifng6gene.sig":"GZMA.sig", "GZMB", "PRF1.y")

# create a dataframe of difference of medians of samples for each cancer type
# i.e for each (gene sig, cancer type) pair, it's value will be 
# median(expression of all samples which are "mesenchymal") - median(expression of all samples which are "epithelial")
reqdata <- wdata %>%
  group_by(ctype) %>%
  summarize(across(ifng6gene.sig:PRF1.y, ~ median(.[class=="mesenchymal"] - median(.[class=="epithelial"]))))
# converting tibble to data.frame
reqdata <- data.frame(reqdata)
rownames(reqdata) <- reqdata$ctype
reqdata <- reqdata[,-1]
reqdata <- data.frame(t(reqdata))
# setting order of cancers as it will be used as it is in corrplot
reqdata <- reqdata %>% select(SKCM, COADREAD, UCEC, SARC, LIHC, LUSC, PRAD, OV, PCPG, GBM, ESCA, BRCA, CESC, TGCT,  LUAD,  PAAD,  LGG, BLCA, KIRP,  HNSC,  STAD,  KIRC)

# for k means
# check for NA
sum(sapply(reqdata, is.na))
# check for infinite
sum(sapply(reqdata, is.infinite))
# check for NaN
sum(sapply(reqdata, is.nan))

# calcuting optimal number of clusters by "Silhouette method"
library(factoextra)
fviz_nbclust(reqdata, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")+
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
FOOk2 <- kmeans(reqdata, centers = 2, nstart = 25)
# we see the distrubution of the data in the 2 clusters
fviz_cluster(FOOk2, data = reqdata, frame.type = "ellipse") + theme_minimal() + ggtitle("k = 2") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

# this is for KMEANS
# making a dataframe of FOOk$cluster
sample.order.data <- data.frame(FOOk2$cluster)
rownames(sample.order.data %>% filter(FOOk2.cluster==1)) -> cluster1Samples
rownames(sample.order.data %>% filter(FOOk2.cluster==2)) -> cluster2Samples
# clusterSamples2 is the order of gene signature columns arranged according to the cluster they fall in(in this case cluster 1 and cluster 2)
clusteredSamples2 <- c(cluster1Samples,rev(cluster2Samples))

# set order of gene signatures
reqdata$genesig <- rownames(reqdata)
reqdata[match(clusteredSamples2, reqdata$genesig),] -> reqdata


reqdata <- data.frame(t(reqdata))
reqdata <- reqdata %>% select(SKCM, COADREAD, UCEC, SARC, LIHC, LUSC, PRAD, OV, PCPG, GBM, ESCA, BRCA, CESC, TGCT,  LUAD,  PAAD,  LGG, BLCA, KIRP,  HNSC,  STAD,  KIRC)

# saving corrplot
png("genesigs.corrplot.bw.png", width=362.47916667, height=170.39166667, units = "mm", res = 118.11)

corrplot(as.matrix(reqdata), is.corr = FALSE, method = "square", bg="white", tl.col = "black")

dev.off()

#######################################################################################################

# preparing a dataframe similar to reqdata, but this time the values will be asterisks(*) corresponding to
# *** = very highly significant
# ** = highly significant
# * = significant

# p-values table of checkpoints--------------------
# funtion
stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}

# base case
spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = "ifng6gene.sig")
spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
pvalues.transposed$V1 -> pvalues.only
colnames(pvalues.transposed) <- c("ifng6gene.sig")
pvalues.transposed$ctype <- rownames(pvalues.transposed)
pvalues.transposed -> pva
stars.pval(pvalues.only) -> pva.dot
pva[,"ifng6gene.sig"] <- pva.dot

# run a for loop to add each gene signature
for (i in c("activated.stroma.sig",
            "stromal.escore.sig",   "macrophages.sig",      "cd8.3gene.sig",
            "cd8.9gene.sig",        "gajeski.sig",          "cytolytic.sig",
            "hypoxia.sig",          "GZMA.sig",             "GZMB",
            "PRF1.y")) {
  spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = i)
  spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
  pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
  pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
  pvalues.transposed$V1 -> pvalues.only
  colnames(pvalues.transposed) <- c(i)
  pvalues.transposed$ctype <- rownames(pvalues.transposed)
  stars.pval(pvalues.only) -> pva.dot
  pvalues.transposed[,i] <- pva.dot
  merge(pva, pvalues.transposed, by.x = "ctype", by.y = "ctype")-> pva
}
rownames(pva) <- pva[,1]
# pva is the final column which we will use (it is equivalent to reqdata)
pva <- pva[,-1]

########################################
# for plotting heatmap of pvalues(asterisks)
pva$ctype <- rownames(pva)

pva.gathered <- pva %>% gather(key = "gene", value = "stars", ifng6gene.sig:PRF1.y)

genesigs.plot.pval <- ggplot() +
  geom_tile(data = pva.gathered, aes(x=factor(ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC")), 
                                     y=factor(gene, levels = clusteredSamples2) , fill=stars),color="white", size = 1) +
  scale_fill_manual(values = c("light grey","dark grey","black","white"))+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 90, vjust = 0.5, hjust=0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0, vjust = 0.5, hjust=1), axis.title = element_blank())+
  scale_x_discrete(position = "top")
#theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
