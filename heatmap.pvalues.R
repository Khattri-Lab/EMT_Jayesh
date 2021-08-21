# #######################################
# # For plotting the heatmaps of pvalues from final.working.data_v3.txt
# # Jayesh Kumar Tiwari
# # 18 June 2021
# ########################################

library(tidyverse)

data <- read.table("final.working.data_v3.txt", header = T, sep = "\t")
wdata <- data %>% filter(class != "intermediate")

stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}

# p-values table of checkpoints--------------------
# base case
spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = "C10orf54")
spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
pvalues.transposed$V1 -> pvalues.only
colnames(pvalues.transposed) <- c("C10orf54")
pvalues.transposed$ctype <- rownames(pvalues.transposed)
pvalues.transposed -> pva
stars.pval(pvalues.only) -> pva.dot
pva[,"C10orf54"] <- pva.dot



for (i in c("CTLA4", "CYBB", "HAVCR2", "ICOS", "LAG3", "PDCD1LG2", "PDCD1", "SIGLEC5", "SIGLEC7", "TNFRSF18", "TNFRSF4", "TNFRSF9", "CD274", "FASLG", "KIR3DL1", "SIGLEC15")) {
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
pva <- pva[,-1]

# getting summary of pva gene wise #######################################
# making empty starter dataframe
stats.pval <- data.frame('gene','high','middle','low')
colnames(stats.pval) <- stats.pval[1,]
isStar3 <- function(X) 
{ X[ ifelse(X=="***", TRUE,FALSE)]
}
isStar2 <- function(X) 
{ X[ ifelse(X=="**", TRUE,FALSE)]
}
isStar1 <- function(X) 
{ X[ ifelse(X=="*", TRUE,FALSE)]
}
for (i in colnames(pva)) {
  print(i)
  print(length(isStar3(pva[,i])))
  stats.pval <- rbind(stats.pval, c(i, length(isStar3(pva[,i])), length(isStar2(pva[,i])), length(isStar1(pva[,i])) ))
}
stats.pval$genetype <- "checkpoint"
# getting summary of pva gene wise ##### ENDED ##################################


pva$ctype <- rownames(pva)

pva.gathered <- pva %>% gather(key = "gene", value = "stars", C10orf54:SIGLEC15)

pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("C10orf54", "VISTA", x)
                }))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CYBB", "NOX2", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("HAVCR2", "TIM3", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("PDCD1LG2", "PD-L2", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("PDCD1", "PD-1", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("TNFRSF18", "GITR", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD274", "PD-L1", x)
}))


checkpoints.plot.pval <- ggplot() +
  geom_tile(data = pva.gathered, aes(x=factor(ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC")), 
                                     y=factor(gene, levels = rev(c("VISTA", "CTLA4", "NOX2", "TIM3", "ICOS", "LAG3", "PD-L2", "PD-1", "SIGLEC5", "SIGLEC7", "GITR", "TNFRSF4", "TNFRSF9", "PD-L1", "FASLG", "KIR3DL1", "SIGLEC15"))) , fill=stars),color="white", size = 1) +
  scale_fill_manual(values = c("light grey","dark grey","black","white"))+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 90, vjust = 0.5, hjust=0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0, vjust = 0.5, hjust=1), axis.title = element_blank())+
  scale_x_discrete(position = "top")
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(plot = checkpoints.plot.pval, filename = "checkpoints.pvalue.final.png", width=221.45625, height=164.57083333, units = "mm", dpi = 300)

# cytokines ---------------------
# base case
spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = "IL10")
spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
pvalues.transposed$V1 -> pvalues.only
colnames(pvalues.transposed) <- c("IL10")
pvalues.transposed$ctype <- rownames(pvalues.transposed)
pvalues.transposed -> pva
stars.pval(pvalues.only) -> pva.dot
pva[,"IL10"] <- pva.dot



for (i in c("TGFB1", "IFNA1", "IFNB1", "IFNG.y", "IL12A", "IL12B", "IL1A", "IL1B", "IL2", "IL3", "IL4", "IL5", "IL6", "IL8", "STAT6", "TNF")) {
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
pva <- pva[,-1]

# getting summary of pva gene wise #######################################
# appending rows to stats.pval dataframe
for (i in colnames(pva)) {
  print(i)
  print(length(isStar3(pva[,i])))
  stats.pval <- rbind(stats.pval, c(i, length(isStar3(pva[,i])), length(isStar2(pva[,i])), length(isStar1(pva[,i])), "cytokine" ))
}

# getting summary of pva gene wise ##### ENDED ##################################


pva$ctype <- rownames(pva)

pva.gathered <- pva %>% gather(key = "gene", value = "stars", IL10:TNF)

pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("IFNG.y", "IFNG", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("TNF", "TNFA", x)
}))


cytokines.plot.pval <- ggplot() +
  geom_tile(data = pva.gathered, aes(x=factor(ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC")), 
                                     y=factor(gene, levels = rev(c("IL10", "TGFB1", "IFNA1", "IFNB1", "IFNG", "IL12A", "IL12B", "IL1A", "IL1B", "IL2", "IL3", "IL4", "IL5", "IL6", "IL8", "STAT6", "TNFA"))) , fill=stars),color="white", size = 1) +
  scale_fill_manual(values = c("light grey","dark grey","black","white"))+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 90, vjust = 0.5, hjust=0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0, vjust = 0.5, hjust=1), axis.title = element_blank())+
  scale_x_discrete(position = "top")
#theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(plot = cytokines.plot.pval, filename = "cytokines.pvalue.final.png", width=221.45625, height=164.57083333, units = "mm", dpi = 300)

# xcell--------------
# base case
spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = "Th1.cells")
spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
pvalues.transposed$V1 -> pvalues.only
colnames(pvalues.transposed) <- c("Th1.cells")
pvalues.transposed$ctype <- rownames(pvalues.transposed)
pvalues.transposed -> pva
stars.pval(pvalues.only) -> pva.dot
pva[,"Th1.cells"] <- pva.dot



for (i in c("Monocytes", "Macrophages.M1", "Macrophages", "B.cells", "CD4..memory.T.cells", "CD4..naive.T.cells", "CD4..T.cells", "CD4..Tcm", "CD4..Tem", "CD8..T.cells", "CD8..naive.T.cells", "CD8..Tcm", "Macrophages.M2", "Memory.B.cells", "Neutrophils", "NK.cells", "Th2.cells","Tregs")) {
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
pva <- pva[,-1]

# getting summary of pva gene wise #######################################
# appending rows to stats.pval dataframe
for (i in colnames(pva)) {
  print(i)
  print(length(isStar3(pva[,i])))
  stats.pval <- rbind(stats.pval, c(i, length(isStar3(pva[,i])), length(isStar2(pva[,i])), length(isStar1(pva[,i])), "xcell" ))
}
stats.pval <- stats.pval[-1,]
stats.pval$high <- as.numeric(as.character(stats.pval$high))
stats.pval$middle <- as.numeric(as.character(stats.pval$middle))
stats.pval$low <- as.numeric(as.character(stats.pval$low))

stats.pval$significant <- stats.pval$high+stats.pval$middle+stats.pval$low

# getting summary of pva gene wise ##### ENDED ##################################


pva$ctype <- rownames(pva)

# reference
# "Th1 cells", "Monocytes", "Macrophages M1", "Macrophages", "B cells", "CD4 memory T cells", "CD4 naive T cells", "CD4 T cells", "CD4 Tcm", "CD4 Tem", "CD8 T cells", "CD8 naive T cells", "CD8 Tcm", "Macrophages M2", "Memory B cells", "Neutrophils", "NK cells", "Th2 cells", "Tregs"
pva.gathered <- pva %>% gather(key = "gene", value = "stars", Th1.cells:Tregs)

pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("Th1.cells", "Th1 cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("Macrophages.M1", "Macrophages M1", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("B.cells", "B cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD4..memory.T.cells", "CD4 memory T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD4..naive.T.cells", "CD4 naive T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD4..T.cells", "CD4 T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD4..Tcm", "Central Memory CD4+ T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD4..Tem", "Effector Memory CD4+ T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD8..T.cells", "CD8 T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD8..naive.T.cells", "CD8 naive T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("CD8..Tcm", "Central memory CD8+ T cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("Macrophages.M2", "Macrophages M2", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("Memory.B.cells", "Memory B cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("NK.cells", "NK cells", x)
}))
pva.gathered <- data.frame(lapply(pva.gathered, function(x) {gsub("Th2.cells", "Th2 cells", x)
}))



xcell.plot.pval <- ggplot() +
  geom_tile(data = pva.gathered, aes(x=factor(ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC")), 
                                     y=factor(gene, levels = rev(c("Th1 cells", "Monocytes", "Macrophages M1", "Macrophages", "B cells", "CD4 memory T cells", "CD4 naive T cells", "CD4 T cells", "Central Memory CD4+ T cells", "Effector Memory CD4+ T cells", "CD8 T cells", "CD8 naive T cells", "Central memory CD8 T cells", "Macrophages M2", "Memory B cells", "Neutrophils", "NK cells", "Th2 cells", "Tregs"))) , fill=stars),color="white", size = 1) +
  scale_fill_manual(values = c("light grey","dark grey","black","white"))+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 90, vjust = 0.5, hjust=0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0, vjust = 0.5, hjust=1), axis.title = element_blank())+
  scale_x_discrete(position = "top")
#theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
xcell.plot.pval

ggsave(plot = xcell.plot.pval, filename = "xcell.pvalue.final.png", width=221.45625, height=164.57083333, units = "mm", dpi = 300)
