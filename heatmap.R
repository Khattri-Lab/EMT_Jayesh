# #######################################
# # For plotting the boxplots and heatmaps of medians from final.working.data_v3.txt
# # Jayesh Kumar Tiwari
# # 18 June 2021
# ########################################

library(tidyverse)

data <- read.table("/mnt/swift/Jayesh/emt/script.jayesh.edited/txt/final.working.data_v3.txt", header = T, sep = "\t")
wdata <- data %>% filter(class != "intermediate")

# setting the orer of genes increasing in median
crap_vec <- c()
for (i in c("SKCM", "COADREAD", "UCEC", "SARC", "LIHC", "LUSC", "PRAD", "OV", "PCPG", "GBM", "ESCA", "BRCA", "CESC", "TGCT",  "LUAD",  "PAAD",  "LGG", "BLCA", "KIRP",  "HNSC",  "STAD",  "KIRC")) {
       crap_vec <- c(crap_vec, median((wdata %>% filter(ctype==i))$score),"\n")
   }
crap_vec <- sort(crap_vec)
for (j in crap_vec){
  for (i in c("SKCM", "COADREAD", "UCEC", "SARC", "LIHC", "LUSC", "PRAD", "OV", "PCPG", "GBM", "ESCA", "BRCA", "CESC", "TGCT",  "LUAD",  "PAAD",  "LGG", "BLCA", "KIRP",  "HNSC",  "STAD",  "KIRC")) {
    if(j == median((wdata %>% filter(ctype==i))$score)){
      cat(i,"")
    }
  }
}

#functions-------------------
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}
po <- position_dodge(width = 0.8)

for (i in c()) {

}
# BETTER ONE IS BELOW THIS
# plot_CD8..Tem.M2 <- ggplot(wdata, aes(x=ctype, y=as.numeric(CD8..Tem.M2), fill=class)) +
#   stat_boxplot(geom = "errorbar", width = 0.3, position = po) +
#   stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", position = po, width=0.68, lwd=0.7) +
#   scale_fill_grey(start = 1, end=.8) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), axis.title.x = element_text(colour = "black", size = 12), axis.title.y = element_text(colour = "black", size=12), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#   labs(title="CD8..Tem M2", x ="Cancer Type", y = "Expression Level")
# plot_CD8..Tem.M2

plot <- ggplot(wdata, aes(x=ctype, y=as.numeric(CD8..Tem), fill=class)) + 
  stat_boxplot(geom = "errorbar", width = 0.3, position = po) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", position = po, width=0.68, lwd=0.7) +
  #scale_fill_grey(start = .4, end=1) +
  scale_fill_manual(values = c("#A1A1A1","#FFFFFF","#555555")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), axis.title.x = element_text(colour = "black", size = 12), axis.title.y = element_text(colour = "black", size=12), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="CD8..Tem", x ="Cancer Type", y = "Expression Level")

#calculating p values
spreaded.dup.data <- wdata %>% tidyr::spread(key = "ctype", value = "CD8..Tem")
spreaded.dup.data %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
pvalues.transposed$V1 -> pvalues.only

pvalues.transposed -> pva

stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}

stars.pval(pvalues.only) -> pva.dot
plot
#adding asterisks to the plot
plot + annotate("text", x=c(1), y=c(0.1), label=pva.dot[1], size=4.5) -> g
g
for (i in c(2:22)) {
  g + annotate("text", x=c(i), y= 0.1, label=pva.dot[i], size=4.5) -> g
}

#adding ___ to the plot
g + annotate("text", x=c(1), y=c(0.1), label="___", size=4.5) -> g
for (i in c(2:22)) {
  g + annotate("text", x=c(i), y= 0.1, label="___", size=4.5) -> g
}

g
ggsave(plot = g, filename = "plot.CD8..Tem.png", width=345.01, height=189.44, units = "mm", dpi = 300)

# heatmap-------------------
library(data.table)
epi_data <- wdata %>% filter(class == "epithelial")
mesen_data <- wdata %>% filter(class == "mesenchymal")

# summarized epi data making --------------------------
summarized_data_epi <- setDT(epi_data)[,list("CD8A"=as.numeric(median(get("CD8A")))), by=ctype]
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("CD8B"=as.numeric(median(get("CD8B")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IFNG.x"=as.numeric(median(get("IFNG.x")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("PRF1"=as.numeric(median(get("PRF1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("C10orf54"=as.numeric(median(get("C10orf54")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("CD274"=as.numeric(median(get("CD274")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("CTLA4"=as.numeric(median(get("CTLA4")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("CYBB"=as.numeric(median(get("CYBB")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("FASLG"=as.numeric(median(get("FASLG")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("HAVCR2"=as.numeric(median(get("HAVCR2")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("ICOS"=as.numeric(median(get("ICOS")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("KIR3DL1"=as.numeric(median(get("KIR3DL1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("LAG3"=as.numeric(median(get("LAG3")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("PDCD1LG2"=as.numeric(median(get("PDCD1LG2")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("PDCD1"=as.numeric(median(get("PDCD1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("SIGLEC15"=as.numeric(median(get("SIGLEC15")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("SIGLEC5"=as.numeric(median(get("SIGLEC5")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("SIGLEC7"=as.numeric(median(get("SIGLEC7")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("TNFRSF18"=as.numeric(median(get("TNFRSF18")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("TNFRSF4"=as.numeric(median(get("TNFRSF4")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("TNFRSF9"=as.numeric(median(get("TNFRSF9")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IFNA1"=as.numeric(median(get("IFNA1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IFNB1"=as.numeric(median(get("IFNB1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IFNG.y"=as.numeric(median(get("IFNG.y")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL10"=as.numeric(median(get("IL10")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL12A"=as.numeric(median(get("IL12A")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL12B"=as.numeric(median(get("IL12B")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL1A"=as.numeric(median(get("IL1A")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL1B"=as.numeric(median(get("IL1B")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL2"=as.numeric(median(get("IL2")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL3"=as.numeric(median(get("IL3")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL4"=as.numeric(median(get("IL4")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL5"=as.numeric(median(get("IL5")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL6"=as.numeric(median(get("IL6")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("IL8"=as.numeric(median(get("IL8")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("STAT6"=as.numeric(median(get("STAT6")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("TGFB1"=as.numeric(median(get("TGFB1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("TNF"=as.numeric(median(get("TNF")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Tregs"=as.numeric(median(get("Tregs")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Th2.cells"=as.numeric(median(get("Th2.cells")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Th1.cells"=as.numeric(median(get("Th1.cells")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("NK.cells"=as.numeric(median(get("NK.cells")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Neutrophils"=as.numeric(median(get("Neutrophils")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Monocytes"=as.numeric(median(get("Monocytes")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("Memory.B.cells"=as.numeric(median(get("Memory.B.cells")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("Macrophages.M1"=as.numeric(median(get("Macrophages.M1")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("Macrophages.M2"=as.numeric(median(get("Macrophages.M2")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("Macrophages"=as.numeric(median(get("Macrophages")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD8..Tem"=as.numeric(median(get("CD8A")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD8..Tcm"=as.numeric(median(get("CD8..Tcm")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD8..naive.T.cells"=as.numeric(median(get("CD8..naive.T.cells")))), by=ctype])
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD8..T.cells"=as.numeric(median(get("CD8..T.cells")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD4..Tem"=as.numeric(median(get("CD4..Tem")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD4..Tcm"=as.numeric(median(get("CD4..Tcm")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD4..T.cells"=as.numeric(median(get("CD4..T.cells")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD4..naive.T.cells"=as.numeric(median(get("CD4..naive.T.cells")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi, setDT(epi_data)[,list("CD4..memory.T.cells"=as.numeric(median(get("CD4..memory.T.cells")))), by=ctype]) 
summarized_data_epi <-merge(summarized_data_epi,  setDT(epi_data)[,list("B.cells"=as.numeric(median(get("B.cells")))), by=ctype])

# summarized mesen data making ---------------------
summarized_data_mesen <- setDT(mesen_data)[,list("CD8A"=as.numeric(median(get("CD8A")))), by=ctype]
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("CD8B"=as.numeric(median(get("CD8B")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IFNG.x"=as.numeric(median(get("IFNG.x")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("PRF1"=as.numeric(median(get("PRF1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("C10orf54"=as.numeric(median(get("C10orf54")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("CD274"=as.numeric(median(get("CD274")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("CTLA4"=as.numeric(median(get("CTLA4")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("CYBB"=as.numeric(median(get("CYBB")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("FASLG"=as.numeric(median(get("FASLG")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("HAVCR2"=as.numeric(median(get("HAVCR2")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("ICOS"=as.numeric(median(get("ICOS")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("KIR3DL1"=as.numeric(median(get("KIR3DL1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("LAG3"=as.numeric(median(get("LAG3")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("PDCD1LG2"=as.numeric(median(get("PDCD1LG2")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("PDCD1"=as.numeric(median(get("PDCD1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("SIGLEC15"=as.numeric(median(get("SIGLEC15")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("SIGLEC5"=as.numeric(median(get("SIGLEC5")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("SIGLEC7"=as.numeric(median(get("SIGLEC7")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("TNFRSF18"=as.numeric(median(get("TNFRSF18")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("TNFRSF4"=as.numeric(median(get("TNFRSF4")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("TNFRSF9"=as.numeric(median(get("TNFRSF9")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IFNA1"=as.numeric(median(get("IFNA1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IFNB1"=as.numeric(median(get("IFNB1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IFNG.y"=as.numeric(median(get("IFNG.y")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL10"=as.numeric(median(get("IL10")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL12A"=as.numeric(median(get("IL12A")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL12B"=as.numeric(median(get("IL12B")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL1A"=as.numeric(median(get("IL1A")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL1B"=as.numeric(median(get("IL1B")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL2"=as.numeric(median(get("IL2")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL3"=as.numeric(median(get("IL3")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL4"=as.numeric(median(get("IL4")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL5"=as.numeric(median(get("IL5")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL6"=as.numeric(median(get("IL6")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("IL8"=as.numeric(median(get("IL8")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("STAT6"=as.numeric(median(get("STAT6")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("TGFB1"=as.numeric(median(get("TGFB1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("TNF"=as.numeric(median(get("TNF")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Tregs"=as.numeric(median(get("Tregs")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Th2.cells"=as.numeric(median(get("Th2.cells")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Th1.cells"=as.numeric(median(get("Th1.cells")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("NK.cells"=as.numeric(median(get("NK.cells")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Neutrophils"=as.numeric(median(get("Neutrophils")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Monocytes"=as.numeric(median(get("Monocytes")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("Memory.B.cells"=as.numeric(median(get("Memory.B.cells")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("Macrophages.M1"=as.numeric(median(get("Macrophages.M1")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("Macrophages.M2"=as.numeric(median(get("Macrophages.M2")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("Macrophages"=as.numeric(median(get("Macrophages")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD8..Tem"=as.numeric(median(get("CD8A")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD8..Tcm"=as.numeric(median(get("CD8..Tcm")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD8..naive.T.cells"=as.numeric(median(get("CD8..naive.T.cells")))), by=ctype])
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD8..T.cells"=as.numeric(median(get("CD8..T.cells")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD4..Tem"=as.numeric(median(get("CD4..Tem")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD4..Tcm"=as.numeric(median(get("CD4..Tcm")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD4..T.cells"=as.numeric(median(get("CD4..T.cells")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD4..naive.T.cells"=as.numeric(median(get("CD4..naive.T.cells")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen, setDT(mesen_data)[,list("CD4..memory.T.cells"=as.numeric(median(get("CD4..memory.T.cells")))), by=ctype]) 
summarized_data_mesen <-merge(summarized_data_mesen,  setDT(mesen_data)[,list("B.cells"=as.numeric(median(get("B.cells")))), by=ctype])

# subtracting mesen -epi dataframe-----------------
# method 1
fs <- summarized_data_epi
fs[1:nrow(fs), 2:ncol(fs)] <- NA

for (i in colnames(summarized_data_epi)) {
  if (i == "ctype"){
    next
  }
  print(i)
  fs[[i]] <- summarized_data_mesen[[i]] - summarized_data_epi[[i]]
}

gathered <- fs %>% gather(key = "genes_and_cells", value = "expression", CD8A:TNF)
gathered.xcell <- fs %>% gather(key = "genes_and_cells", value = "expression", Tregs:Macrophages,CD8..Tcm:B.cells)
ggplot() + 
  geom_tile(data = gathered, 
            aes(x=ctype, y=genes_and_cells, fill=expression)) +
  scale_fill_gradient2(high = "red",mid = "white", low = "blue") +
  # geom_tile(data = subset(ghdata, expression > 5.0), 
  #           aes(x=barcode, y=Gene, color=Expression)) +
  scale_colour_manual(name = "Expression", values = c('black')) +
  theme(text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot() + 
  geom_tile(data = gathered.xcell, 
            aes(x=ctype, y=genes_and_cells, fill=expression)) +
  scale_fill_gradient2(high = "red",mid = "white", low = "blue") +
  # geom_tile(data = subset(ghdata, expression > 5.0), 
  #           aes(x=barcode, y=Gene, color=Expression)) +
  scale_colour_manual(name = "Expression", values = c('black')) +
  theme(text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#################### EXPERIMENT ############################
# clustering (hclust) the genes separately and cells separately
fs %>% select(CD8A:TNF) -> genes_data
rownames(genes_data) <- fs$ctype
genes_data <- data.frame(t(genes_data))
colnames(genes_data) <- fs$ctype
FOOcl <- genes_data
FOOcl <- FOOcl[,-23]
FOOcl <- mutate_all(FOOcl, function(x) as.numeric(as.character(x)))
FOOcl <- as.matrix(FOOcl)

# for hclust
FOOcl <- dist(FOOcl)
FOOcl <- hclust(FOOcl)
plot(FOOcl)

# for k means
# check for NA
sum(sapply(FOOcl, is.na))
# check for infinite
sum(sapply(FOOcl, is.infinite))
# check for NaN
sum(sapply(FOOcl, is.nan))

# calcuting optimal number of clusters by "Silhouette method"
library(factoextra)
fviz_nbclust(FOOcl, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")+
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
FOOk <- kmeans(FOOcl, centers = 2, nstart = 25)
# we see the distrubution of the data in the 2 clusters
fviz_cluster(FOOk, data = FOOcl, frame.type = "ellipse") + theme_minimal() + ggtitle("k = 2") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

# this is for KMEANS
# making a dataframe of FOOk$cluster
sample.order.data <- data.frame(FOOk$cluster)
rownames(sample.order.data %>% filter(FOOk.cluster==1)) -> cluster1Samples
rownames(sample.order.data %>% filter(FOOk.cluster==2)) -> cluster2Samples
clusteredSamples <- c(cluster1Samples,cluster2Samples)

# below is for HCLUST
crap_vec <- c()
for (variable in FOOcl$order) {
  crap_vec <- c(crap_vec, rownames(genes_data)[variable])
}
rownames(genes_data) -> genes_data$gene_type
genes_data[match(crap_vec, genes_data$gene_type),] -> genedatasorted
genedatasorted <- genedatasorted[,-23]
genedatasorted <- data.frame(t(genedatasorted))
genedatasorted$ctype <- rownames(genedatasorted)
gathered <- genedatasorted %>% gather(key = "genes_and_cells", value = "expression", C10orf54:TNF)
# for hclust
gathered$genes_and_cells <- factor(gathered$genes_and_cells, levels = crap_vec)
# for kmeans
gathered$genes_and_cells <- factor(gathered$genes_and_cells, levels = clusteredSamples)

ggplot() + 
  geom_tile(data = gathered, 
            aes(x=ctype, y=genes_and_cells, fill=expression)) +
  scale_fill_gradient2(high = "red",mid = "white", low = "blue") +
  # geom_tile(data = subset(ghdata, expression > 5.0), 
  #           aes(x=barcode, y=Gene, color=Expression)) +
  scale_colour_manual(name = "Expression", values = c('black')) +
  theme(text = element_text(size = 10))

M <- M[,-39]
M <- data.frame(t(M))
rownames(M) -> M$gene_type
M[match(clusteredSamples, M$gene_type),] -> M
M <- M[,-23]
M <- as.matrix(M)
corrplot(t(M), is.corr = FALSE, method = "square")

# plotting corrplot for CD8 molecules
cd8 <- data.frame(t(M)) %>% select(CD8A, PRF1, IFNG.x, CD8B)
colnames(cd8) <- c("CD8A","PRF1","IFNG","CD8B")
cd8 <- t(cd8)
cd8 <- data.frame(cd8) %>% select(SKCM, COADREAD, UCEC, SARC, LIHC, LUSC, PRAD, OV, PCPG, GBM, ESCA, BRCA, CESC, TGCT,  LUAD,  PAAD,  LGG, BLCA, KIRP,  HNSC,  STAD,  KIRC)
cd8.plot <- corrplot(as.matrix(cd8), is.corr = FALSE, method = "square", tl.col = "black")

# plotting corrplot for checkpoint molecules
checkpoints <- data.frame(t(M)) %>% select(C10orf54, CTLA4, CYBB, HAVCR2, ICOS, LAG3, PDCD1LG2, PDCD1, SIGLEC5, SIGLEC7, TNFRSF18, TNFRSF4, TNFRSF9, CD274, FASLG, KIR3DL1, SIGLEC15)
colnames(checkpoints) <-             c("VISTA", "CTLA4", "NOX2", "TIM3", "ICOS", "LAG3", "PD-L2", "PD-1", "SIGLEC5", "SIGLEC7", "GITR", "TNFRSF4", "TNFRSF9", "PD-L1", "FASLG", "KIR3DL1", "SIGLEC15")
checkpoints <- t(checkpoints)
checkpoints <- data.frame(checkpoints) %>% select(SKCM , PRAD, COADREAD, PAAD, BLCA, UCEC, LGG, STAD, TGCT, LUAD, BRCA, CESC, OV, ESCA, LUSC, GBM, PCPG, LIHC, HNSC, SARC, KIRP, KIRC)
col <- colorRampPalette(c("light grey", "grey", "dark grey", "black"))
corrplot(as.matrix(checkpoints), is.corr = FALSE, method = "shade", col=gray(16:0 / 16), bg="white", tl.col = "black")

png("checkpoints.corrplot.bw.png", width=362.47916667, height=170.39166667, units = "mm", res = 118.11)

corrplot(as.matrix(checkpoints), is.corr = FALSE, method = "square", col=gray(16:0 / 16), bg="white", tl.col = "black")

# corrplot(as.matrix(checkpoints), is.corr = FALSE, method = "square", tl.col = "black")

dev.off()

# plotting corrplot for cytokines
cytokines <- data.frame(t(M)) %>% select(IL10, TGFB1, IFNA1, IFNB1, IFNG.y, IL12A, IL12B, IL1A, IL1B, IL2, IL3, IL4, IL5, IL6, IL8, STAT6, TNF)
colnames(cytokines) <- c("IL10", "TGFB1", "IFNA1", "IFNB1", "IFNG", "IL12A", "IL12B", "IL1A", "IL1B", "IL2", "IL3", "IL4", "IL5", "IL6", "CXCL8", "STAT6", "TNFA")
cytokines <- t(cytokines)
cytokines <- data.frame(cytokines) %>% select(SKCM , PRAD, COADREAD, PAAD, BLCA, UCEC, LGG, STAD, TGCT, LUAD, BRCA, CESC, OV, ESCA, LUSC, GBM, PCPG, LIHC, HNSC, SARC, KIRP, KIRC)


png("cytokines.corrplot.bw.png", width=362.47916667, height=170.39166667, units = "mm", res = 118.11)


corrplot(as.matrix(cytokines), is.corr = FALSE, method = "square",col=gray(16:0 / 16), bg="white", tl.col = "black")
# corrplot(as.matrix(cytokines), is.corr = FALSE, method = "square", tl.col = "black")

dev.off()

# for xCell
fs %>% select(Tregs:B.cells) -> xcell_data
rownames(xcell_data) <- fs$ctype
xcell_data <- mutate_all(xcell_data, function(x) as.numeric(as.character(x)))
xcell_data <- xcell_data[,-11]
xcell_data <- data.frame(t(xcell_data))
colnames(xcell_data) <- fs$ctype
# for k means
# check for NA
sum(sapply(xcell_data, is.na))
# check for infinite
sum(sapply(xcell_data, is.infinite))
# check for NaN
sum(sapply(xcell_data, is.nan))

# calcuting optimal number of clusters by "Silhouette method"
# library(factoextra)
fviz_nbclust(xcell_data, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")+
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
FOOk2 <- kmeans(xcell_data, centers = 2, nstart = 25)
# we see the distrubution of the data in the 2 clusters
fviz_cluster(FOOk2, data = xcell_data, frame.type = "ellipse") + theme_minimal() + ggtitle("k = 2") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

# this is for KMEANS
# making a dataframe of FOOk$cluster
sample.order.data <- data.frame(FOOk2$cluster)
rownames(sample.order.data %>% filter(FOOk2.cluster==1)) -> cluster1Samples
rownames(sample.order.data %>% filter(FOOk2.cluster==2)) -> cluster2Samples
clusteredSamples2 <- c(cluster1Samples,rev(cluster2Samples))


rownames(xcell_data) -> xcell_data$cell_type
xcell_data[match(clusteredSamples2, xcell_data$cell_type),] -> xcell_data
xcell_data <- xcell_data[,-23]
xcell_data <- as.matrix(xcell_data)
xcell_data <- data.frame(xcell_data) %>% select(SKCM , PRAD, COADREAD, PAAD, BLCA, UCEC, LGG, STAD, TGCT, LUAD, BRCA, CESC, OV, ESCA, LUSC, GBM, PCPG, LIHC, HNSC, SARC, KIRP, KIRC)
xcell_data <- data.frame(t(xcell_data)) %>% select(clusteredSamples2)
colnames(xcell_data) <- c("Th1 cells", "Monocytes", "Macrophages M1", "Macrophages", "B cells", "CD4 memory T cells", "CD4 naive T cells", "CD4 T cells", "Central Memory CD4+ T cells", "Effector Memory CD4+ T cells", "CD8 T cells", "CD8 naive T cells", "Central Memory CD8+ T cells", "Macrophages M2", "Memory B cells", "Neutrophils", "NK cells", "Th2 cells", "Tregs")

#### doing the below with previously made xcell_data dataframe
xcell_data <- data.frame(t(xcell_data)) %>% select(SKCM , PRAD, COADREAD, PAAD, BLCA, UCEC, LGG, STAD, TGCT, LUAD, BRCA, CESC, OV, ESCA, LUSC, GBM, PCPG, LIHC, HNSC, SARC, KIRP, KIRC)

png("xcell.corrplot.bw.png", width=362.47916667, height=170.39166667, units = "mm", res = 118.11)

corrplot(as.matrix(xcell_data), is.corr = FALSE, method = "square",col=gray(16:0 / 16), bg="white", tl.col = "black")
# corrplot(as.matrix(xcell_data), is.corr = FALSE, method = "square", tl.col = "black")

dev.off()

# scatter plot of scores and cancer----------------
ggplot(wdata) + 
  geom_point(aes(x=ctype, y=score, color = ctype))
  # geom_jitter(aes(x=ctype, y=score, color=ctype), width = 0.1)

ggplot(data = wdata %>% filter(ctype == "BLCA")) + geom_point(aes(x=Sample.ID, y=score))

data %>% filter(ctype == "SKCM") -> cancer
cancer <- cancer[order(cancer$score), ]
cancer$Sample.ID <- factor(cancer$Sample.ID, levels = cancer$Sample.ID)
SKCM <- ggplot(data = cancer, aes(x=Sample.ID, y=score)) +
  geom_point(stat='identity', size=0.5) +
  geom_point(aes(x=as.integer(length(Sample.ID)/2), y=median(score)), shape=23, fill="white", color="black", size=1.5)+
  annotate("text", x=c(as.integer(length(cancer$Sample.ID)/4)), y= quantile(cancer$score,0.25), label="---", size=7, color = "black")+
  annotate("text", x=c(as.integer(3*length(cancer$Sample.ID)/4)), y= quantile(cancer$score,0.75), label="---", size=7, color = "black")+
  theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  xlab("SKCM")+ylim(-5,5)
cat("SKCM",":",median(cancer$score))

for (i in c("COADREAD", "BRCA" ,"CESC", "ESCA","GBM", "HNSC", "KIRC", "KIRP", "LGG", "LIHC", 'LUAD', "LUSC", "OV", "PAAD", "PCPG", "PRAD", "SARC", "BLCA", "STAD", "TGCT", "UCEC")){
  data %>% filter(ctype == i) -> cancer
  cancer <- cancer[order(cancer$score), ]
  cancer$Sample.ID <- factor(cancer$Sample.ID, levels = cancer$Sample.ID)
  assign(i, ggplot(data = cancer, aes(x=Sample.ID, y=score)) +
           geom_point(stat='identity', size=0.5) +
           geom_point(aes(x=as.integer(length(Sample.ID)/2), y=median(score)), shape=23, fill="white", color="black", size=1.5)+
           annotate("text", x=c(as.integer(length(cancer$Sample.ID)/4)), y= quantile(cancer$score,0.25), label="---", size=7, color = "black")+
           annotate("text", x=c(as.integer(3*length(cancer$Sample.ID)/4)), y= quantile(cancer$score,0.75), label="---", size=7, color = "black")+
           theme(axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())+
           xlab(i) + ylim(-5,5)
  )
  cat(i,":",median(cancer$score)," ")
}

library(patchwork)
#BLCA+BRCA+CESC+COADREAD+ESCA+GBM+HNSC+KIRC+KIRP+LGG+LIHC+LUAD+LUSC+OV+PAAD+PCPG+PRAD+SARC+SKCM+STAD+TGCT+UCEC+plot_layout(nrow = 1)

# divergent.dotplot <- SKCM + COADREAD+ UCEC+ SARC+ LIHC+
# LUSC+ PRAD+ OV+ PCPG+ GBM+ ESCA+
# BRCA+ CESC+ TGCT+  LUAD+  PAAD+  LGG+ BLCA+
# KIRP+  HNSC+  STAD+  KIRC + plot_layout(nrow = 1)

for (i in c("COADREAD","ESCA","TGCT","UCEC","BRCA","BLCA","STAD","KIRP","SKCM","KIRC","PRAD","PAAD","LUSC","OV","GBM","PCPG","LGG","LIHC","SARC","LUAD","CESC","HNSC")) {
  print(median((wdata %>% filter(ctype==i))$score))
}

#COADREAD + ESCA +TGCT+ UCEC+ BRCA+BLCA+STAD+KIRP+SKCM+KIRC+PRAD+PAAD+LUSC+OV+GBM+PCPG+LGG+LIHC+SARC+LUAD+CESC+HNSC+ plot_layout(nrow = 1)

divergent.dotplot <- SKCM + PRAD + COADREAD + PAAD + BLCA + UCEC + LGG + STAD + TGCT + LUAD + BRCA + CESC + OV + ESCA + LUSC + GBM + PCPG + LIHC + HNSC + SARC + KIRP + KIRC+ plot_layout(nrow = 1)

height = 362.47916667
width = 170.39166667
ggsave(plot = divergent.dotplot, filename = "divergent.dotplot.final.png", width=362.47916667, height=170.39166667, units = "mm", dpi = 300)
