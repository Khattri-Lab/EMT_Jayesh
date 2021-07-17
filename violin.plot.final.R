#########################
# plot violin plots for the genes and cells
# Jayesh Kumar Tiwari
# 17 July 2021
#########################

# set home directory for txt files
ome <- "/mnt/swift/Jayesh/script.jayesh.edited/txt/"

library(dgof)
library(tidyverse)
library(ggsci)

# p-values to stars function
stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}


# read data
data <- read.table(paste0(home,'final.working.data_v3.txt',sep=''), header = T, sep = '\t')
data.copy <- data

initial_violin_plotter <- function(gene){
  # @ gene: the gene or cell which you want to plot
  # @ returns: a list containing: 
  #             1) a raw plot to decide the y limits and further add p-value
  #             2) a processed that dataframe that will be use by the final_plotter function
  #             3) gene name. It will be used in the function final_plotter
  #            This list is being outputted so that the user does not have to pass
  #             so many arguments in the final_plotter function and the objects are also preserved
  # for testing let's say gene is IL10
  # gene <- "IL10"
  data %>% select(Sample.ID, score, class, ctype, all_of(gene)) -> data.copy
  data.copy$class <- as.factor(data.copy$class)
  data.copy <- data.copy %>% filter(class != "intermediate")
  # view command insid function because data.copy is a local variable, and out side
  # the function it will not be saved
  View(data.copy)
  
  # plot without p-values
  p <- ggplot(data.copy, aes(x=factor(ctype, 
                                  levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC")), 
                         y= get(gene), group=interaction(ctype, class)))+ 
    geom_violin(scale = 'width', 
                position = position_dodge(width = .9), 
                aes(fill = class), size = 0.25) + 
    geom_boxplot(width=0.1, 
                 position = position_dodge(width = .9), 
                 outlier.shape = NA , fill='black') + 
    stat_summary(fun.y=median, 
                 geom="point", 
                 size=1, 
                 color="white", 
                 position = position_dodge(width = .9)) + 
    scale_fill_tron() + 
    # YOU NEED TO ADD y limits after looking at p, so just returning p so that
    # one can view p then add y limits by themself
    # ylim(-1,4) + # for LAG 3 only
    # ylim(-1.6,7) + # for SIGLEC 7 only
    # ylim(-1.3, 5) + # for TNFRSF 4 only
    theme_bw() + xlab("Cancer")+ylab(gene)+
    theme(axis.text.x = element_text(color = "black", 
                                     face = "plain", 
                                     angle = 90, 
                                     vjust = 0.5, 
                                     hjust=1),
          axis.text.y = element_text(color = "black", 
                                     face = "plain", 
                                     angle = 0, 
                                     vjust = 0.5, 
                                     hjust=1), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  # p
  # return a list containing the ggplot object and the processed dataframe
  rlist <- list("plot"=p, "sorted_data"=data.copy, "gene"=gene)
  return(rlist)
}

final_plotter <- function(p, data.copy, gene ,ylimit){
  # @param p: ggplot object obtained from initial_violin_plotter
  # @param data.copy: a processed dataframe having Sample.ID, score, class, ctype and gene as columns
  # @param gene: the name of gene or cell
  # @param ylimit: vector having lower and upper y limit
  
  p <- p+ylim(ylimit)
  
  spreaded.data.copy <- data.copy %>% tidyr::spread(key = "ctype", value = gene)
  spreaded.data.copy %>% dplyr::summarise(across((BLCA:UCEC), ~(ks.test(.[class == "epithelial"], .[class == "mesenchymal"])$p.value))) -> pvalues
  pvalues %>% t() %>% as.data.frame() -> pvalues.transposed
  pvalues.transposed$V1 <- round(pvalues.transposed$V1,digits = 5)
  pvalues.transposed$V1 -> pvalues.only
  pvalues.transposed -> pva  
  stars.pval(pvalues.only) -> pva.dot
  
  f <- p + annotate('text', x = c(1:22),
                              y=c(ylimit[2]),
                    label=c(pva.dot[19], pva.dot[17], pva.dot[4], pva.dot[15], pva.dot[1],
                            pva.dot[22], pva.dot[10], pva.dot[20], pva.dot[21], pva.dot[12],
                            pva.dot[2], pva.dot[3], pva.dot[14], pva.dot[5], pva.dot[13],
                            pva.dot[6], pva.dot[16], pva.dot[11], pva.dot[7], pva.dot[18],
                            pva.dot[9], pva.dot[8]))
  pva.dota <- c('___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___','___')
  #h <-  f +annotate('text', x=c(1), y=c(.3), label=pva.dota[1])+ annotate('text', x=c(2), y=c(.35), label=pva.dota[2])+annotate('text', x=c(3), y=c(.28), label=pva.dota[3])+annotate('text', x=c(4), y=c(.3), label=pva.dota[4])+annotate('text', x=c(5), y=c(.3), label=pva.dota[5])+annotate('text', x=c(6), y=c(.45), label=pva.dota[6])+annotate('text', x=c(7), y=c(.35), label=pva.dota[7])+annotate('text', x=c(8), y=c(.52), label=pva.dota[8])+annotate('text', x=c(9), y=c(.5), label=pva.dota[9])+annotate('text', x=c(10), y=c(.3), label=pva.dota[10])+annotate('text', x=c(11), y=c(.45), label=pva.dota[11])+annotate('text', x=c(12), y=c(.5), label=pva.dota[12])+annotate('text', x=c(13), y=c(.5), label=pva.dota[13])+annotate('text', x=c(14), y=c(.3), label=pva.dota[14])+annotate('text', x=c(15), y=c(.35), label=pva.dota[15])+annotate('text', x=c(16), y=c(.3), label=pva.dota[16])+annotate('text', x=c(17), y=c(.17), label=pva.dota[17])+annotate('text', x=c(18), y=c(.4), label=pva.dota[18])+annotate('text', x=c(19), y=c(.3), label=pva.dota[19])+annotate('text', x=c(20), y=c(.3), label=pva.dota[20])+annotate('text', x=c(21), y=c(.48), label=pva.dota[21])+annotate('text', x=c(22), y=c(.27), label=pva.dota[22])+annotate('text', x=c(23), y=c(.27), label=pva.dota[23])
  h <- f + annotate('text', x = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),
                    y=c(ylimit[2]),
                    label=c(pva.dota[1], pva.dota[2], pva.dota[3], pva.dota[4], pva.dota[5],
                            pva.dota[6], pva.dota[7], pva.dota[8], pva.dota[9], pva.dota[10],
                            pva.dota[11], pva.dota[12], pva.dota[13], pva.dota[14], pva.dota[15],
                            pva.dota[16], pva.dota[17], pva.dota[18], pva.dota[19], pva.dota[20],
                            pva.dota[21], pva.dota[22]))

}

########################### PLOT THE GENES and CELLS ################
# for IL10--------------------
# plot initial
raw_list <- initial_violin_plotter("IL10")
# view the raw plot and look at data.copy's gene column to determine 
# upper and lower y limits for the final plot
raw_list$plot
# set y limit
ylimit = c(-1.21,5)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
# view final plot
final_plot
# save the plot
ggsave(paste("./plt/IL10.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 34 rows containing non-finite values (stat_ydensity). 
# 2: Removed 34 rows containing non-finite values (stat_boxplot). 
# 3: Removed 34 rows containing non-finite values (stat_summary).

# for TGFB1--------------------------
raw_list <- initial_violin_plotter("TGFB1")
raw_list$plot
ylimit = c(-2.1,6)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/TGFB1.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 8 rows containing non-finite values (stat_ydensity). 
# 2: Removed 8 rows containing non-finite values (stat_boxplot). 
# 3: Removed 8 rows containing non-finite values (stat_summary).


# for CD4 naive T cells-----------------
raw_list <- initial_violin_plotter("CD4..naive.T.cells")
raw_list$plot
ylimit = c(-0.01,0.5)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/CD4.naive.T.cell.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 12 rows containing non-finite values (stat_ydensity). 
# 2: Removed 12 rows containing non-finite values (stat_boxplot). 
# 3: Removed 12 rows containing non-finite values (stat_summary).

# for CTLA4----------------------
raw_list <- initial_violin_plotter("CTLA4")
raw_list$plot
ylimit = c(-1.09,6)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/CTLA4.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 23 rows containing non-finite values (stat_ydensity). 
# 2: Removed 23 rows containing non-finite values (stat_boxplot). 
# 3: Removed 23 rows containing non-finite values (stat_summary).

# for CTLA4----------------------
raw_list <- initial_violin_plotter("CTLA4")
raw_list$plot
ylimit = c(-1.09,6)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/CTLA4.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 23 rows containing non-finite values (stat_ydensity). 
# 2: Removed 23 rows containing non-finite values (stat_boxplot). 
# 3: Removed 23 rows containing non-finite values (stat_summary).

# for PD1----------------------
raw_list <- initial_violin_plotter("PDCD1")
raw_list$plot
ylimit = c(-1,5.6)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot <- final_plot+ylab("PD-1")
ggsave(paste("./plt/PD1.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 29 rows containing non-finite values (stat_ydensity). 
# 2: Removed 29 rows containing non-finite values (stat_boxplot). 
# 3: Removed 29 rows containing non-finite values (stat_summary).

# for NOX2----------------------
raw_list <- initial_violin_plotter("CYBB")
raw_list$plot
ylimit = c(-1.23,6)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot <- final_plot+ylab("NOX2")
ggsave(paste("./plt/NOX2.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 17 rows containing non-finite values (stat_ydensity). 
# 2: Removed 17 rows containing non-finite values (stat_boxplot). 
# 3: Removed 17 rows containing non-finite values (stat_summary).

# for SIGLEC5----------------------
raw_list <- initial_violin_plotter("SIGLEC5")
raw_list$plot
ylimit = c(-1.2,5)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/SIGLEC5.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 23 rows containing non-finite values (stat_ydensity). 
# 2: Removed 23 rows containing non-finite values (stat_boxplot). 
# 3: Removed 23 rows containing non-finite values (stat_summary).

# for MACROPHAGES----------------------
raw_list <- initial_violin_plotter("Macrophages")
raw_list$plot
ylimit = c(-0.01,0.4)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/Macrophages.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 5 rows containing non-finite values (stat_ydensity). 
# 2: Removed 5 rows containing non-finite values (stat_boxplot). 
# 3: Removed 5 rows containing non-finite values (stat_summary).

# for Macrophages.M1----------------------
raw_list <- initial_violin_plotter("Macrophages.M1")
raw_list$plot
ylimit = c(-0.01,0.28)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/Macrophages.M1.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)

# for Macrophages.M2----------------------
raw_list <- initial_violin_plotter("Macrophages.M2")
raw_list$plot
ylimit = c(-0.01,0.14)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot
ggsave(paste("./plt/Macrophages.M2.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 31 rows containing non-finite values (stat_ydensity). 
# 2: Removed 31 rows containing non-finite values (stat_boxplot). 
# 3: Removed 31 rows containing non-finite values (stat_summary).

# for CD8..Tem----------------------
raw_list <- initial_violin_plotter("CD8..Tem")
raw_list$plot
summary(raw_list$sorted_data)
median(raw_list$sorted_data$CD8..Tem)
mean(raw_list$sorted_data$CD8..Tem)
sd(raw_list$sorted_data$CD8..Tem)
ylimit = c(-0.01,0.05)
final_plot <- final_plotter(raw_list$plot, raw_list$sorted_data, raw_list$gene, ylimit)
final_plot+ylab("Effector memory CD8+ T cells")
ggsave(paste("./plt/CD8.Tem.violin.plot.with.pvalues.final",".png",sep = ""), final_plot, width = (3508/300), height = (701/300), dpi = 300)
# for above: 
# Warning messages:
# 1: Removed 83 rows containing non-finite values (stat_ydensity). 
# 2: Removed 83 rows containing non-finite values (stat_boxplot). 
# 3: Removed 83 rows containing non-finite values (stat_summary).
