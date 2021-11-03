#########################################
# for plotting the pca of the 22 cancer types
# Jayesh Kumar Tiwari
# Modified: November 4, 2021

# Important: Skip to line 101 in the code, as the upper part is just experimental work
# and is kept for archival purposes, it is not useful in plotting the final figure.
##########################################


library(rgl)
require(SciViews)
library (plotrix)
library(corrplot)
require(ggplot2)
require(reshape)
require("gridExtra")

# read from del.data
pdata <- data.frame(t(pca_data.copy))
pdata.pca <- pcomp(pdata)
pdata_pca = cbind(cbind(pdata, pdata.pca$scores), genes = rownames(pdata))

screeplot(pdata.pca)
pca1 <- lala
screeplot(pca1)
summary(pdata.pca)

Loadings <- as.data.frame(pdata.pca$loadings[,1:3])
hm = cbind(cancer = rownames(Loadings), Loadings)
hmm = melt(hm, id = c("cancer"))
ggplot(hmm, aes(variable, cancer)) + geom_tile(aes(fill = value),
                                            colour = "white") + scale_fill_gradient(low = "red",
                                                                                    high = "green")

plot(pdata.pca, which = "correlations") 

# look at PC3
plot(pdata.pca, which = "correlations", choices = c(1,3))

# look at PC2 and PC3
plot(pdata.pca, which = "correlations", choices = c(2,3))

library(factoextra)
plot3d(scores(pdata.pca)[, 1:3])
lala <- prcomp(pdata, scale. = T)
fviz_pca_ind(lala,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(lala,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_biplot(lala, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

library(rgl)
text3d(pdata.pca$scores[,1:3], texts = colnames(pdata))
pca3d(lala, group = colnames(pdata))
options(rgl.printRglwidget = TRUE)

plot3d(lala, group= colnames(pdata))
rgl.viewpoint( theta = -58, phi = 20, fov = 60, zoom = 0.55, 
               scale = par3d("scale"), interactive = TRUE, 
               type = c("userviewpoint", "modelviewpoint"))
snapshotPCA3d(file="snapshotPCA3d.bmp")




######################### trying mclust ######################
# install.packages("mclust")
library(mclust)
clPairs(pdata, rownames(pdata))

BIC <- mclustBIC(pdata)
plot(BIC)

summary(BIC)

mod1 <- Mclust(pdata, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")

table(rownames(pdata), mod1$classification)

library(scatterplot3d)
scatterplot3d(lala$x[,1], lala$x[,2], lala$x[,3], main = "PCA 3D plot")










##################################################### THIS IS THE STARTING POINT OF THE REAL CODE ###########################################################



library(tidyverse)
# read final data
data <- read.table("final.working.data_v4.txt", header = T, sep = '\t')
# for testing again (edit made on 3 October 2021)
testing.data <- read.table("/mnt/swift/Jayesh/emt/script.jayesh.edited/txt/final.working.data_v4.txt", header = T, sep = '\t')
# output: 4028 obs of 1024 variables

# we only require sample id, ctype, class, checkpoints, cytokines, xcells, and gene signatures
data <- data %>% select(Sample.ID, ctype, class, 
                        C10orf54, CTLA4, CYBB.x, HAVCR2, ICOS.x, LAG3.x, PDCD1LG2, PDCD1, SIGLEC5, SIGLEC7, TNFRSF18, TNFRSF4, TNFRSF9, CD274, FASLG.x, KIR3DL1, SIGLEC15,
                        IL10, TGFB1, IFNA1, IFNB1, IFNG.x, IL12A, IL12B, IL1A, IL1B.x, IL2, IL3, IL4, IL5, IL6, IL8, STAT6, TNF,
                        Th1.cells, Monocytes, Macrophages, B.cells, CD4..memory.T.cells, CD4..naive.T.cells, CD4..T.cells, CD4..Tcm, CD4..Tem, CD8..T.cells, CD8..naive.T.cells, CD8..Tcm, Memory.B.cells, Neutrophils, NK.cells, Th2.cells, Tregs,
                        ifng6gene.sig, activated.stroma.sig, macrophages.sig, cd8.3gene.sig, cd8.9gene.sig, gajeski.sig, cytolytic.sig, hypoxia.sig)
# testing
testing.data <- testing.data %>% select(Sample.ID, ctype, class, 
                        C10orf54, CTLA4, CYBB.x, HAVCR2, ICOS.x, LAG3.x, PDCD1LG2, PDCD1, SIGLEC5, SIGLEC7, TNFRSF18, TNFRSF4, TNFRSF9, CD274, FASLG.x, KIR3DL1, SIGLEC15,
                        IL10, TGFB1, IFNA1, IFNB1, IFNG.x, IL12A, IL12B, IL1A, IL1B.x, IL2, IL3, IL4, IL5, IL6, IL8, STAT6, TNF,
                        Th1.cells, Monocytes, Macrophages, B.cells, CD4..memory.T.cells, CD4..naive.T.cells, CD4..T.cells, CD4..Tcm, CD4..Tem, CD8..T.cells, CD8..naive.T.cells, CD8..Tcm, Memory.B.cells, Neutrophils, NK.cells, Th2.cells, Tregs,
                        ifng6gene.sig, activated.stroma.sig, macrophages.sig, cd8.3gene.sig, cd8.9gene.sig, gajeski.sig, cytolytic.sig, hypoxia.sig)
# output: 4028 obs. of 62 variables

# creating dataframe of difference of medians i.e. 
# (median(expression of mesenchymal samples (for a particular cancer) - median(expression of epithelial samples (for a particular cancer))))
diffOfMedians <- data %>%
  group_by(ctype) %>%
  dplyr::summarize(across(C10orf54:hypoxia.sig, ~ median(.[class=="mesenchymal"] - median(.[class=="epithelial"]))))
# convert tiblle to dataframe
diffOfMedians <- data.frame(diffOfMedians)
# set cancer type as rownames and remove it from the dataframe
rownames(diffOfMedians) <- diffOfMedians$ctype
diffOfMedians <- diffOfMedians[,-1]

# testing
testing.diffOfMedians <- testing.data %>%
  group_by(ctype) %>%
  dplyr::summarize(across(C10orf54:hypoxia.sig, ~ median(.[class=="mesenchymal"] - median(.[class=="epithelial"]))))
# output: 22 obs. of 60 variables

# convert tiblle to dataframe
testing.diffOfMedians <- data.frame(testing.diffOfMedians)
# set cancer type as rownames and remove it from the dataframe
rownames(testing.diffOfMedians) <- testing.diffOfMedians$ctype
testing.diffOfMedians <- testing.diffOfMedians[,-1]
# output: 22 obs. of 59 variables

# doing pca
# pca_res <- prcomp(diffOfMedians, scale=TRUE)

# got this error:
# Error in prcomp.default(diffOfMedians, scale = TRUE) : 
#   cannot rescale a constant/zero column to unit variance

# and saw that IL3 (column 28) had zero in all rows, so removing that
testing.diffOfMedians <- testing.diffOfMedians[,-28]
# output: 22 obs. of 58 variables

# again doing pca
pca_res <- prcomp(diffOfMedians, scale. = T)

# testing
testing.pca_res <- prcomp(testing.diffOfMedians, scale. = T)
# output: List of 5

# look at the summary
summary(pca_res)
# Output:
# Importance of components:
#                         PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8
# Standard deviation     4.1018 2.7050 2.31163 2.19027 1.81739 1.76210 1.59434 1.49066
# Proportion of Variance 0.2901 0.1262 0.09213 0.08271 0.05695 0.05353 0.04383 0.03831
# Cumulative Proportion  0.2901 0.4162 0.50837 0.59108 0.64802 0.70156 0.74538 0.78370
#                         PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16
# Standard deviation     1.41585 1.35444 1.24314 1.16967 1.14128 1.05266 0.92378 0.84925
# Proportion of Variance 0.03456 0.03163 0.02664 0.02359 0.02246 0.01911 0.01471 0.01243
# Cumulative Proportion  0.81826 0.84989 0.87653 0.90012 0.92258 0.94168 0.95640 0.96883
#                         PC17    PC18    PC19   PC20   PC21      PC22
# Standard deviation     0.72961 0.68264 0.56241 0.5491 0.4377 5.308e-16
# Proportion of Variance 0.00918 0.00803 0.00545 0.0052 0.0033 0.000e+00
# Cumulative Proportion  0.97801 0.98604 0.99150 0.9967 1.0000 1.000e+00

names(pca_res)
# [1] "sdev"     "rotation" "center"   "scale"    "x"

# testing
summary(testing.pca_res)
# output:
# Importance of components:
#   PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8
# Standard deviation     4.1018 2.7050 2.31163 2.19027 1.81739 1.76210 1.59434 1.49066
# Proportion of Variance 0.2901 0.1262 0.09213 0.08271 0.05695 0.05353 0.04383 0.03831
# Cumulative Proportion  0.2901 0.4162 0.50837 0.59108 0.64802 0.70156 0.74538 0.78370
# PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16
# Standard deviation     1.41585 1.35444 1.24314 1.16967 1.14128 1.05266 0.92378 0.84925
# Proportion of Variance 0.03456 0.03163 0.02664 0.02359 0.02246 0.01911 0.01471 0.01243
# Cumulative Proportion  0.81826 0.84989 0.87653 0.90012 0.92258 0.94168 0.95640 0.96883
# PC17    PC18    PC19   PC20   PC21      PC22
# Standard deviation     0.72961 0.68264 0.56241 0.5491 0.4377 5.308e-16
# Proportion of Variance 0.00918 0.00803 0.00545 0.0052 0.0033 0.000e+00
# Cumulative Proportion  0.97801 0.98604 0.99150 0.9967 1.0000 1.000e+00

names(testing.pca_res)
# output: 
# [1] "sdev"     "rotation" "center"   "scale"    "x"  


# x object has PCs of our interest
pca_res$x[1:5,1:3]

# testing
testing.pca_res$x[1:5,1:3]
# output:
# PC1        PC2         PC3
# BLCA      -0.5917695  1.2155644  0.31322189
# BRCA       0.6043745 -0.6804462  2.32865993
# CESC       4.3010029  0.9675006  0.29343463
# COADREAD -10.8101295  2.3848535 -0.02425138
# ESCA      -4.2172438  3.5132498 -2.08324974

var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
var_explained[1:5]

# testing
testing.var_explained <- testing.pca_res$sdev^2/sum(testing.pca_res$sdev^2)
# output: num [1:22] 0.2901 0.1262 0.0921 0.0827 0.0569 ...
testing.var_explained[1:5]
# output: [1] 0.29008223 0.12615167 0.09213145 0.08271186 0.05694674

# plotting screeplot
screeplot(pca_res)

# testing
screeplot(testing.pca_res)
# output: saved as testing_output/screeplot_testing.pca_res.png

# first two PCs
pca_res$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(size=4) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

# testing
testing.pca_res$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(size=4) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")
# output: saved as testing_output/bw_pc1_pc2_unlabelled.png
# comment(by Jayesh): not a useful figure

fit <- hclust(dist(pca_res$x[,1:3]), method="complete") # 1:3 -> based on 3 components
groups <- cutree(fit, k=2)                            # k=2 -> 2 groups

# testing
testing.fit <- hclust(dist(testing.pca_res$x[,1:3]), method="complete") # 1:3 -> based on 3 components
# output: List of 7
testing.groups <- cutree(testing.fit, k=2)                            # k=2 -> 2 groups
# output: Named int [1:22] 1 2 2 1 1 1 2 2 2 2 ...

library(rgl)
plotPCA <- function(x, nGroup) {
  n <- ncol(x) 
  if(!(n %in% c(2,3))) { # check if 2d or 3d
    stop("x must have either 2 or 3 columns")
  }
  
  fit <- hclust(dist(x), method="complete") # cluster
  groups <- cutree(fit, k=nGroup)
  
  if(n == 3) { # 3d plot
    open3d(windowRect=c(100,100,700,700))
    plot3d(x, col=groups, type="s", size=1, axes=F)
    text3d(x, texts = rownames(x))
    axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
    grid3d("x")
    grid3d("y")
    grid3d("z")
  } else { # 2d plot
    maxes <- apply(abs(x), 2, max)
    rangeX <- c(-maxes[1], maxes[1])
    rangeY <- c(-maxes[2], maxes[2])
    plot(x, col=groups, pch=19, xlab=colnames(x)[1], ylab=colnames(x)[2], xlim=rangeX, ylim=rangeY)
    lines(c(0,0), rangeX*2)
    lines(rangeY*2, c(0,0))
    text2D(pca_res$x[,1]+0.2,pca_res$x[,2]+0.5,  labels = rownames(pca_res$x),
           add = TRUE, colkey = FALSE, cex = 0.5)
    
  }
}

plotPCA(pca_res$x[,1:2], 2)

# testing
# Note: Line 269 is a function, so need not be included in testing, can be just run from the same
plotPCA(testing.pca_res$x[,1:2], 2)
# output: saved as plotPCA_fxn_output.png
#
# Error in text2D(pca_res$x[, 1] + 0.2, pca_res$x[, 2] + 0.5, labels = rownames(pca_res$x),  : 
# could not find function "text2D"
#
# comment(by Jayesh): Can act as validator for the final pca plot


library(plot3D)
scatter3D(pca_res$x[,1],pca_res$x[,2],pca_res$x[,3], phi = 0, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed")
# Add text

# testing
testing.scatter3D(testing.pca_res$x[,1],testing.pca_res$x[,2],testing.pca_res$x[,3], phi = 0, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed")
# output:
# Error in testing.scatter3D(testing.pca_res$x[, 1], testing.pca_res$x[,  : 
# could not find function "testing.scatter3D"

# comment(by Jayesh): threw error and not useful anyways

text3D(pca_res$x[,1],pca_res$x[,2],pca_res$x[,3],  labels = rownames(pca_res$x),
       add = TRUE, colkey = FALSE, cex = 0.5)
playwidget()
rgl.postscript( "plot", fmt = "eps", drawText = TRUE )
snapshot3d( filename = "3dplot.png", 
            fmt = "png")
typeof(pca_res$x)
rglwidget()

# testing
text3D(testing.pca_res$x[,1],testing.pca_res$x[,2],testing.pca_res$x[,3],  labels = rownames(testing.pca_res$x),
       add = TRUE, colkey = FALSE, cex = 0.5)
playwidget()
rgl.postscript( "plot", fmt = "eps", drawText = TRUE )
snapshot3d( filename = "3dplot.png", 
            fmt = "png")
typeof(pca_res$x)
rglwidget()

# output:
# Error in cbind(as.vector(x), as.vector(y), as.vector(z), 1) %*% plist$mat : 
#   requires numeric/complex matrix/vector arguments

# comment(by Jayesh): was not used finally, so doesn't matter



pca_res.dataframe <- data.frame(pca_res$x)
pca_res.dataframe$ctype <- rownames(pca_res.dataframe)

# testing
testing.pca_res.dataframe <- data.frame(testing.pca_res$x)
# output: 22 obs. of 22 variables
testing.pca_res.dataframe$ctype <- rownames(testing.pca_res.dataframe)
# output: 22 obs. of 23 variables


library(ggrepel)
library(ggsci)
ggplot(data = pca_res.dataframe, aes(x= PC1, y = PC2, color = PC1)) + geom_point() + 
  # scale_color_gradient2(midpoint=-5, low="red", mid = "orange",
  #                      high="blue", space ="Lab" ) + 
  scale_color_gradient2(midpoint=-5, low="grey", mid = "dark grey",
                       high="black", space ="Lab" ) +

  # scale_color_gradientn(colours = rainbow(5)) +
  # scale_color_gradient2(low="red",mid = "orange",  high="blue") +
  # geom_text(aes(label=rownames(pca_res.dataframe)),hjust=0, vjust=0)
  geom_label_repel(aes(label = rownames(pca_res.dataframe)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()


# testing
ggplot(data = testing.pca_res.dataframe, aes(x= PC1, y = PC2, color = PC1)) + geom_point() + 
  # scale_color_gradient2(midpoint=-5, low="red", mid = "orange",
  #                      high="blue", space ="Lab" ) + 
  scale_color_gradient2(midpoint=-5, low="grey", mid = "dark grey",
                        high="black", space ="Lab" ) +
  
  # scale_color_gradientn(colours = rainbow(5)) +
  # scale_color_gradient2(low="red",mid = "orange",  high="blue") +
  # geom_text(aes(label=rownames(pca_res.dataframe)),hjust=0, vjust=0)
  geom_label_repel(aes(label = rownames(pca_res.dataframe)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()
# output:
# Warning message:
#   ggrepel: 2 unlabeled data points (too many overlaps). Consider increasing max.overlaps

# saved as ggplot.pca_res.png


library(ggrepel)
pca.plot <- ggplot(data = pca_res.dataframe, aes(x= PC1, y = PC2, color = PC1>-0.6)) + geom_point() + 
  scale_colour_manual(name = 'PC1 > -0.6', values = setNames(c('dark green','red'),c(T, F)))+
  geom_text_repel(aes(label = rownames(pca_res.dataframe)), hjust = -0.3,  vjust = .8, segment.color = NA)+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0), axis.title = element_blank())+
  theme_classic()

# testing
testing.pca.plot <- ggplot(data = testing.pca_res.dataframe, aes(x= PC1, y = PC2, color = PC1>-0.6)) + geom_point() + 
  scale_colour_manual(name = 'PC1 > -0.6', values = setNames(c('dark green','red'),c(T, F)))+
  geom_text_repel(aes(label = rownames(testing.pca_res.dataframe)), hjust = -0.3,  vjust = .8, segment.color = NA)+
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain", angle = 0),axis.text.y = element_text(color = "black", size = 12, face = "plain", angle = 0), axis.title = element_blank())+
  theme_classic()
testing.pca.plot

# output
# saved as final_pca_plot.png

ggsave(plot = pca.plot, filename = "all.cancers.pca.png", width=342.37083333, height=162.71875, units = "mm", dpi = 300)

