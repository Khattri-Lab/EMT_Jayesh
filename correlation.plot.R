#########################################################
# scatter correlation plots for EMT score vs genes
# Jayesh Kumar Tiwari
# 29 July 2021
#####################################################

# set home directory
home <- "/mnt/swift/Jayesh/emt/script.jayesh.edited/"

# load libraries
library(tidyverse)
library(ggpubr)

# load data
data <- read.table(paste0(home,"txt/final.working.data_v3.txt"), header = T, sep = '\t')
# for gene signatures
data <- read.table(paste0(home,"txt/final.working.data_v4.txt"), header = T, sep = '\t')

# setting the right order of cancers
data$ctype <- factor(data$ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC"))

# plot correlation plot
# Add regression lines
# stat_cor reference (https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html)
corMacrophages.M1 <- ggplot(data, aes(x=score, y=Macrophages.M1)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  # below commented code adds an x and y axis to the plot
  # geom_hline(yintercept=0)+
  # geom_vline(xintercept=0)+
  # scale_x_continuous(expand=c(0,0))+
  # scale_y_continuous(expand=c(0,0))+
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -6,
    size = 3
  )


corTGFB1 <- ggplot(data, aes(x=score, y=TGFB1)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )


corIL10 <- ggplot(data, aes(x=score, y=IL10)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(panel.border = element_blank(),
        strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )


corMacrophages <- ggplot(data, aes(x=score, y=Macrophages)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )

corSIGLEC7 <- ggplot(data, aes(x=score, y=SIGLEC7)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )


corCD89gene <- ggplot(data, aes(x=score, y=cd8.9gene.sig)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )


corActivatedStroma <- ggplot(data, aes(x=score, y=activated.stroma.sig)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = -3,
    size = 3
  )

# save plot
ggsave("correlation.plot.Macrophages.M1.png", corMacrophages.M1, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.TGFB1.png", corTGFB1, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.IL10.png", corIL10, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.Macrophages.png", corMacrophages, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.SIGLEC7.png", corSIGLEC7, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.CD89gene.png", corCD89gene, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.ActivatedStroma.png", corActivatedStroma, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
