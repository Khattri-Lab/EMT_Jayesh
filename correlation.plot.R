#########################################################
# scatter correlation plots for EMT score vs genes
# Jayesh Kumar Tiwari
# 29 July 2021
#####################################################

# set home directory
home <- "/mnt/swift/Jayesh/emt/script.jayesh.edited/"

# load libraries
library(tidyverse)

# load data
data <- read.table(paste0(home,"txt/final.working.data_v3.txt"), header = T, sep = '\t')

# plot correlation plot
# Add regression lines
corMacrophages.M1 <- ggplot(data, aes(x=score, y=Macrophages.M1)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

corTGFB1 <- ggplot(data, aes(x=score, y=TGFB1)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

corIL10 <- ggplot(data, aes(x=score, y=IL10)) +
  geom_point(aes(color=class)) + 
  geom_smooth(method=lm) +
  facet_wrap(~ctype)+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))


# save plot
ggsave("correlation.plot.Macrophages.M1.png", corMacrophages.M1, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.TGFB1.png", , width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
ggsave("correlation.plot.IL10.png", , width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
