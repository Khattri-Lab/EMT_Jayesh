#######################################################
# prepare a summary table of number of samples in each class and cancer type
# Jayesh Kumar Tiwari
# 17 July 2021
########################################################

# set home directory
home <- "/mnt/swift/Jayesh/emt/script.jayesh.edited/txt/"
library(tidyverse)

data <- read.table(paste0(home,"final.working.data_v3.txt",sep=""), header = T, sep = "\t")

crap <- c()
crapn <- c()
crapf <- c()
crapd <- c()
crapr <- c()
# measuring spread of cancer types
for (i in unique(data$ctype)) {
  crapn <- rbind(crapn, i)
  crap <- rbind(crap, IQR((data %>% filter(ctype==i))$score))
  crapr <- rbind(crapr, range((data %>% filter(ctype==i))$score))
  crapd <- rbind(crapd, range((data %>% filter(ctype==i))$score)[2]-range((data %>% filter(ctype==i))$score)[1])
  crapf <- cbind(crapn, crap, crapr, crapd)
}

crapf <- data.frame(crapf)
colnames(crapf) <- c("ctype","IQR","min","max","max-min")
write.table(crapf, file = "~/Dropbox/data.spread.txt", row.names = F, sep = '\t')

epithelial <- c()
mesenchymal <- c()
intermediate <- c()
total <- c()
count =0
for(i in unique(data$ctype)){
  count = count +1
  print(i)
  epithelial <- c(epithelial, length((data %>% filter(class=="epithelial" & ctype==i))$Sample.ID))
  mesenchymal <- c(mesenchymal, length((data %>% filter(class=="mesenchymal" & ctype==i))$Sample.ID))
  intermediate <- c(intermediate, length((data %>% filter(class=="intermediate" & ctype==i))$Sample.ID))
  total <- c(total, sum(length((data %>% filter(class=="epithelial" & ctype==i))$Sample.ID),
                        length((data %>% filter(class=="mesenchymal" & ctype==i))$Sample.ID),
                        length((data %>% filter(class=="intermediate" & ctype==i))$Sample.ID)))
  # check if sum is alright
  print(sum(length((data %>% filter(class=="epithelial" & ctype==i))$Sample.ID),
            length((data %>% filter(class=="mesenchymal" & ctype==i))$Sample.ID),
            length((data %>% filter(class=="intermediate" & ctype==i))$Sample.ID)))
  # check if total number of cancers are 22 at the end
  print(count)
}

summary_data <- data.frame("Abbreviation"=unique(data$ctype), 
                           "Epithelial"=epithelial, 
                           "Mesenchymal"=mesenchymal,
                           "Intermediate"=intermediate,
                           "Total"=total)
write.table(summary_data, file = "./txt/Data.Summary.txt", sep = '\t')
