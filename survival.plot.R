##################################
# survival analysis of all cancer types
# Jayesh Kumar tiwari
# 26 June 2021
##################################

# set home directory for text files
home <- "/mnt/swift/Jayesh/emt/script.jayesh.edited/txt/"

# loading libraries
library(survival)
library(survminer)
library(tidyverse)

################## UPDATE #######################
# load the below file only and skip to 'creating an OS fit' section
data.survival <- read.table(paste0(home,"final.working.data_v3.with.survival.txt",sep=''), header = T, sep = '\t')

################# Skip the below lines and go to 'creating an OS fit' section ############
# loading survival data
load("/mnt/swift/Jayesh/emt/script.jayesh.edited/txt/clinical_survival_pancancer_atlas (1).RData")
# above dataframe is stored as dat

dat$bcr_patient_barcode <- gsub('-','.', dat$bcr_patient_barcode)
dat$bcr_patient_barcode <- paste0(dat$bcr_patient_barcode,".01")

# loading final.working.data_v3.txt
data <- read.table(paste0(home,"final.working.data_v3.txt",sep=''), header = T, sep = '\t')

data.survival <- merge(data, dat, by.x = "Sample.ID", by.y = "bcr_patient_barcode")
# removing intermediate samples
data.survival <- data.survival %>% filter(class != "intermediate") %>% select(Sample.ID, class, ctype, OS, OS.time, PFI, PFI.time, everything())

# # write this file as final.working.data_v3.with.survival.txt
# write.table(data.survival, file = paste0(home,"final.working.data_v3.with.survival.txt",sep=''), sep = '\t')

# creating an OS fit----------------------------------------------------
# fix order of cancers
data.survival$ctype <- factor(data.survival$ctype, levels = c("SKCM" , "PRAD", "COADREAD", "PAAD", "BLCA", "UCEC", "LGG", "STAD", "TGCT", "LUAD", "BRCA", "CESC", "OV", "ESCA", "LUSC", "GBM", "PCPG", "LIHC", "HNSC", "SARC", "KIRP", "KIRC"))

fitOS <- survfit(Surv(OS.time, OS) ~ class, data = data.survival)
print(fitOS)
# creating a fit for PFI
fitPFI <- survfit(Surv(PFI.time, PFI) ~ class, data = data.survival)
print(fitPFI)

# plotting survival curves for all cancer types
# plot.os.all <- ggsurvplot(fitOS,
#            pval = TRUE, conf.int = TRUE,
#            # risk.table = TRUE, # Add risk table
#            # risk.table.col = "strata", # Change risk table color by groups
#            #linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"),
#            title="Overall Survival across all cancer types",
#            facet.by = "ctype")
# 
# plot.pfi.all <-  ggsurvplot(fitPFI,
#            pval = TRUE, conf.int = TRUE,
#            # risk.table = TRUE, # Add risk table
#            # risk.table.col = "strata", # Change risk table color by groups
#            #linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"),
#            title="Progression Free Survival across all cancer types",
#            facet.by = "ctype")


plot.os.all <- ggsurvplot(fit = fitOS, size = 2,  pval = T, title='Overall Survival across all cancer types', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green'), facet.by = "ctype") #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("./plt/OS.all_cancer_types.png", plot.os.all, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)


plot.pfi.all <- ggsurvplot(fit = fitPFI, size = 2,  pval = T, title='Progression Free Survival across all cancer types', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green'), facet.by = "ctype") #+
# labs(title = 'Overall Survival across all cancer types')
ggsave("./plt/PFI.all_cancer_types.png", plot.pfi.all, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)


# KIRC and LGG have significant p-values so make individual survival for them
fitOSkirc <- survfit(Surv(OS.time, OS) ~ class, data = data.survival%>%filter(ctype=="KIRC"))
fitOSlgg <- survfit(Surv(OS.time, OS) ~ class, data = data.survival%>%filter(ctype=="LGG"))

fitPFIkirc <- survfit(Surv(PFI.time, PFI) ~ class, data = data.survival%>%filter(ctype=="KIRC"))
fitPFIlgg <- survfit(Surv(PFI.time, PFI) ~ class, data = data.survival%>%filter(ctype=="LGG"))
fitPFIkirp <- survfit(Surv(PFI.time, PFI) ~ class, data = data.survival%>%filter(ctype=="KIRP"))


# change the fitOS thing for different plots
ggsurvplot(fitOSlgg,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           title="Overall Survival across all cancer types",
           facet.by = "ctype")

ggsurvplot(fitPFIkirp,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("black", "dark grey"),
           title="Progression Free Survival across all cancer types",
           facet.by = "ctype")


# kirc and lgg
plot.os.kirc <- ggsurvplot(fit = fitOSkirc, size = 2,  pval = T, title='Overall Survival in Kidney renal clear cell carcinoma', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green')) #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("OS.kirc.png", plot.os.kirc$plot, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)

plot.os.lgg <- ggsurvplot(fit = fitOSlgg, size = 2,  pval = T, title='Overall Survival in Brain Lower Grade Glioma', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green')) #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("OS.lgg.png", plot.os.lgg$plot, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)

# pfi plots
plot.pfi.kirc <- ggsurvplot(fit = fitPFIkirc, size = 2,  pval = T, title='Progression Free Survival in Kidney renal clear cell carcinoma', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green')) #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("PFI.kirc.png", plot.pfi.kirc$plot, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)

plot.pfi.lgg <- ggsurvplot(fit = fitPFIlgg, size = 2,  pval = T, title='Progression Free Survival in Brain Lower Grade Giloma', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green')) #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("PFI.lgg.png", plot.pfi.lgg$plot, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)

plot.pfi.kirp <- ggsurvplot(fit = fitPFIkirp, size = 2,  pval = T, title='Progression Free Survival in Kidney renal papillary cell carcinoma', ggtheme = theme_classic2(base_size=10), palette = c('red', 'blue', 'green')) #+
#labs(title = 'Overall Survival across all cancer types')
ggsave("PFI.kirp.png", plot.pfi.kirp$plot, width=345.01666667, height=189.44166667, units = "mm", dpi = 300)
