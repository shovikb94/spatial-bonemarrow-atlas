# Script for plotting CFU assay data from SB68 - experiments performed by SB and MD 

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

cfu <- read_csv("CFU.csv")

cfu_long <- pivot_longer(cfu, cols = 2:4, names_to = "Sample", values_to = "Relative CFU")
cfu_long$Sample <- factor(cfu_long$Sample)
cfu_long$`Cell Type` <- factor(cfu_long$`Cell Type`, levels = c("Fibro", "Thy1", "Osteo", "Adipo"))
p1 <- cfu_long %>% ggplot(aes(x = `Cell Type`, y = `Relative CFU`)) + geom_boxplot() + geom_point() + theme_classic()
p1
ggsave(p1, device = "pdf", height = 7, width = 5, units = "in", filename = "~/Documents/NBM_Microenvironment/SB68_MSC_Sorting_FunctionalAssay/Data_Summary/cfu_boxplot.pdf")

fibro = as.numeric(cfu[1, 2:4])
thy1 = as.numeric(cfu[2, 2:4])
osteo = as.numeric(cfu[3, 2:4])
adipo = as.numeric(cfu[4, 2:4])

t.test(fibro, thy1)
t.test(fibro, osteo)
t.test(fibro, adipo)
