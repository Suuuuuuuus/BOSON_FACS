args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]
plot <- args[2]
indir <- args[3]
plotdir <- args[4]

library(qqman)
library(data.table)
library(tidyverse)
library(CMplot)
source("../scripts/kactk.R")

#setwd("Desktop/Boson FACS/boson_vcf/")
#marker = "CD3+CD57+"

res_file = paste0(indir , marker, ".assoc.linear")

df = fread(res_file)
res = data.frame(df$SNP, df$CHR, df$BP, df$P)
names(res)<-c("SNP", "CHR", "BP", "P")

if (plot == "TRUE"){

  setwd(plotdir)

  CMplot(res,plot.type = "q",threshold = 0.05, file = "jpg", file.name = marker, dpi = 300)
  CMplot(res,plot.type = "m", LOG10=TRUE, amplify = T, col = c("grey30", "grey60"), 
         file = "jpg", file.name = marker, 
         signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"), dpi = 300)
}

est = estlambda(df$P)$estimate
res = paste0(marker, "\t", est)
write(res, file = "../../boson_vcf/genomic_inflation_factor_all.txt", append = TRUE)

#setwd(dirname(normalizePath(getwd())))
