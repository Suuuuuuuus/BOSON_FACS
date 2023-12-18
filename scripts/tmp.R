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


##
