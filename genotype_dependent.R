#!/usr/bin/env Rscript
library(data.table)

pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
gendat=read.csv("../data/genotype.csv")
table(rownames(pheno) %in% gendat$ID)
# TRUE 
# 3240
idx=match(rownames(pheno),gendat$ID)
pheno$HBVgenotype=gendat$HBV_genotype[idx]
table(pheno$HBVgenotype)
# 0    B   BC    C    X 
# 140 1510   83  835  672 