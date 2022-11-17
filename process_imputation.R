#!/usr/bin/env Rscript
#generate phenotype data used for GWAS
setwd("/data/DCEGLeiSongData/Kevin/code")
pheno = read.table("../data/ana_final_VL_new.csv",sep=",",header=T)

#used to select 3240 samples
tmp = data.frame(fam=0,iid=pheno$subject_id)
write.table(tmp,file = "../data/samples.txt",row.names = F,col.names = F,quote = F)

rownames(pheno) = pheno$subject_id
colnames(pheno) = gsub("^X","",colnames(pheno))
pheno = pheno[,c(paste0("PC",1:5),"SEX","AGE.grp","HBVDNA.grp")]
pheno$SEX = factor(pheno$SEX)
pheno$AGE.grp = factor(pheno$AGE.grp)
pheno$HBVDNA.grp = factor(pheno$HBVDNA.grp,ordered = T)
#save as Rdata
save(pheno,file="../result/pheno.RData")
