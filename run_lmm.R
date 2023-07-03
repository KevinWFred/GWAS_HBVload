#!/usr/bin/env Rscript
setwd("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")
selsnps=read.csv("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/suggestiveSNPs.csv")
idxs=seq(1,nrow(selsnps),10)
if (max(idxs)<nrow(selsnps)) idxs=c(idxs,nrow(selsnps)+1)
tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/lmm.R",length(idxs)-1),nstart=idxs[-length(idxs)],nend=idxs[2:length(idxs)]-1,
               outprefix=paste0("../result/lmm/lmm",1:(length(idxs)-1)))
write.table(tmp,file="log/lmm.swarm",row.names = F,col.names = F,sep=" ",quote=F)
"swarm -f lmm.swarm -g 16 --module R/4.2.0 --time=2-10:00:00 --gres=lscratch:16"

tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/lmm.R",length(idxs)-1),nstart=idxs[-length(idxs)],nend=idxs[2:length(idxs)]-1,
               outprefix=paste0("../result/lmm_hbeag_neg/lmm",1:(length(idxs)-1)),hbeag="negative")
write.table(tmp,file="log/lmm_hbeag_neg.swarm",row.names = F,col.names = F,sep=" ",quote=F)
"swarm -f lmm_hbeag_neg.swarm -g 16 --module R/4.2.0 --time=2-10:00:00 --gres=lscratch:16"