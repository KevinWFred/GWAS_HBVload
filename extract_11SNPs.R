#!/usr/bin/env Rscript
library(data.table)
candidatesnps=c("rs6906021","rs148385846","rs187933986","rs76199331","rs805139",
                "rs9276417","rs11137882","rs80206952","rs79734975",
                "rs140224846","rs79351470")
candidateids=c("chr6:32658534:T:C","chr10:12762680:T:A","chr1:219552375:T:C",
               "chr2:54437895:T:C","chr5:9317648:C:T","chr6:32742935:G:C",
               "chr9:78863789:C:T","chr11:3030019:C:T","chr15:27792026:T:C",
               "chr15:28940922:T:C","chr18:70206830:G:A")
write.table(candidateids,file="../result/elevensnpsid.txt",row.names = F,col.names = F,quote=F)

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
cmd=paste0(plink," --bfile /data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed --extract ../result/elevensnpsid.txt --recode A-transpose --out ../result/elevensnps")
system(cmd)
snpdat=as.data.frame(fread("../result/elevensnps.traw"))
all(snpdat$SNP == candidateids)
idx=match(candidateids,snpdat$SNP)
snpdat=snpdat[idx,]
snpdat$SNP=candidatesnps
colnames(snpdat)=gsub("^0_","",colnames(snpdat))
snpdat=snpdat[,c(2,5:ncol(snpdat))]
colnames(snpdat)[which(colnames(snpdat)=="COUNTED")]="A1"
colnames(snpdat)[which(colnames(snpdat)=="ALT")]="A2"
write.csv(snpdat,file="../result/elevensnps.csv",row.names = F)
#compute Effect allele frequency (clearance vs. persistence)
seroclearance=read.csv("../data/id_seroclearance.csv")
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
idx=match(rownames(pheno),seroclearance$id_new)
pheno$sero=seroclearance$HBsAg_seroclearance[idx]
#AF1 is clearance
AF1=AF0=rep(NA,nrow(snpdat))
idx0=which(colnames(snpdat) %in% rownames(pheno)[pheno$sero==0])
idx1=which(colnames(snpdat) %in% rownames(pheno)[pheno$sero==1])
for (i in 1:nrow(snpdat))
{
  AF0[i]=sum(snpdat[i,idx0])/2/length(idx0)
  AF1[i]=sum(snpdat[i,idx1])/2/length(idx1)
}
AF1
# [1] 0.46822430 0.03644860 0.03084112 0.02990654 0.01962617
# [6] 0.03177570 0.01869159 0.02429907 0.03738318 0.04018692
# [11] 0.02149533
AF0
# [1] 0.429014901 0.054221854 0.012003311 0.010140728 0.008071192
# [6] 0.014486755 0.007243377 0.009933775 0.015107616 0.020695364
# [11] 0.007450331