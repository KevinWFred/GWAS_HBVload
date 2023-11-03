#!/usr/bin/env Rscript
library(data.table)
hbv_gwas=as.data.frame(fread("../result/sumstat.cvs",sep=","))

resultfolder="/data/BB_Bioinformatics/Kevin/HBV_GWAS/result/lmm/"
allfiles = list.files(resultfolder,pattern = "lmm*")
alljobs = as.numeric(gsub("lmm","",allfiles))
if (length(unique(alljobs))!=max(alljobs)) warning("Some results are missing")
lmmres = NULL
i=1
for (myfile in paste0(resultfolder,allfiles))
{
  if (i %% 100==0) cat(i,'..')
  if(file.exists(myfile))
  {
    tmp = fread(myfile,header = T)
    lmmres = rbind(lmmres,tmp)
  }
  i = i+1
}
#align alleles
lmmtraw=as.data.frame(fread("../result/ordinalr_selsnps.traw")) 
idx=match(lmmres$SNP,lmmtraw$SNP)
lmmtraw=lmmtraw[idx,]
all(lmmtraw$SNP==lmmres$SNP)
lmmres$a1=lmmtraw$COUNTED
lmmres$a2=lmmtraw$ALT
bim=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed.bim"))
idx=match(lmmres$SNP,bim$V2)
idx1=which(lmmtraw$COUNTED==bim$V6[idx]) #need to flip (ordinal results are V5)
tmp=lmmres$a1[idx1]
lmmres$a1[idx1]=lmmres$a2[idx1]
lmmres$a2[idx1]=tmp
lmmres$beta[idx1]=-lmmres$beta[idx1]
write.csv(lmmres[,c(1,2,6,7,3:5)],file="../result/HBVload_LMM_sumstat.csv",row.names=F,quote=F)
lmmres1=lmmres[!is.na(lmmres$P),]
lmmres1$P[which(lmmres1$P==0)]=10^(-200) #8


idx=match(lmmres1$SNP,hbv_gwas$SNP)
allres=cbind.data.frame(hbv_gwas[idx,1:6],log_beta=hbv_gwas$Beta[idx],log_P=hbv_gwas$P[idx],
                        lmmres1[,c(3,5)])

png(filename="../result/scater_logis_lmm_hbv_pvalue.png",res=100)
par(mar=c(4.5,5,2,1))
ylim=range(c(-log10(hbv_gwas$P[idx]),-log10(lmmres1$P)))
plot(-log10(hbv_gwas$P[idx]),-log10(lmmres1$P),xlab="-log10(Logistic p-value)",ylab="-log10(LMM pvalue)",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_lmm_hbv_pvalue1.png",res=100)
par(mar=c(4.5,5,2,1))
ylim=c(0,17)
plot(-log10(hbv_gwas$P[idx]),-log10(lmmres1$P),xlab="-log10(Logistic p-value)",ylab="-log10(LMM pvalue)",cex.lab=1.5,cex.axis=1.5,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_lmm_beta.png",res=100)
par(mar=c(4.5,5,2,1))
ylim=range(c(hbv_gwas$Beta[idx],lmmres1$beta))
plot(hbv_gwas$Beta[idx],lmmres1$beta,xlab="Logistic beta",ylab="LMM beta",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
