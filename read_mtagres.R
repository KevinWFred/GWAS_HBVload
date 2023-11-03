#!/usr/bin/env Rscript
setwd("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")
library(data.table)
sero_mtag=as.data.frame(fread("../result/mtag_hbvload_sero_trait_1.txt"))
hbv_mtag=as.data.frame(fread("../result/mtag_hbvload_sero_trait_2.txt"))
sero_gwas=as.data.frame(fread("../result/HBsAg_seroclearance_sumstat.csv"))
dim(sero_gwas) #5636103
idx=which((sero_gwas$REF %in% c("A","T","G","C")) & (sero_gwas$ALT %in% c("A","T","G","C")) )
sero_gwas=sero_gwas[idx,]
dim(sero_gwas) #5272982
hbv_gwas=as.data.frame(fread("../result/sumstat.cvs",sep=","))
idx=which((hbv_gwas$EFF %in% c("A","T","G","C")) & (hbv_gwas$REF %in% c("A","T","G","C")) )
hbv_gwas=hbv_gwas[idx,]
all(sero_gwas$ID==hbv_gwas$SNP) #T
all(sero_gwas$rsid==hbv_gwas$rsid,na.rm = T)
idx=match(sero_mtag$SNP,sero_gwas$rsid)
sero_gwas=sero_gwas[idx,]
hbv_gwas=hbv_gwas[idx,]

idx=which(sero_gwas$ALT==sero_mtag$A2)
sero_mtag$mtag_beta[idx]=-sero_mtag$mtag_beta[idx]
sero_mtag$mtag_z[idx]=-sero_mtag$mtag_z[idx]
tmp=sero_mtag$A2[idx]
sero_mtag$A2[idx]=sero_mtag$A1[idx]
sero_mtag$A1[idx]=tmp

idx=which(hbv_gwas$EFF==hbv_mtag$A2)
hbv_mtag$mtag_beta[idx]=-hbv_mtag$mtag_beta[idx]
hbv_mtag$mtag_z[idx]=-hbv_mtag$mtag_z[idx]
tmp=hbv_mtag$A2[idx]
hbv_mtag$A2[idx]=hbv_mtag$A1[idx]
hbv_mtag$A1[idx]=tmp

write.table(sero_mtag[,c(1:5,9,10,12)],file="../result/HBsAg_seroclearance_MTAG.csv",row.names = F,quote=F)
write.table(hbv_mtag[,c(1:5,9,10,12)],file="../result/HBVload_MTAG.csv",row.names = F,quote=F)



png(filename="../result/scater_logis_mtag_sero_pvalue.png",res=100)
ylim=range(c(-log10(sero_gwas$P),-log10(sero_mtag$mtag_pval)))
plot(-log10(sero_gwas$P),-log10(sero_mtag$mtag_pval),xlab="-log10(Logistic p-value)",ylab="-log10(MTAG pvalue)",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_mtag_sero_beta.png",res=100)
ylim=range(c(log(sero_gwas$OR),sero_mtag$mtag_beta))
plot(log(sero_gwas$OR),sero_mtag$mtag_beta,xlab="Logistic beta",ylab="MTAG beta",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_mtag_hbv_pvalue.png",res=100)
ylim=range(c(-log10(hbv_gwas$P),-log10(hbv_mtag$mtag_pval)))
plot(-log10(hbv_gwas$P),-log10(hbv_mtag$mtag_pval),xlab="-log10(Ordinal logistic p-value)",ylab="-log10(MTAG pvalue)",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_mtag_hbv_beta.png",res=100)
ylim=range(c(hbv_gwas$Beta,hbv_mtag$mtag_beta))
plot(hbv_gwas$Beta,hbv_mtag$mtag_beta,xlab="Ordinal logistic beta",ylab="MTAG beta",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()

sero_mtag=as.data.frame(fread("../result/mtag_qHBsAg_sero_trait_1.txt"))
qHBsAg_mtag=as.data.frame(fread("../result/mtag_qHBsAg_sero_trait_2.txt"))
write.table(sero_mtag[,c(1:5,9,10,12)],file="../result/HBsAg_seroclearance_mtag.csv",row.names = F,quote=F)
write.table(qHBsAg_mtag[,c(1:5,9,10,12)],file="../result/qHBsAg_mtag.csv",row.names = F,quote=F)

tmp1=fread("../result/HBsAg_seroclearance_mtag.csv")
tmp2=fread("../result/HBsAg_seroclearance_MTAG.csv")
plot(-log10(tmp1$mtag_pval),-log10(tmp2$mtag_pval))
plot(tmp1$mtag_beta,tmp2$mtag_beta)

hbeag_mtag=as.data.frame(fread("../result/mtag_hbvload_hbeag_trait_1.txt"))
sum(hbeag_mtag$mtag_pval<5e-8) #13
hbvload_mtag=as.data.frame(fread("../result/mtag_hbvload_hbeag_trait_2.txt"))
sum(hbvload_mtag$mtag_pval<5e-8) #13
write.table(hbeag_mtag[,c(1:5,9,10,12)],file="../result/HBeAg_mtag3.csv",row.names = F,quote=F)
write.table(hbvload_mtag[,c(1:5,9,10,12)],file="../result/HBV_mtag3.csv",row.names = F,quote=F)

hbeag_mtag=as.data.frame(fread("../result/mtag_hbvload_hbeag_qHBsAg_trait_2.txt"))
sum(hbeag_mtag$mtag_pval<5e-8) #22
hbvload_mtag=as.data.frame(fread("../result/mtag_hbvload_hbeag_qHBsAg_trait_1.txt"))
sum(hbvload_mtag$mtag_pval<5e-8) #22
qHBsAg_mtag=as.data.frame(fread("../result/mtag_hbvload_hbeag_qHBsAg_trait_3.txt"))
sum(qHBsAg_mtag$mtag_pval<5e-8) #22
write.table(hbeag_mtag[,c(1:5,9,10,12)],file="../result/HBeAg_mtag4.csv",row.names = F,quote=F)
write.table(hbvload_mtag[,c(1:5,9,10,12)],file="../result/HBV_mtag4.csv",row.names = F,quote=F)
write.table(qHBsAg_mtag[,c(1:5,9,10,12)],file="../result/qHBsAg_mtag4.csv",row.names = F,quote=F)
hbvload=as.data.frame(fread("../result/sumstat.cvs"))
sum(hbvload$P<5e-8) #17
hbeag=as.data.frame(fread("../result/hbeag_sumstat.csv"))
sum(hbeag$P<5e-8) #12
qhbsag=as.data.frame(fread("../result/qHBsAg_sumstat.csv"))
sum(qhbsag$P<5e-8)
