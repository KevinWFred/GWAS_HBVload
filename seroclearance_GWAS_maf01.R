#/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

seroclearance=read.csv("../data/id_seroclearance.csv")
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
idx=match(rownames(pheno),seroclearance$id_new)
pheno$sero=seroclearance$HBsAg_seroclearance[idx]

#add HBV genotype
gendat=read.csv("../data/genotype.csv")
table(rownames(pheno) %in% gendat$ID)
# TRUE 
# 3240
idx=match(rownames(pheno),gendat$ID)
pheno$HBVgenotype=gendat$HBV_genotype[idx]

phenotypeB=rownames(pheno)[pheno$HBVgenotype=="B"]
phenotypeC=rownames(pheno)[pheno$HBVgenotype=="C"]

#add ALT
metadat=read.csv("../data/HBVmeta.csv")
idx=match(rownames(pheno),metadat$ID)
pheno$ALT=metadat$XALT[idx]

#form phenotype data
tmp=data.frame(FID=0,IID=rownames(pheno),seroclearance=pheno$sero)
tmp$seroclearance=tmp$seroclearance+1
tmp$seroclearance[which(is.na(tmp$seroclearance))]=-9
#write.table(tmp,file="../result/sero_pheno.txt",row.names = F,col.names = T,sep=" ",quote=F)
#form covariate data
tmp=data.frame(FID=0,IID=rownames(pheno),pheno[,1:7])
levels(tmp$SEX)=c("F","M")
levels(tmp$AGE.grp)
levels(tmp$AGE.grp)=c("Age2","Age3","Age4","Age5U")
colnames(tmp)[which(colnames(tmp)=="AGE.grp")]="AGE"
#write.table(tmp,file="../result/sero_cov.txt",row.names = F,col.names = T,sep=" ",quote=F)

tmp=data.frame(FID=0,IID=rownames(pheno),pheno$sero)
tmp=tmp[!is.na(tmp$pheno.sero),]
#write.table(tmp[,1:2],file="../result/sero_samples.txt",row.names = F,col.names = T,sep=" ",quote=F)

#Run GWAS
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic hide-covar --ci 0.95 --out ../result/seroclearance_maf01")
system(cmd)
tmp=as.data.frame(fread("../result/seroclearance_maf01.seroclearance.glm.logistic.hybrid"))
tmp=tmp[tmp$V8=="ADD",]
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/seroclearance_maf01")
system(cmd)

logistres=as.data.frame(fread("../result/seroclearance_maf01.seroclearance.glm.logistic"))

addfreq=function(dat=logistres)
{
  freqdat=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed.frq"))
  idx=match(dat$ID,freqdat$SNP)
  dat$MAF=freqdat$MAF[idx]
  return(dat)
}
logistres=addfreq()
oldlogistres=as.data.frame(fread("../result/seroclearance.seroclearance.glm.logistic"))
library(tidyverse)
allres=inner_join(oldlogistres,logistres,by="ID")
table(allres$P.x==allres$P.y)
logistres1=as.data.frame(fread("../result/seroclearance_maf01.seroclearance.glm.logistic.hybrid"))
all(logistres$ID==logistres1$ID)
idx=which(is.na(logistres$P))
logistres[idx,1:15]=logistres1[idx,-7]
#check one example
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_samples.txt --snp chr6:31297537:T:C --recode A-transpose --out test")
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_samples.txt --snp chr17:57564600:G:A --recode A-transpose --out test")
system(cmd)
tmp=as.data.frame(fread("test.traw"))
colnames(tmp)=gsub("0_","",colnames(tmp))
tmp=tmp[,7:ncol(tmp)]
all(colnames(tmp) %in% rownames(pheno))
idx=match(colnames(tmp),rownames(pheno))
boxplot(unlist(tmp[1,])~pheno$sero[idx])
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  pvalue=pvalue[!is.na(pvalue)]
  n=length(pvalue)
  
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  
  abline(0,1,lty=2)
  chisq <- qchisq(1-pvalue,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
  
}
png("../result/sero_maf01_qqplot.png",res=100)
qqplot(logistres$P) #1.01
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/sero_maf01_manhattan.png",res=100, width=1200)
manhattan(logistres[!is.na(logistres$P),], chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
logistres$rsid=my_snps$RefSNP_id[idx]
logistres$a1=logistres$ALT
logistres$a2=logistres$REF
idx=which(logistres$A1==logistres$REF)
logistres$a1[idx]=logistres$REF[idx]
logistres$a2[idx]=logistres$ALT[idx]
#tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
tmp=logistres[,c(3,17,1,2,16,18,19,9,10,14)]
#tmp=logistres[,c(3,17,1,2,18,19,10,11,15)]
colnames(tmp)[8]="BETA"
tmp$BETA=log(tmp$BETA)
colnames(tmp)[9]="SE"
tmp$OR95=logistres$U95
tmp$OR05=logistres$L95
write.csv(tmp,file="../result/HBsAg_seroclearance_maf01_sumstat.csv",row.names = F,quote=F)

load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/refgeneshg38codinggenes.RData")
topsnps=logistres[logistres$P<5e-5,]
colnames(topsnps)[1:2]=c("CHR","BP")
topsnps$gene=NA
topsnps$closestgene=NA
topsnps$gene_dist=NA
gr_topsnps=GRanges(seqnames = topsnps$CHR,ranges = IRanges(topsnps$BP,width = 1))
gr_allgenes=GRanges(seqnames = allgenes$chr,ranges = IRanges(start=allgenes$start,end=allgenes$end))
for (i in 1:nrow(topsnps))
{
  tmp=distance(gr_allgenes,gr_topsnps[i])
  idx=which.min(tmp)
  if (tmp[idx]>0)
  {
    if (allgenes$start[idx]>topsnps$BP[i])
    {
      topsnps$gene[i]=paste0(allgenes$gene[idx-1],";",allgenes$gene[idx])
    }
    if (allgenes$end[idx]<topsnps$BP[i])
    {
      topsnps$gene[i]=paste0(allgenes$gene[idx],";",allgenes$gene[idx+1])
    }
  }else
  {
    topsnps$gene[i]=allgenes$gene[idx]
  }
  topsnps$closestgene[i]=allgenes$gene[idx]
  topsnps$gene_dist[i]=min(tmp,na.rm = T)
}
head(topsnps)
tmp=data.frame(topsnps[,c(3,16,1,2,4,5,9,10,14,17,18,19)])
tmp=tmp[order(tmp$P),]
write.csv(tmp,file="../result/seroclearance_topsnps.csv",row.names = F)
tmp1=read.csv("../result/seroclearance_topsnps.csv")
# on HBeAg negative samples-----------------
seroclearance=read.csv("../data/id_seroclearance.csv")
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
idx=match(rownames(pheno),seroclearance$id_new)
pheno$sero=seroclearance$HBsAg_seroclearance[idx]
hbeagdat=read.csv("../data/GLMdata.csv")
tmp=hbeagdat$id_new[hbeagdat$hbeag==0]
idx=match(rownames(pheno),hbeagdat$id_new)
all(rownames(pheno)==hbeagdat$id_new[idx])
pheno$bheag=hbeagdat$hbeag[idx]
table(pheno$sero,pheno$bheag,useNA="ifany")
pheno=pheno[pheno$bheag==0,]

#form phenotype data
tmp=data.frame(FID=0,IID=rownames(pheno),seroclearance=pheno$sero)
tmp$seroclearance=tmp$seroclearance+1
tmp$seroclearance[which(is.na(tmp$seroclearance))]=-9
write.table(tmp,file="../result/sero_hbeag_pheno.txt",row.names = F,col.names = T,sep=" ",quote=F)
#form covariate data
tmp=data.frame(FID=0,IID=rownames(pheno),pheno[,1:7])
levels(tmp$SEX)=c("F","M")
levels(tmp$AGE.grp)
levels(tmp$AGE.grp)=c("Age2","Age3","Age4","Age5U")
colnames(tmp)[which(colnames(tmp)=="AGE.grp")]="AGE"
write.table(tmp,file="../result/sero_hbeag_cov.txt",row.names = F,col.names = T,sep=" ",quote=F)

tmp=data.frame(FID=0,IID=rownames(pheno),pheno$sero)
tmp=tmp[!is.na(tmp$pheno.sero),]
write.table(tmp[,1:2],file="../result/sero_hbeag_samples.txt",row.names = F,col.names = T,sep=" ",quote=F)

#Run GWAS
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_samples.txt --pheno ../result/sero_hbeag_pheno.txt --covar  ../result/sero_hbeag_cov.txt --logistic hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_maf01")
system(cmd)

cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_samples.txt --pheno ../result/sero_hbeag_pheno.txt --covar  ../result/sero_hbeag_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_maf01")
system(cmd)
logistres=as.data.frame(fread("../result/seroclearance_hbeag_maf01.seroclearance.glm.logistic"))
table(is.na(logistres$P))
logistres=addfreq()
logistres1=as.data.frame(fread("../result/seroclearance_hbeag_maf01.seroclearance.glm.logistic.hybrid"))
library(tidyverse)
allres=inner_join(oldlogistres,logistres,by="ID")
table(allres$P.x==allres$P.y)
all(logistres$ID==logistres1$ID)
idx=which(is.na(logistres$P))
logistres[idx,1:15]=logistres1[idx,-7]

qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  
  abline(0,1,lty=2)
  chisq <- qchisq(1-pvalue,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
  
}
png("../result/sero_hbeag_maf01_qqplot.png",res=100)
qqplot(logistres$P) #1.02
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/sero_hbeag_maf01_manhattan.png",res=100, width=1200)
manhattan(logistres, chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
logistres$rsid=my_snps$RefSNP_id[idx]
logistres$a1=logistres$ALT
logistres$a2=logistres$REF
idx=which(logistres$A1==logistres$REF)
logistres$a1[idx]=logistres$REF[idx]
logistres$a2[idx]=logistres$ALT[idx]
#tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
#tmp=logistres[,c(3,17,1,2,18,19,10,11,15)]
tmp=logistres[,c(3,17,1,2,16,18,19,9,10,14)]
colnames(tmp)[8]="BETA"
tmp$BETA=log(tmp$BETA)
colnames(tmp)[9]="SE"
tmp$OR95=logistres$U95
tmp$OR05=logistres$L95
write.csv(tmp,file="../result/HBsAg_seroclearance_HBeAg_negative_maf01_sumstat.csv",row.names = F,quote=F)

# load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/refgeneshg38codinggenes.RData")
# topsnps=logistres[logistres$P<5e-5,]
# colnames(topsnps)[1:2]=c("CHR","BP")
# topsnps$gene=NA
# topsnps$closestgene=NA
# topsnps$gene_dist=NA
# gr_topsnps=GRanges(seqnames = topsnps$CHR,ranges = IRanges(topsnps$BP,width = 1))
# gr_allgenes=GRanges(seqnames = allgenes$chr,ranges = IRanges(start=allgenes$start,end=allgenes$end))
# for (i in 1:nrow(topsnps))
# {
#   tmp=distance(gr_allgenes,gr_topsnps[i])
#   idx=which.min(tmp)
#   if (tmp[idx]>0)
#   {
#     if (allgenes$start[idx]>topsnps$BP[i])
#     {
#       topsnps$gene[i]=paste0(allgenes$gene[idx-1],";",allgenes$gene[idx])
#     }
#     if (allgenes$end[idx]<topsnps$BP[i])
#     {
#       topsnps$gene[i]=paste0(allgenes$gene[idx],";",allgenes$gene[idx+1])
#     }
#   }else
#   {
#     topsnps$gene[i]=allgenes$gene[idx]
#   }
#   topsnps$closestgene[i]=allgenes$gene[idx]
#   topsnps$gene_dist[i]=min(tmp,na.rm = T)
# }
# head(topsnps)
# tmp=data.frame(topsnps[,c(3,16,1,2,4,5,9,10,14,17,18,19)])
# tmp=tmp[order(tmp$P),]
# write.csv(tmp,file="../result/seroclearance_hbeag_topsnps.csv",row.names = F)
# tmp1=as.data.frame(fread("../result/seroclearance_hbeag_topsnps.csv"))

#hbeag negative, HBVgenotype B
tmp=data.frame(FID=0,IID=rownames(pheno),pheno$sero)
tmp=tmp[!is.na(tmp$pheno.sero),]
tmp=tmp[tmp$IID %in% phenotypeB,]
write.table(tmp[,1:2],file="../result/sero_hbeag_neg_typeB_samples.txt",row.names = F,col.names = T,sep=" ",quote=F)

#Run GWAS 124 cases and 1061 controls
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_neg_typeB_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_negative_typeB_maf01")
system(cmd)
tmp=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeB_maf01.seroclearance.glm.logistic.hybrid"))
tmp=tmp[tmp$V8=="ADD",]
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_neg_typeB_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_negative_typeB_maf01")
system(cmd)

logistres=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeB_maf01.seroclearance.glm.logistic"))

logistres=addfreq()
# oldlogistres=as.data.frame(fread("../result/seroclearance.seroclearance.glm.logistic"))
# library(tidyverse)
# allres=inner_join(oldlogistres,logistres,by="ID")
# table(allres$P.x==allres$P.y)
logistres1=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeB_maf01.seroclearance.glm.logistic.hybrid"))
all(logistres$ID==logistres1$ID)
idx=which(is.na(logistres$P))
logistres[idx,1:15]=logistres1[idx,-7]

png("../result/seroclearance_hbeag_negative_typeB_maf01_qqplot.png",res=100)
qqplot(logistres$P) #1.03
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/seroclearance_hbeag_negative_typeB_maf01_manhattan.png",res=100, width=1200)
manhattan(logistres[!is.na(logistres$P),], chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
logistres$rsid=my_snps$RefSNP_id[idx]
logistres$a1=logistres$ALT
logistres$a2=logistres$REF
idx=which(logistres$A1==logistres$REF)
logistres$a1[idx]=logistres$REF[idx]
logistres$a2[idx]=logistres$ALT[idx]
#tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
tmp=logistres[,c(3,17,1,2,16,18,19,9,10,14)]
#tmp=logistres[,c(3,17,1,2,18,19,10,11,15)]
colnames(tmp)[8]="BETA"
tmp$BETA=log(tmp$BETA)
colnames(tmp)[9]="SE"
tmp$OR95=logistres$U95
tmp$OR05=logistres$L95
write.csv(tmp,file="../result/HBsAg_seroclearance_hbeag_negative_typeB_maf01_sumstat.csv",row.names = F,quote=F)


#hbeag negative, HBVgenotype C
tmp=data.frame(FID=0,IID=rownames(pheno),pheno$sero)
tmp=tmp[!is.na(tmp$pheno.sero),]
tmp=tmp[tmp$IID %in% phenotypeC,]
write.table(tmp[,1:2],file="../result/sero_hbeag_neg_typeC_samples.txt",row.names = F,col.names = T,sep=" ",quote=F)

#Run GWAS 116 cases and 670 controls
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_neg_typeC_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_negative_typeC_maf01")
system(cmd)
tmp=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeB_maf01.seroclearance.glm.logistic.hybrid"))
tmp=tmp[tmp$V8=="ADD",]
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_hbeag_neg_typeC_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/seroclearance_hbeag_negative_typeC_maf01")
system(cmd)

logistres=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeC_maf01.seroclearance.glm.logistic"))

logistres=addfreq()
# oldlogistres=as.data.frame(fread("../result/seroclearance.seroclearance.glm.logistic"))
# library(tidyverse)
# allres=inner_join(oldlogistres,logistres,by="ID")
# table(allres$P.x==allres$P.y)
logistres1=as.data.frame(fread("../result/seroclearance_hbeag_negative_typeC_maf01.seroclearance.glm.logistic.hybrid"))
all(logistres$ID==logistres1$ID)
idx=which(is.na(logistres$P))
logistres[idx,1:15]=logistres1[idx,-7]

png("../result/seroclearance_hbeag_negative_typeC_maf01_qqplot.png",res=100)
qqplot(logistres$P) #1.03
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/seroclearance_hbeag_negative_typeC_maf01_manhattan.png",res=100, width=1200)
manhattan(logistres[!is.na(logistres$P),], chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
logistres$rsid=my_snps$RefSNP_id[idx]
logistres$a1=logistres$ALT
logistres$a2=logistres$REF
idx=which(logistres$A1==logistres$REF)
logistres$a1[idx]=logistres$REF[idx]
logistres$a2[idx]=logistres$ALT[idx]
#tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
tmp=logistres[,c(3,17,1,2,16,18,19,9,10,14)]
#tmp=logistres[,c(3,17,1,2,18,19,10,11,15)]
colnames(tmp)[8]="BETA"
tmp$BETA=log(tmp$BETA)
colnames(tmp)[9]="SE"
tmp$OR95=logistres$U95
tmp$OR05=logistres$L95
write.csv(tmp,file="../result/HBsAg_seroclearance_hbeag_negative_typeC_maf01_sumstat.csv",row.names = F,quote=F)

#type and alt
#form sample file
tmp=data.frame(FID=0,IID=rownames(pheno),pheno$sero)
tmp=tmp[!is.na(tmp$pheno.sero),]
tmp=tmp[tmp$IID %in% c(phenotypeB,phenotypeC),] #only include typeB and typeC
write.table(tmp[,1:2],file="../result/sero_type_alt_samples.txt",row.names = F,col.names = T,sep=" ",quote=F)

#form covariate data
tmp=data.frame(FID=0,IID=rownames(pheno),pheno[,c(1:7,10,11)])
levels(tmp$SEX)=c("F","M")
levels(tmp$AGE.grp)
levels(tmp$AGE.grp)=c("Age2","Age3","Age4","Age5U")
colnames(tmp)[which(colnames(tmp)=="AGE.grp")]="AGE"
write.table(tmp,file="../result/sero_type_alt_cov.txt",row.names = F,col.names = T,sep=" ",quote=F)

#Run GWAS 243 cases and 1893 controls
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_type_alt_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_type_alt_cov.txt --logistic hide-covar --covar-variance-standardize --ci 0.95 --out ../result/seroclearance_type_alt_maf01")
system(cmd)
cmd=paste0(plink2," --bfile ",prefix," --keep ../result/sero_type_alt_samples.txt --pheno ../result/sero_pheno.txt --covar  ../result/sero_type_alt_cov.txt --logistic no-firth hide-covar --covar-variance-standardize --ci 0.95 --out ../result/seroclearance_type_alt_maf01")
system(cmd)

logistres=as.data.frame(fread("../result/seroclearance_type_alt_maf01.seroclearance.glm.logistic"))

logistres=addfreq()
# oldlogistres=as.data.frame(fread("../result/seroclearance.seroclearance.glm.logistic"))
# library(tidyverse)
# allres=inner_join(oldlogistres,logistres,by="ID")
# table(allres$P.x==allres$P.y)
logistres1=as.data.frame(fread("../result/seroclearance_type_alt_maf01.seroclearance.glm.logistic.hybrid"))
all(logistres$ID==logistres1$ID)
idx=which(is.na(logistres$P))
logistres[idx,1:15]=logistres1[idx,-7]

png("../result/seroclearance_type_alt_maf01_qqplot.png",res=100)
qqplot(logistres$P) # 1.01
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/seroclearance_type_alt_maf01_manhattan.png",res=100, width=1200)
manhattan(logistres[!is.na(logistres$P),], chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
logistres$rsid=my_snps$RefSNP_id[idx]
logistres$a1=logistres$ALT
logistres$a2=logistres$REF
idx=which(logistres$A1==logistres$REF)
logistres$a1[idx]=logistres$REF[idx]
logistres$a2[idx]=logistres$ALT[idx]
#tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
tmp=logistres[,c(3,17,1,2,16,18,19,9,10,14)]
#tmp=logistres[,c(3,17,1,2,18,19,10,11,15)]
colnames(tmp)[8]="BETA"
tmp$BETA=log(tmp$BETA)
colnames(tmp)[9]="SE"
tmp$OR95=logistres$U95
tmp$OR05=logistres$L95
write.csv(tmp,file="../result/HBsAg_seroclearance_type_alt_maf01_sumstat.csv",row.names = F,quote=F)
