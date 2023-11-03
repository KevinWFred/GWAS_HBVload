#!/usr/bin/env Rscript
library(data.table)
plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"
qHBsAg=read.csv("../data/Baseline_qHBsAg.csv")
ann=read.csv("../data/ana_final_VL_new.csv")
tmp=intersect(qHBsAg$chip_id,ann$chip_id)
idx1=match(tmp,ann$chip_id)
idx2=match(tmp,qHBsAg$chip_id)
qHBsAg$subjectid=NA
qHBsAg$subjectid[idx2]=ann$subject_id[idx1]
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
idx=match(rownames(pheno),qHBsAg$subjectid)
all(rownames(pheno)==qHBsAg$subjectid[idx])
pheno$qHBsAg=log10(qHBsAg$baseline_HBsAg[idx])
seroclearance=read.csv("../data/id_seroclearance.csv")
idx=match(rownames(pheno),seroclearance$id_new)
pheno$sero=seroclearance$HBsAg_seroclearance[idx]
boxplot(pheno$qHBsAg~pheno$sero)
hist(pheno$qHBsAg)
#form phenotype data
tmp=data.frame(FID=0,IID=rownames(pheno),qHBsAg=pheno$qHBsAg)
#write.table(tmp,file="../result/qHBsAg_pheno.txt",row.names = F,col.names = T,sep=" ",quote=F)
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed"
cmd=paste0(plink2," --bfile ",prefix," --pheno ../result/qHBsAg_pheno.txt --covar  ../result/sero_cov.txt --glm hide-covar --ci 0.95 --out ../result/qHBsAg_maf01")
system(cmd)
res=as.data.frame(fread("../result/qHBsAg_maf01.qHBsAg.glm.linear"))
addfreq=function(dat=logistres)
{
  freqdat=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed.frq"))
  idx=match(dat$ID,freqdat$SNP)
  dat$MAF=freqdat$MAF[idx]
  return(dat)
}
res=addfreq(dat=res)
oldres=as.data.frame(fread("../result/qHBsAg.qHBsAg.glm.linear"))
allres=inner_join(oldres,res,by="ID")
table(allres$P.x==allres$P.y)
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
png("../result/qHBsAg_maf01_qqplot.png",res=100)
qqplot(res$P) #1.04
dev.off()
library(qqman)
colnames(res)[1]="CHR"
png("../result/qHBsAg_maf01_manhattan.png",res=100, width=1200)
manhattan(res, chr="CHR", bp="POS", snp="ID", p="P", suggestiveline = -log10(5e-6),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

.libPaths(c("/data/wangx53",.libPaths()))
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = res$CHR, pos = res$POS)

## query the genome with positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)
tmp=intersect(paste0(res$CHR,":",res$POS,":",res$REF,":",res$ALT),
              paste0(my_snps$seqnames,":",my_snps$pos,":",my_snps$ref_allele,":",my_snps$alt_alleles))

idx = match(paste0(res$CHR,":",res$POS),paste0(my_snps$seqnames,":",my_snps$pos))
res$rsid=my_snps$RefSNP_id[idx]
res$a1=res$ALT
res$a2=res$REF
idx=which(res$A1==res$REF)
res$a1[idx]=res$REF[idx]
res$a2[idx]=res$ALT[idx]
tmp=res[,c(3,16,1,2,17,18,9,10,14)]
tmp=res[,c(3,17,1,2,16,18,19,9,10,14)]
tmp$OR95=res$U95
tmp$OR05=res$L95
write.csv(tmp,file="../result/qHBsAg_maf01_sumstat.csv",row.names = F,quote=F)
tmp=fread("../result/qHBsAg_sumstat.csv")
tmp1=fread("../result/HBsAg_seroclearance_sumstat.csv")
cor(tmp$BETA,log(tmp1$OR))
