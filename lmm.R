#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
nstart=as.integer(args[1])
nend=as.integer(args[2])
outprefix=as.character(args[3])
hbeag="all"
if(length(args)==4)
  hbeag=as.character(args[4])
print(paste0(nstart,"-",nend))
#extract data
library(data.table)
# plink="/usr/local/apps/plink/1.9/plink"
# plink2="/usr/local/apps/plink/2.3-alpha/plink2"
setwd("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")
selsnps=read.csv("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/suggestiveSNPs.csv") #p<5e-5
# write.table(selsnps$SNP,file="../result/ordinalr_selsnps.txt",row.names = F,col.names = F,quote=F)
# prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
# cmd=paste0(plink2," --bfile ",prefix," --extract ../result/ordinalr_selsnps.txt --recode A-transpose --out ../result/ordinalr_selsnps")
# system(cmd)

library(ordinal)
gt=as.data.frame(fread("../result/ordinalr_selsnps.traw"))
# bim=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed.bim"))
# idx=match(gt$SNP,bim$V2)
# all(gt$COUNTED==bim$V6[idx])
gt1=gt[,7:ncol(gt)]
rownames(gt1)=gt$SNP
colnames(gt1)=gsub("0_","",colnames(gt1))
HBeAgtable=read.table("../result/sero_hbeag_pheno.txt",header=T)
HBeAgtable=HBeAgtable[HBeAgtable$seroclearance %in% c(1,2),]
new=read.csv("../data/GLMdata.csv")
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
all(new$id_new %in% rownames(pheno))
idx=match(new$id_new,rownames(pheno))
new$SEX=pheno$SEX[idx]
new$SEX<-as.factor(new$SEX)
new$measure_age_cat<-ifelse(new$measure_age<=39, 1, ifelse(new$measure_age<=49,2,ifelse(new$measure_age<=59,3,ifelse(new$measure_age<=79,4,5))))
new$measure_age_cat<-as.factor(new$measure_age_cat)
new$hbvdnacat<-as.ordered(new$hbvdnacat)
if (hbeag=="negative")
{
  new=new[new$id_new %in% HBeAgtable$IID,]
}
Sys.time()

res=data.frame(SNP=rownames(gt1)[nstart:nend],rsid=selsnps$rsid[nstart:nend],P=NA,SE=NA,beta=NA)
idx=match(new$id_new,colnames(gt1))
set.seed(1000)
for (i in nstart:nend)
{
  cat(i,'..')
  new$snp=as.numeric(unlist(gt1[i,idx]))
  m = tryCatch(
    expr = {
      clmm(hbvdnacat ~ measure_age_cat+SEX+PC1+PC2+PC3+PC4+PC5+snp+(1|id_new), data=new, Hess =
             TRUE, na.action=na.omit, nAGQ = 10, control = clmm.control(maxIter = 200, maxLineIter = 200))
    },
    error = function(e){ 
      return(NULL)
    }
  )
  if (!is.null(m))
  {
    tmp=summary(m)$coefficients
    idx2=which(rownames(tmp)=="snp")
    if (length(idx2)>0)
    {
      res$P[i-nstart+1]=tmp[idx2,4]
      res$SE[i-nstart+1]=tmp[idx2,2]
      res$beta[i-nstart+1]=tmp[idx2,1]
    }
  }
}
write.table(res,file=outprefix,quote=F,row.names = F,sep="\t")
Sys.time()
print("done")

# new=read.csv("../data/GLMdata.csv")
# pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
# load(pheno.Rdata.file)
# all(new$id_new %in% rownames(pheno))
# idx=match(new$id_new,rownames(pheno))
# new$SEX=pheno$SEX[idx]
# new$SEX<-as.factor(new$SEX)
# new$measure_age_cat<-ifelse(new$measure_age<=39, 1, ifelse(new$measure_age<=49,2,ifelse(new$measure_age<=59,3,ifelse(new$measure_age<=79,4,5))))
# new$measure_age_cat<-as.ordered(new$measure_age_cat)
# new$hbvdnacat<-as.ordered(new$hbvdnacat)
# new$rs13146927_A<-as.numeric(new$rs13146927_A)
# new$rs28780111_T<-as.numeric(new$rs28780111_T)
# new$rs28771426_A<-as.numeric(new$rs28771426_A)
# new$rs510205_G  <-as.numeric(new$rs510205_G)
# new$rs7475986_G <-as.numeric(new$rs7475986_G)
# new$rs9561485_T <-as.numeric(new$rs9561485_T)
# new$rs4410197_C <-as.numeric(new$rs4410197_C)
# new$rs1466298_T <-as.numeric(new$rs1466298_T)
# Sys.time()
# m <- clmm(hbvdnacat ~ measure_age_cat+SEX+PC1+PC2+PC3+PC4+PC5+rs13146927_A+(1|id_new), data=new, Hess =
#             TRUE, na.action=na.omit, nAGQ = 10, control = clmm.control(maxIter = 200, maxLineIter = 200))
# summary(m)
# Sys.time()

# [1] "2023-04-12 20:00:19 EDT"
# > m <- clmm(hbvdnacat ~ measure_age_cat+SEX+PC1+PC2+PC3+PC4+PC5+rs13146927_A+(1|id_new), data=new, Hess =
#               +             TRUE, na.action=na.omit, nAGQ = 10, control = clmm.control(maxIter = 200, maxLineIter = 200))
# 
# > summary(m)
# Cumulative Link Mixed Model fitted with the adaptive Gauss-Hermite 
# quadrature approximation with 10 quadrature points 
# 
# formula: hbvdnacat ~ measure_age_cat + SEX + PC1 + PC2 + PC3 + PC4 + PC5 +  
#   rs13146927_A + (1 | id_new)
# data:    new
# 
# link  threshold nobs  logLik    AIC      niter       max.grad cond.H 
# logit flexible  13025 -16062.84 32155.69 2206(28347) 1.01e-03 5.9e+04
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# id_new (Intercept) 12.32    3.51    
# Number of groups:  id_new 3240 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# measure_age_cat.L  -1.46797    0.08260 -17.772  < 2e-16 ***
#   measure_age_cat.Q  -0.11486    0.05158  -2.227  0.02595 *  
#   measure_age_cat.C  -0.04848    0.04037  -1.201  0.22985    
# SEX1                0.32070    0.14123   2.271  0.02316 *  
#   PC1                -1.43887    4.13005  -0.348  0.72755    
# PC2                11.02178    4.14836   2.657  0.00789 ** 
#   PC3                 1.45554    3.92734   0.371  0.71092    
# PC4               -10.51670    3.79678  -2.770  0.00561 ** 
#   PC5                -6.34511    4.10310  -1.546  0.12200    
# rs13146927_A        0.61518    0.14967   4.110 3.95e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Threshold coefficients:
#   Estimate Std. Error z value
# 1|2  -1.8387     0.1216 -15.119
# 2|3   0.9141     0.1215   7.521
# 3|4   2.9053     0.1253  23.186
# 4|5   4.3824     0.1299  33.745
# (58255 observations deleted due to missingness)
# > Sys.time()
# [1] "2023-04-12 20:48:35 EDT"

# Sys.time()
# m <- clmm2(hbvdnacat ~ measure_age_cat+SEX+PC1+PC2+PC3+PC4+PC5+rs13146927_A+(1|id_new), data=new, Hess =
#             TRUE, na.action=na.omit, nAGQ = 10, control = clm2.control(maxIter = 200, maxLineIter = 200,method = "Newton"))
# summary(m)
# Sys.time()
# imp<-jomo.clmm(hbvdnacat ~ measure_age_cat+SEX+PC1+PC2+PC3+PC4+PC5+rs13146927_A+(1|id_new), data=new, nburn=1000, nbetween=1000, nimp=2)
# Sys.time()
