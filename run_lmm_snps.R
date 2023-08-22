#!/usr/bin/env Rscript
library(data.table)

plink2="/usr/local/apps/plink/2.3-alpha/plink2"
#to get the genotype data
selsnps=c("rs148385846", "rs446717", "rs440640", "rs6906021","rs7838283","rs74424600","rs2031551")
#to findout snpID
# sumstat=as.data.frame(fread("../result/sumstat.cvs"))
# idx=match(selsnps,sumstat$rsid)
# tmp=sumstat$SNP[idx]
# write.table(tmp,file="../result/lmm_selsnps.txt",row.names = F,col.names = F,quote=F)
# prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
# cmd=paste0(plink2," --bfile ",prefix," --extract ../result/lmm_selsnps.txt --recode A-transpose --out ../result/lmm_selsnps")
# system(cmd)

gt=as.data.frame(fread("../result/lmm_selsnps.traw"))
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

Sys.time()

res=data.frame(SNP=rownames(gt1),rsid=NA,P=NA,SE=NA,beta=NA)
idx=match(rownames(gt1),sumstat$SNP)
res$rsid=sumstat$rsid[idx]
idx=match(new$id_new,colnames(gt1))
set.seed(1000)
for (i in 1:nrow(res))
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
      res$P[i]=tmp[idx2,4]
      res$SE[i]=tmp[idx2,2]
      res$beta[i]=tmp[idx2,1]
    }
  }
}
write.table(res,file="../result/lmm_selsnps.txt",quote=F,row.names = F,sep="\t")
Sys.time()

