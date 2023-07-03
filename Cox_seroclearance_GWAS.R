#!/usr/bin/env Rscript
#I would like to conduct a comparison of estimates, including effect size and p-value, for logistic regression and Cox regression models when the outcome is HBsAg clearance. In the Cox regression, we will use time-in-study as the time scale, adjusted for baseline age, sex, and top 5PCs.
require(survival)
#coxph(Surv(time_entry_sc, HBsAg_seroclearance) ~  Baseline_AGECAT + SEX + PC1 + PC2 + PC3 + PC4 + PC5+SNP_of_interest, data=GLMdata)
setwd("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")

#logistic model
# plink="/usr/local/apps/plink/1.9/plink"
# plink2="/usr/local/apps/plink/2.3-alpha/plink2"
dat=read.csv("../data/GLMdata.csv")
load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData")
idx=match(rownames(pheno),dat$id_new)
all(rownames(pheno)==dat$id_new[idx])
pheno$HBsAg_seroclearance=dat$HBsAg_seroclearance[idx]
pheno$entry_date=as.numeric(as.character(dat$entry_date[idx]))
pheno$sc_date=as.numeric(as.character(dat$sc_date[idx]))
# 
# #form phenotype data
# tmp=data.frame(FID=0,IID=rownames(pheno),seroclearance=pheno$HBsAg_seroclearance)
# tmp$seroclearance=tmp$seroclearance+1
# tmp$seroclearance[which(is.na(tmp$seroclearance))]=-9
# write.table(tmp,file="../result/HBsAg_sero_pheno.txt",row.names = F,col.names = T,sep=" ",quote=F)
# #form covariate data
# tmp=data.frame(FID=0,IID=rownames(pheno),pheno[,1:7])
# levels(tmp$SEX)=c("F","M")
# levels(tmp$AGE.grp)
# levels(tmp$AGE.grp)=c("Age2","Age3","Age4","Age5U")
# colnames(tmp)[which(colnames(tmp)=="AGE.grp")]="AGE"
# write.table(tmp,file="../result/HBsAg_sero_cov.txt",row.names = F,col.names = T,sep=" ",quote=F)
# 
# prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
# cmd=paste0(plink2," --bfile ",prefix," --pheno ../result/HBsAg_sero_pheno.txt --covar  ../result/HBsAg_sero_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/HBsAg_seroclearance")
# system(cmd)
# 
# logistres=as.data.frame(fread("../result/HBsAg_seroclearance.seroclearance.glm.logistic"))
# logistres=logistres[logistres$TEST=="ADD",]
# 
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
# genome <- BSgenome.Hsapiens.UCSC.hg38
# all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
# seqlevelsStyle(genome) <- "NCBI"
# 
# ## construct a GPos object containing all the positions we're interested in
# positions <- GPos(seqnames = logistres$CHR, pos = logistres$POS)
# 
# ## query the genome with positions
# my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
# 
# ## this gives us a GPos object
# my_snps = as.data.frame(my_snps)
# 
# idx = match(paste0(logistres$CHR,":",logistres$POS),paste0(my_snps$seqnames,":",my_snps$pos))
# logistres$rsid=my_snps$RefSNP_id[idx]
# tmp=logistres[,c(3,16,1,2,4,5,9,10,14)]
# write.csv(tmp,file="../result/seroclearance_sumstat.csv",row.names = F,quote=F)

#Cox model
args = commandArgs(trailingOnly=TRUE)
SNP_TRAW_FILE = args[1]
OUT_FILE = args[2]
hbeag="all"
if(length(args)==3)
hbeag = args[3]
library(survival)
coxmodel=function(snp.traw.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/processed.traw__1.gz", out.file = "../result/cox_sero/processed.traw__1.txt",
                  HBeAg="all")
{
  print(snp.traw.file)
  snp.traw = read.delim(file = snp.traw.file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # CHR     SNP     (C)M    POS     COUNTED ALT     1000_CHB.AxiomGT1_1000_CHB.AxiomGT1
  snp.info = snp.traw[, c("CHR", "SNP", "POS", "COUNTED", "ALT")]
  colnames(snp.info) = c("CHR", "SNP", "BP", "EFF", "REF")
  
  snp.data = snp.traw[, -(1:6)]
  col.names = gsub("^0_","",colnames(snp.data))
  
  colnames(snp.data) = col.names
  idx = match(rownames(pheno), col.names)
  snp.data = t(snp.data[, idx])
  colnames(snp.data) = snp.info$SNP
  
  res = snp.info
  res$P=res$coef=res$se=NA
  
  HBeAgtable=read.table("../result/sero_hbeag_pheno.txt",header=T)
  HBeAgtable=HBeAgtable[HBeAgtable$seroclearance %in% c(1,2),]
  
  for(i in 1:ncol(snp.data)){
    cur.data = data.frame(SNP = snp.data[, i], pheno, check.names = FALSE, stringsAsFactors = FALSE)
    if (HBeAg=="negative") #HBeAg negative
    {
      cur.data=cur.data[rownames(cur.data) %in% HBeAgtable$IID,]
    }
    cur.data=cur.data[!is.na(cur.data$HBsAg_seroclearance),]
    
    m1 = tryCatch(
      expr = {
        coxph(Surv(as.numeric(entry_date),as.numeric(sc_date), HBsAg_seroclearance) ~  AGE.grp + SEX + PC1 + PC2 + PC3 + PC4 + PC5+ SNP, cur.data)
      },
      error = function(e){ 
        return(NULL)
      }
    )
    
    if (!is.null(m1))
    {
      coeff=summary(m1)$coefficients
      idx=which(rownames(coeff)=="SNP")
      if(length(idx)>0)
      {
        res$coef[i]=coeff[idx,1]
        res$se[i]=coeff[idx,3]
        res$P[i]=coeff[idx,5]
      }
    }
  }
  write.table(res, file = out.file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  #return(res)  
}
coxmodel(snp.traw.file = SNP_TRAW_FILE, out.file = OUT_FILE, HBeAg =hbeag)
print("done")
# allfiles=list.files("/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/","*.gz")
# tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/Cox_seroclearance_GWAS.R",length(allfiles)),
#                snptrawfile=paste0("/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/",allfiles),
#                outfile=paste0("../result/cox_sero/",gsub(".gz",".txt",allfiles,fixed=T)))
# write.table(tmp,file="log/cox_sero.swarm",row.names = F,col.names = F,sep="\t",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HBV_GWAS/code/log/cox_sero.swarm -g 16 --module R/4.2.0 --time=10:00:00

allfiles=list.files("/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/","*.gz")
tmp=data.frame(code=rep("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/Cox_seroclearance_GWAS.R",length(allfiles)),
               snptrawfile=paste0("/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/",allfiles),
               outfile=paste0("../result/cox_sero_hbeag_neg/",gsub(".gz",".txt",allfiles,fixed=T)),
               hbeag="negative")
write.table(tmp,file="log/cox_sero_hbeag_neg.swarm",row.names = F,col.names = F,sep="\t",quote=F)
#swarm -f /data/BB_Bioinformatics/Kevin/HBV_GWAS/code/log/cox_sero_hbeag_neg.swarm -g 16 --module R/4.2.0 --time=10:00:00