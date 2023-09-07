#!/usr/bin/env Rscript

#My idea is to conduct association analyses using HBsAg seroclearance as an outcome, just for those 6 SNPs, and focus on individuals aged 30-39, 40-49, 50-59, and 60+, separately. This way, we will have 4 ORs for each SNP. Next, we can compare the ORs among the age groups to identify potential differences in effects.

library(data.table)
hbv_gwas=as.data.frame(fread("../result/sumstat.cvs",sep=","))

plink="/usr/local/apps/plink/1.9.0-beta4.4/plink"
selsnps=c("rs148385846", "rs446717", "rs440640", "rs6906021","rs7838283","rs74424600","rs2031551")
idx=match(selsnps,hbv_gwas$rsid)
snpid=hbv_gwas$SNP[idx]

allgt=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed.traw"))
idx=match(snpid,allgt$SNP)
selgt=allgt[idx,]
rownames(selgt)=selsnps
snp.info = selgt[, c("CHR", "SNP", "POS", "COUNTED", "ALT")]
colnames(snp.info) = c("CHR", "SNP", "BP", "EFF", "REF")
selgt=selgt[,7:ncol(selgt)]
colnames(selgt)=gsub("^0_","",colnames(selgt))

pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
all(rownames(pheno) %in% colnames(selgt))
pheno=pheno[match(colnames(selgt),rownames(pheno)),]

hbeagdat=read.csv("../data/GLMdata.csv")
tmp=hbeagdat$id_new[hbeagdat$hbeag==0]
idx=match(rownames(pheno),hbeagdat$id_new)
all(rownames(pheno)==hbeagdat$id_new[idx])
pheno$bheag=hbeagdat$hbeag[idx]
table(pheno$sero,pheno$bheag,useNA="ifany")
phenoN=pheno[pheno$bheag==0,]

orininal_lr=function(opt="all")
{
  library(MASS)
  uniq_agegrp=names(table(pheno$AGE.grp)) #"(29,39]" "(39,49]" "(49,59]" "(59,79]"
  res=list()
  for(i in 1:length(uniq_agegrp))
  {
    if (opt=="all")
    {
      mypheno=pheno[pheno$AGE.grp==uniq_agegrp[i],]
    }else
    {
      mypheno=phenoN[phenoN$AGE.grp==uniq_agegrp[i],]
    }
    
    mygt=selgt[,match(rownames(mypheno),colnames(selgt))]
    snpres=data.frame(OR=rep(NA,nrow(mygt)),P=NA,SE=NA)
    rownames(snpres)=rownames(mygt)
    for (j in 1:nrow(mygt))
    {
      mypheno$snp=unlist(mygt[j,])
      m1 = tryCatch(
        expr = {
          polr(formula = HBVDNA.grp ~ ., data = mypheno, Hess = TRUE)
        },
        error = function(e){ 
          return(NULL)
        }
      )
      if (!is.null(m1))
      {
        ctable = coef(summary(m1))
        ctable = cbind(ctable, "p_value" = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)
        snpres$OR[j]=exp(ctable["snp","Value"])
        snpres$P[j]=ctable["snp","p_value"]
        snpres$SE[j]=ctable["snp","Std. Error"]
      }
    }
    res[[i]]=snpres
  }
  return(res)
}
res=orininal_lr()
resN=orininal_lr(opt="negative")
save(res,resN,file="../result/check_candidateSNPs.RData")
res[[1]][,1:2]
#                OR           P
# rs148385846 1.6254391 0.009812330
# rs446717    0.8709710 0.075844323
# rs440640    0.8709710 0.075844323
# rs6906021   0.7807152 0.001041712
# rs7838283   1.1844467 0.031591409
# rs74424600  1.2104920 0.223396427
# rs2031551   0.7150151 0.004865881
res[[2]][,1:2]
#                OR           P
# rs148385846 1.8025461 0.002439019
# rs446717    0.8103246 0.014426388
# rs440640    0.8103246 0.014426388
# rs6906021   0.8079475 0.010445984
# rs7838283   1.0366745 0.668996671
# rs74424600  1.4650960 0.014468517
# rs2031551   0.7708683 0.033793783
res[[3]][,1:2]
#                OR           P
# rs148385846 1.6552719 0.007761907
# rs446717    0.9142188 0.300951809
# rs440640    0.9142188 0.300951809
# rs6906021   0.8426789 0.040449962
# rs7838283   1.2664949 0.005396014
# rs74424600  1.3122901 0.113680973
# rs2031551   0.6988454 0.004633582
res[[4]][,1:2]
#                OR           P
# rs148385846 1.0862821 0.7685196
# rs446717    0.8582245 0.3091392
# rs440640    0.8582245 0.3091392
# rs6906021   0.7836792 0.1096412
# rs7838283   0.9860227 0.9308202
# rs74424600  1.3009657 0.4495705
# rs2031551   0.7237060 0.1845914

plotOR=function()
{
  library(ggplot2)
  agegrp=c("30_40","40_50","50_60","60plus")
  snps=c(rep(selsnps,length(agegrp)))
  agegrps=rep(agegrp,each=length(selsnps))
  ORs=NULL
  for (i in 1:length(res))
  {
    ORs=c(ORs,res[[i]]$OR)
  }
  data <- data.frame(agegroup=agegrps,snp=snps,OR=ORs)
  # Grouped
  source("../../PRS_EASLC/code/theme_publication.R")
  pdf("../result/check_candidateSNPs_barplot.pdf",width = 12)
  ggplot(data, aes(fill=agegroup, y=OR, x=snp)) + 
    geom_bar(position = position_dodge(0.9), width=0.8,stat="identity")+theme_Publication()
  dev.off()
  
  ORNs=NULL
  for (i in 1:length(resN))
  {
    ORNs=c(ORNs,resN[[i]]$OR)
  }
  data <- data.frame(agegroup=agegrps,snp=snps,OR=ORNs)
  # Grouped
  source("../../PRS_EASLC/code/theme_publication.R")
  pdf("../result/check_candidateSNPs_hbeagN_barplot.pdf",width = 12)
  ggplot(data, aes(fill=agegroup, y=OR, x=snp)) + 
    geom_bar(position = position_dodge(0.9), width=0.8,stat="identity")+theme_Publication()
  dev.off()
  
}