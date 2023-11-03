#!/usr/bin/env Rscript
#read and visualize GWAS results
#check imputed data
#setwd("/data/DCEGLeiSongData/Kevin/HBVloadGwas/code")
library(data.table)

resultfolder = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result_maf01/"
allfiles = list.files(resultfolder,pattern = "processed.*.txt")
alljobs = unlist(strsplit(allfiles,"__"))
alljobs = alljobs[seq(2,length(alljobs),2)]
alljobs = unlist(strsplit(alljobs,".",fixed=T))
alljobs = as.numeric(alljobs[seq(1,length(alljobs),2)])
if (length(unique(alljobs))!=max(alljobs)) warning("Some results are missing")

allres = NULL
i=1
for (myfile in paste0(resultfolder,allfiles))
{
  if (i %% 100==0) cat(i,'..')
  tmp = fread(myfile,header = T)
  allres = rbind(allres,tmp)
  i = i+1
}

idx = order(allres$CHR,allres$BP)
allres = allres[idx,]

addfreq=function(dat=allres)
{
  freqdat=as.data.frame(fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed_new/processed.frq"))
  idx=match(dat$SNP,freqdat$SNP)
  dat$MAF=freqdat$MAF[idx]
  return(dat)
}
allres=addfreq()
#write.table(allres,file="../result/ordinal_lr_maf01_result.txt",row.names = F,sep="\t",quote=F)
allres=as.data.frame(fread("../result/ordinal_lr_maf01_result.txt",header=T))
oldres=as.data.frame(fread("../result/ordinal_lr_result.txt",header = T))
sig_res_old=oldres[oldres$P<5e-8,]
quantile(allres$P)
# 0%          25%          50%          75%         100% 
# 5.232458e-10 2.405516e-01 4.917590e-01 7.458708e-01 9.999999e-01
allres[which.min(allres$P),] #this is due to 10pcs used
# CHR               SNP       BP EFF REF MAF.cases MAF.ctrls       OR
# 2937630   6 chr6:31347583:G:A 31347583   A   G 0.2992223 0.2434128 1.375967
# OR_CI_%95_Low OR_CI_%95_High         SE            P
# 2937630      1.246876       1.518423 0.05026389 2.158462e-10

findsnpingzfiles=function(snp="chr6:31347583:G:A",infolder="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/")
{
   allfiles=list.files(infolder,"*.gz")
   for (i in 1:length(allfiles))
   {
     tmp=as.data.frame(fread(paste0(infolder,allfiles[i])))
     if (sum(tmp$SNP==snp)>0)
     {
       print(allfiles[i]) #processed.traw__216.gz
       break()
     }
   }
   return(allfiles[i])
}

# #processed.traw__294.gz 
# findsnpingzfiles(snp="chr6:31347583:G:A",infolder="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited_maf01/") 
# snp.traw.file="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited/processed.traw__216.gz"
# snp.traw.file="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited_maf01/processed.traw__294.gz"
# 
# snp="chr6:31347583:G:A"
# pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
# load(pheno.Rdata.file)
# 
# snp.traw = read.delim(file = snp.traw.file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# # CHR     SNP     (C)M    POS     COUNTED ALT     1000_CHB.AxiomGT1_1000_CHB.AxiomGT1
# snp.info = snp.traw[, c("CHR", "SNP", "POS", "COUNTED", "ALT")]
# colnames(snp.info) = c("CHR", "SNP", "BP", "EFF", "REF")
# which(snp.info$SNP==snp) #7948,7630
# snp.data = snp.traw[, -(1:6)]
# col.names = gsub("^0_","",colnames(snp.data))
# colnames(snp.data) = col.names
# idx = match(rownames(pheno), col.names)
# snp.data = t(snp.data[, idx])
# colnames(snp.data) = snp.info$SNP
# library(MASS)
# i=7948
# i=7630
# col.names.res = c("MAF.cases", "MAF.ctrls", "OR", "OR_CI_%95_Low", "OR_CI_%95_High", "SE", "P")
# res = matrix(NA, 1, length(col.names.res))
# colnames(res) = col.names.res
# cur.data = data.frame(SNP = snp.data[, i], pheno, check.names = FALSE, stringsAsFactors = FALSE)
# idx.cases = cur.data$HBVDNA.grp != 1
# x = cur.data[idx.cases, "SNP"]
# res[1, "MAF.cases"] = sum(x, na.rm = TRUE)/(2*sum(!is.na(x)))
# 
# idx.ctrls = !idx.cases
# x = cur.data[idx.ctrls, "SNP"]
# res[1, "MAF.ctrls"] = sum(x, na.rm = TRUE)/(2*sum(!is.na(x)))
# m1 = tryCatch(
#   expr = {
#     polr(formula = HBVDNA.grp ~ ., data = cur.data, Hess = TRUE)
#   },
#   error = function(e){ 
#     return(NULL)
#   }
# )
# ctable = coef(summary(m1))
# ctable = cbind(ctable, "p_value" = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)
# res[1, "SE"] = ctable["SNP", "Std. Error"]
# res[1, "P"] = ctable["SNP", "p_value"]
# or = exp(cbind(OR = coef(m1), ci = confint.default(m1)))
# res[1, c("OR", "OR_CI_%95_Low", "OR_CI_%95_High")] = or["SNP", ]

# MAF.cases MAF.ctrls       OR OR_CI_%95_Low OR_CI_%95_High         SE            P
# [1,] 0.2992223 0.2434128 1.364149      1.236833       1.504571 0.04998912 5.232458e-10

sum(allres$P<5e-8) #17
sum(allres$P<0.05/nrow(allres)) #8-->7
idx = which(allres$P<5e-8)
sig_res=allres[idx,]
table(sig_res$SNP %in% oldres$SNP)
table(sig_res$SNP %in% sig_res_old$SNP)
table(sig_res$SNP==sig_res_old$SNP)
which(sig_res$OR!=sig_res_old$OR)
idx=match(sig_res$SNP,oldres$SNP)
all(sig_res$EFF==oldres$EFF[idx])
all(sig_res$MAF.ctrls==oldres$MAF.ctrls[idx])
cor(sig_res$OR,oldres$OR[idx])
table(sig_res_old$SNP %in% sig_res$SNP)

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"
positions <- GPos(seqnames = allres$CHR, pos = allres$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
my_snps = as.data.frame(my_snps)
tmp1=paste0(allres$CHR,"_",allres$BP)
tmp2=paste0(my_snps$seqnames,"_",my_snps$pos)
idx = match(tmp1,tmp2)
allres$rsid=my_snps$RefSNP_id[idx]
load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/info.RData")
idx=match(allres$SNP,allinfo$SNP)
allres$Genotyped=1-as.numeric(allinfo$Genotyped[idx])
allres$BETA=log(allres$OR)
colnames(allres)[c(9,10)]=c("OR05","OR95")
colnames(allres)[c(4,5)]=c("a1","a2")
allres$ID=allres$SNP
allres1=allres[,c("ID","rsid","CHR","BP","MAF","a1","a2","BETA","OR05","OR95","SE","P","Genotyped")]
fwrite(allres1,file="../result/HBV_GWAS_maf01_sumstat.csv",row.names = F) #

#check info
impfolder = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/"
allinfo = NULL
for (i in 1:22)
{
  cat(i,'..')
  tmp = fread(input=paste0(impfolder,"chr",i,".info.gz"))
  # tmp1 = unlist(strsplit(tmp$SNP,":"))
  # tmp$BP = tmp1[seq(2,length(tmp1),4)]
  # tmp1 = tmp[tmp$MAF>0.05 & tmp$Rsq>0.3,]
  # tmp2 = unlist(strsplit(tmp1$SNP,":"))
  # tmp1$BP = tmp2[seq(2,length(tmp2),4)]
  tmp = tmp[,c(1,5,7,8)]
  tmp$Genotyped[which(tmp$Genotyped=="Imputed")]=1
  tmp$Genotyped[which(tmp$Genotyped=="Genotyped")]=0
  allinfo = rbind(allinfo,tmp)
}

#save(allinfo,allres,file="../result/info_res.RData")
load("../result/info_res.RData")

sum(allinfo$MAF>0.05)

sum(allinfo$Rsq>0.3)
sum(allinfo$MAF>0.05 & allinfo$Rsq>0.3)

#annotate the significant res
library("biomaRt")
#use GRCH38
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://useast.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
snpmart = useEnsembl(biomart="snp",host="https://useast.ensembl.org",dataset="hsapiens_snp")
#add rsid
sig_res$rsid=NA
for (i in 1:nrow(sig_res))
{
  cat(i,'..')
  myallels = c(sig_res$EFF[i],sig_res$REF[i])
  tmp=getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(sig_res$CHR[i],sig_res$BP[i],sig_res$BP[i]),mart = snpmart)
  if (nrow(tmp)>0)
  {
    for (j in 1:nrow(tmp))
    {
      alllels=unlist(strsplit(tmp$allele[j],"/",fixed=T))
      if (all(myallels %in% alllels))
      {
        sig_res$rsid[i]=tmp$refsnp_id[j]
        break
      }
    }  
  }
}
sig_res$gene=NA
for (i in 1:nrow(sig_res))
{
  cat(i,'..')

  tmp=getBM(attributes=c('hgnc_symbol','refseq_mrna'), filters = c('chromosome_name',"start","end"), values =list(sig_res$CHR[i],sig_res$BP[i],sig_res$BP[i]),mart = mart)
  if (nrow(tmp)>0)
  {
    #print(tmp$hgnc_symbol)
    sig_res$gene[i]=unique(tmp$hgnc_symbol)
  }
}  

sig_res0=sig_res[,c(2,1,3,4,5,13,12,8,14)]
#                    SNP CHR       BP EFF REF       rsid            P        OR      gene
# 1:   chr6:31297537:T:C   6 31297537   C   T  rs2524123 3.899683e-08 1.2860103 LINC02571
# 2:   chr6:31340108:G:A   6 31340108   A   G  rs1634734 2.234796e-08 0.7722487      <NA>
# 3:   chr6:31344281:C:G   6 31344281   G   C  rs2394977 2.501494e-08 0.7730776      <NA>
# 4:   chr6:31347583:G:A   6 31347583   A   G rs28771426 5.232458e-10 1.3641491      <NA>
# 5:   chr6:31348083:A:G   6 31348083   G   A  rs2596521 4.394575e-08 1.2863100      <NA>
# 6:   chr6:31351614:A:C   6 31351614   A   C  rs9266063 7.721952e-09 1.3169083      <NA>
# 7:   chr6:31351920:A:G   6 31351920   G   A  rs9391847 1.573517e-08 1.2963025      <NA>
# 8:   chr6:32591896:T:G   6 32591896   G   T rs34434863 3.353909e-08 0.7030943      <NA>
# 9:   chr6:32616916:C:G   6 32616916   G   C   rs510205 1.610256e-09 0.6873741      <NA>
# 10:   chr6:32662585:C:G   6 32662585   C   G  rs9274172 4.427464e-09 0.7656533  HLA-DQB1
# 11:   chr6:32668956:G:A   6 32668956   G   A  rs3135000 1.405080e-08 0.7744863      <NA>
# 12:   chr6:32686725:A:G   6 32686725   G   A  rs9275183 3.295882e-09 0.6933525      <NA>
# 13:   chr6:32698396:A:G   6 32698396   G   A  rs9275318 4.228326e-09 0.6943827      <NA>
# 14:   chr6:32698518:A:G   6 32698518   G   A  rs9275319 4.228326e-09 0.6943827      <NA>
# 15:   chr6:32700634:G:A   6 32700634   A   G  rs9275373 1.825946e-08 0.6137521      <NA>
# 16:  chr6:32700658:A:AG   6 32700658  AG   A       <NA> 3.368788e-09 0.6924236      <NA>
# 17: chr6:32704570:AAT:A   6 32704570   A AAT       <NA> 4.380638e-08 0.7192649      <NA>
sig_res0$P=format(sig_res0$P,digits=3)
sig_res0$OR=format(sig_res0$OR,digits=3)
rownames(sig_res0)=NULL
sig_res0[,1:8]
# SNP CHR       BP EFF REF       rsid        P    OR
# 1    chr6:31297537:T:C   6 31297537   C   T  rs2524123 3.90e-08 1.286
# 2    chr6:31340108:G:A   6 31340108   A   G  rs1634734 2.23e-08 0.772
# 3    chr6:31344281:C:G   6 31344281   G   C  rs2394977 2.50e-08 0.773
# 4    chr6:31347583:G:A   6 31347583   A   G rs28771426 5.23e-10 1.364
# 5    chr6:31348083:A:G   6 31348083   G   A  rs2596521 4.39e-08 1.286
# 6    chr6:31351614:A:C   6 31351614   A   C  rs9266063 7.72e-09 1.317
# 7    chr6:31351920:A:G   6 31351920   G   A  rs9391847 1.57e-08 1.296
# 8    chr6:32591896:T:G   6 32591896   G   T rs34434863 3.35e-08 0.703
# 9    chr6:32616916:C:G   6 32616916   G   C   rs510205 1.61e-09 0.687
# 10   chr6:32662585:C:G   6 32662585   C   G  rs9274172 4.43e-09 0.766
# 11   chr6:32668956:G:A   6 32668956   G   A  rs3135000 1.41e-08 0.774
# 12   chr6:32686725:A:G   6 32686725   G   A  rs9275183 3.30e-09 0.693
# 13   chr6:32698396:A:G   6 32698396   G   A  rs9275318 4.23e-09 0.694
# 14   chr6:32698518:A:G   6 32698518   G   A  rs9275319 4.23e-09 0.694
# 15   chr6:32700634:G:A   6 32700634   A   G  rs9275373 1.83e-08 0.614
# 16  chr6:32700658:A:AG   6 32700658  AG   A       <NA> 3.37e-09 0.692
# 17 chr6:32704570:AAT:A   6 32704570   A AAT       <NA> 4.38e-08 0.719

chr6info = fread(input=paste0(impfolder,"chr",6,".info.gz"))
idx = match(sig_res$SNP,chr6info$SNP)
sig_res$Genotyped=chr6info$Genotyped[idx]
#for zoomlocus plot
tmp1=min(sig_res$BP)
tmp2=max(sig_res$BP)
idx=which(allres$CHR==6 & allres$BP>tmp1-2000000 & allres$BP<tmp2+2000000)
sig_res2=allres[idx,]
#sig_res1=sig_res1[sig_res1$P<0.1,]
sig_res2$rsid=NA
dbsnp=fread(input="../data/dbsnp155_6_29298314_34704290")
for (i in 1:nrow(sig_res2))
{
  myallels = c(sig_res2$EFF[i],sig_res2$REF[i])
  idx = which(dbsnp$chromEnd==sig_res2$BP[i])
  if (length(idx)>0)
  for (j in idx)
  {
    alllels = c(dbsnp$ref[j],unlist(strsplit(dbsnp$alts[j],","))) 
    if (all(myallels %in% alllels))
    {
      sig_res2$rsid[i] = dbsnp$name[j]
      break
    }
  }
}

# too slow
# for (i in 1:nrow(sig_res2))
# {
#   if (i %%1000==0) 
#   {
#     cat(i,'..')
#     save(sig_res2,file="../result/sig_res2all.RData")
#   }  
#     
#   myallels = c(sig_res2$EFF[i],sig_res2$REF[i])
#   tmp=getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_end'), filters = c('chr_name',"start","end"), values =list(sig_res2$CHR[i],sig_res2$BP[i],sig_res2$BP[i]),mart = snpmart)
#   if (nrow(tmp)>0)
#   {
#     for (j in 1:nrow(tmp))
#     {
#       alllels=unlist(strsplit(tmp$allele[j],"/",fixed=T))
#       if (all(myallels %in% alllels))
#       {
#         sig_res2$rsid[i]=tmp$refsnp_id[j]
#         break
#       }
#     }  
#   }
# }
save(sig_res2,file="../result/sig_res2all.RData")
write.table(sig_res2,file="../result/sig_res2.txt",row.names = F,quote=F,sep="\t")
#https://statgen.github.io/localzoom/
#cat ../result/sig_res2.txt|bgzip -c > ../result/summstats.sorted.tab.gz && tabix -s1 -b 3 -e3 --skip-lines 1 -f ../result/summstats.sorted.tab.gz

#zoomlocus plot
#module load locuszoom/1.3

#change coordinate from hg38 to hg19
# library(rtracklayer)
# library(GenomicRanges)
# chain=import.chain("../tools/hg38ToHg19.over.chain")
# gr_dat=GRanges(seqnames = paste0("chr",sig_res1$CHR),ranges=IRanges(start=sig_res1$BP,width=1))
# tmp=liftOver(gr_dat,chain)
# sig_res1$hg19BP=NA
# for (i in 1:length(tmp))
# {
#   tmp1=unlist(tmp[i])
#   if (length(tmp1)>0)
#   {
#     if (length(tmp1)==1)
#     {
#       sig_res1$hg19BP[i]=start(tmp1)
#     }else
#     {
#       warning(paste0(i," snp has problem"))
#     }
#   }
# }
sig_res1 = sig_res2[!is.na(sig_res2$rsid),]
idx = which(duplicated(sig_res1$rsid))
idx1= which(sig_res1$rsid %In% sig_res1$rsid[idx])
sig_res1 = sig_res1[-idx1,]
METALdat=data.frame(MarkerName=sig_res1$rsid,Pvalue=sig_res1$P,ref=sig_res1$REF,
                    alt=sig_res1$EFF,stringsAsFactors = F)
METALfile=paste0("../result/sig_snps_metalfile.txt")
write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
title=""
cmd=paste0("locuszoom"," --metal ",METALfile," --flank 2000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs28771426 --pvalcol Pvalue --plotonly --prefix zoomlocus --ld-vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz ","ymax=",10," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
cmd
cmd=paste0("locuszoom"," --metal ",METALfile," --flank 2000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs28780111 --pvalcol Pvalue --plotonly --prefix zoomlocus --ld-vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz ","ymax=",10," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
cmd
cmd=paste0("locuszoom"," --metal ",METALfile," --flank 2000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs510205 --pvalcol Pvalue --plotonly --prefix zoomlocus --ld-vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz ","ymax=",10," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
cmd

# library(GenomicRanges)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(org.Hs.eg.db)
# library(annotate)
# library(stringr)
# 
# 
# # Make an header for the data.frame "myIntervals" containing 
# # the coordinates files
# header <- c("Chromosome", "Start", "End")
# myIntervals <- read.table("coordinates.txt", header = TRUE, sep = "\t")
# 
# # The function "makeGRangesFromDataFrame" from the library 
# # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
# intervals = GenomicRanges::makeGRangesFromDataFrame(myIntervals)
# txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# 
# #previous results based on genotyped data
# oldres=fread(input="../result/HLA_GWAS.txt.gz")
# library(rtracklayer)
# library(GenomicRanges)
# chain=import.chain("../tools/hg38ToHg19.over.chain")
# #for the genotyped oe
# gr_dat=GRanges(seqnames = "chr6",ranges=IRanges(start=31297537,width=1))
# tmp=liftOver(gr_dat,chain)
# start(tmp) #31265314, this is the same one in the old result

qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
  
{
  
  n=length(pvalue)
  
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  
  abline(0,1,lty=2)
  
}
png("../result/qqplot.png",res=100)
qqplot(allres$P)
dev.off()
#inflation factor
chisq <- qchisq(1-allres$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda

# chisq = qchisq(allres$P,1,lower.tail=FALSE);
# lambda <- median(chisq) / qchisq(0.5,1)

library(qqman)
png("../result/manhattan.png",res=100, width=1200)
manhattan(allres, chr="CHR", bp="BP", snp="SNP", p="P", suggestiveline = -log10(5e-5),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
ensembl <- useEnsembl(biomart = "snp")
datasets <- listDatasets(ensembl)

filters = listFilters(ensembl)

#3 snps in chr6
idx1 = which(allres$P<1e-7 & allres$CHR==6)
tmp1=allres[idx1,]
tmp1_1=tmp1[1:30,]
tmp1_2=tmp1[31:99,]
topsnps=NULL
topsnps=rbind(topsnps,tmp1[1,])
topsnps=rbind(topsnps,tmp1_1[which.min(tmp1_1$P),])
topsnps=rbind(topsnps,tmp1_2[which.min(tmp1_2$P),])

idx = which(allres$CHR==4)
tmp=allres[idx,]
topsnps=rbind(topsnps,tmp[which.min(tmp$P),])

idx = which(allres$CHR==10)
tmp=allres[idx,]
topsnps=rbind(topsnps,tmp[which.min(tmp$P),])

idx = which(allres$CHR==13)
tmp=allres[idx,]
topsnps=rbind(topsnps,tmp[which.min(tmp$P),])

idx = which(allres$CHR==19)
tmp=allres[idx,]
topsnps=rbind(topsnps,tmp[which.min(tmp$P),])
# 
# 
# #add rsid
# topsnps$rsid=NA
# for (i in c(1:7,9:nrow(topsnps)))
# {
#   cat(i,'..')
#   myallels = c(topsnps$EFF[i],topsnps$REF[i])
#   tmp = tryCatch(
#     expr = {
#       getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_end'), filters = c('chr_name',"start","end"), values =list(topsnps$CHR[i],topsnps$BP[i],topsnps$BP[i]),mart = snpmart)
#     },
#     error = function(e){ 
#       return(NULL)
#     }
#   )
#   #tmp=getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_end'), filters = c('chr_name',"start","end"), values =list(topsnps$CHR[i],topsnps$BP[i],topsnps$BP[i]),mart = snpmart)
#   if (length(tmp)>0)
#   {
#     for (j in 1:nrow(tmp))
#     {
#       alllels=unlist(strsplit(tmp$allele[j],"/",fixed=T))
#       if (all(myallels %in% alllels))
#       {
#         topsnps$rsid[i]=tmp$refsnp_id[j]
#         break
#       }
#     }  
#   }
# }
# topsnps$gene=NA
# for (i in 1:nrow(topsnps))
# {
#   cat(i,'..')
#   
#   tmp=getBM(attributes=c('hgnc_symbol','refseq_mrna'), filters = c('chromosome_name',"start","end"), values =list(topsnps$CHR[i],topsnps$BP[i],topsnps$BP[i]),mart = mart)
#   if (nrow(tmp)>0)
#   {
#     print(tmp$hgnc_symbol)
#     topsnps$gene[i]=paste0(unique(tmp$hgnc_symbol),collapse = ";")
#   }
# } 

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
# library(annotate)
# library(stringr)
# myIntervals = topsnps
# intervals = GRanges(seqnames = paste0("chr",myIntervals$CHR),ranges = IRanges(start=myIntervals$BP,width = 1))
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = genes(txdb)
# annotateIntervals <- function(intervals, txdb)
#   {
#     stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
#     anno = genes(txdb)
#     olaps = findOverlaps(intervals, anno)
#     mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
#     intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
#     intervals$gene_id =S4Vectors::splitAsList(mcols(olaps)$gene_id, intervals_factor)
#     intervals
#   }        
# 
# myAnnotation <- as.data.frame(annotateIntervals(intervals, txdb))
# myDf_master <- data.frame()
# 
# # Now we want Hugo gene names in our annotations! 
# #So, for each annotated interval get hugo gene names...
# for (i in 1:length(myAnnotation$gene_id)) {
#   # if the gene list is not empty...
#   if(length(c(na.omit(myAnnotation$gene_id[i])[[1]])) != 0) {
#     # annotate the interval and copy into a myDf data.frame
#     myDf <- data.frame(myAnnotation$seqnames[i], myAnnotation$start[i], 
#                        myAnnotation$end[i], toString(unname(getSYMBOL(c(na.omit(myAnnotation$gene_id[i])[[1]]), data='org.Hs.eg'))))
#     # append tge myDF annotations with rbind into the myDf_master
#     myDf_master <- rbind(myDf_master, myDf)
#   }
# }

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
# # extract the list of all gene_Ids from the txdb object
# genes = as.data.frame(genes(txdb))
# genes = genes[genes$seqnames %in% c(paste0("chr",c(1:22,"X","Y"))),]
# genes$seqnames = factor(genes$seqnames,levels = paste0("chr",c(1:22,"X","Y")))
# idx=order(genes$seqnames,genes$start)
# genes=genes[idx,]
# genes$symbol=NA
# for (i in 1:nrow(genes))
# {
#   genes$symbol[i]=getSYMBOL(as.character(genes$gene_id[i]),data='org.Hs.eg')
# }

library(data.table)
genes=fread("grep -v '^#' ../../tools/GRCh38_latest_genomic.gff" ,sep="\t",fill=T)
genes=genes[genes$V3=="gene",]
allgenes=data.frame(matrix(NA,nrow=nrow(genes),ncol=7))
colnames(allgenes)=c("NC","gene","chr","start","end","geneid","type")
allgenes$NC=genes$V1
allgenes$start=genes$V4
allgenes$end=genes$V5
for ( i in 1:nrow(allgenes))
{
  if (i %%2000==0) cat(i,'..')
  tmp=unlist(strsplit(genes$V9[i],";"))
  idx=which(grepl("Name",tmp))
  allgenes$gene[i]=gsub("Name=","",tmp[idx])
  idx=which(grepl("gene_biotype",tmp))
  allgenes$type[i]=gsub("gene_biotype=","",tmp[idx])
  idx=which(grepl("Dbxref=GeneID:",tmp))
  tmp1=unlist(strsplit(tmp[idx],','))
  idx=which(grepl("Dbxref=GeneID:",tmp1))
  allgenes$geneid[i]=gsub("Dbxref=GeneID:","",tmp1[idx])
}
allgenes=allgenes[allgenes$type=="protein_coding",]
tmp=names(table(allgenes$NC)[1:24])
allgenes=allgenes[allgenes$NC %in% tmp,]
for (i in 1:length(tmp))
{
  idx=which(allgenes$NC==tmp[i])
  tmp1=gsub("NC_","",tmp[i])
  tmp2=unlist(strsplit(tmp1,'.',fixed=T))[1]
  tmp2=as.numeric(tmp2)
  allgenes$chr[idx]=tmp2
}
save(allgenes,file="../result/refgeneshg38codinggenes.RData")

sum(allgenes$geneid %in% genes$gene_id)
load("../result/refgeneshg38codinggenes.RData")

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

## this step is optional!
## here we just simplify the names of the objects, making the code neater
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38

## By default the genome we're using follows the UCSC convention for
## naming chromosome e.g. "chr8".  This step changes that to match our
## SNP data which uses NCBI naming e.g. "8"
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = topsnps$CHR, pos = topsnps$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)
idx = match(topsnps$BP,my_snps$pos)
topsnps$rsid=my_snps$RefSNP_id[idx]
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
topsnps
# CHR                SNP        BP EFF REF  MAF.cases  MAF.ctrls        OR OR_CI_%95_Low OR_CI_%95_High         SE            P       rsid
# 1:   6  chr6:30752534:C:T  30752534   T   C 0.38538682 0.41969887 0.7792844     0.7121179      0.8527861 0.04598689 5.866219e-08 rs28780111
# 2:   6  chr6:31347583:G:A  31347583   A   G 0.29922227 0.24341280 1.3641491     1.2368329      1.5045709 0.04998912 5.232458e-10 rs28771426
# 3:   6  chr6:32616916:C:G  32616916   G   C 0.13917315 0.18005019 0.6873741     0.6085562      0.7764001 0.06213853 1.610256e-09   rs510205
# 4:   4 chr4:157814171:G:A 157814171   A   G 0.12095784 0.08908407 1.4128771     1.2330442      1.6189378 0.06946152 6.497140e-07 rs13146927
# 5:  10 chr10:25753815:C:G  25753815   G   C 0.06344658 0.04705144 1.6111940     1.3347172      1.9449408 0.09605077 6.839426e-07  rs7475986
# 6:  13 chr13:93995797:C:T  93995797   T   C 0.22329104 0.18757842 1.2918834     1.1617107      1.4366422 0.05418849 2.288548e-06  rs9561485
# 7:  19 chr19:30883447:T:C  30883447   T   C 0.14326648 0.18255960 0.7318874     0.6462965      0.8288133 0.06345434 8.701048e-07  rs1466298
# gene closestgene gene_dist
# 1:        FLOT1;IER3        IER3      7986
# 2:       HLA-C;HLA-B       HLA-B      6291
# 3: HLA-DRB1;HLA-DQA1    HLA-DQA1     20489
# 4:    GASK1B;TMEM144      GASK1B    310302
# 5:      ENKUR;GPR158      GPR158    151585
# 6:              GPC6        GPC6         0
# 7:      ZNF536;TSHZ3      ZNF536    169860
topsnps1=allres[allres$P<5e-5,]
positions <- GPos(seqnames = topsnps1$CHR, pos = topsnps1$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)
idx = match(topsnps1$BP,my_snps$pos)
topsnps1$rsid=my_snps$RefSNP_id[idx]
topsnps1$gene=NA
topsnps1$closestgene=NA
topsnps1$gene_dist=NA
gr_topsnps1=GRanges(seqnames = topsnps1$CHR,ranges = IRanges(topsnps1$BP,width = 1))
gr_allgenes=GRanges(seqnames = allgenes$chr,ranges = IRanges(start=allgenes$start,end=allgenes$end))
for (i in 1:nrow(topsnps1))
{
  tmp=distance(gr_allgenes,gr_topsnps1[i])
  idx=which.min(tmp)
  if (tmp[idx]>0)
  {
    if (allgenes$start[idx]>topsnps1$BP[i])
    {
      topsnps1$gene[i]=paste0(allgenes$gene[idx-1],";",allgenes$gene[idx])
    }
    if (allgenes$end[idx]<topsnps1$BP[i])
    {
      topsnps1$gene[i]=paste0(allgenes$gene[idx],";",allgenes$gene[idx+1])
    }
  }else
  {
    topsnps1$gene[i]=allgenes$gene[idx]
  }
  topsnps1$closestgene[i]=allgenes$gene[idx]
  topsnps1$gene_dist[i]=min(tmp,na.rm = T)
}
write.table(topsnps1,file="../result/suggestiveSNPs.csv",sep=",",row.names = F,quote=F)

#adjust for 10pcs
resultfolder = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pc10/"
allfiles = list.files(resultfolder,pattern = "processed.*.txt")
alljobs = unlist(strsplit(allfiles,"__"))
alljobs = alljobs[seq(2,length(alljobs),2)]
alljobs = unlist(strsplit(alljobs,".",fixed=T))
alljobs = as.numeric(alljobs[seq(1,length(alljobs),2)])
if (length(unique(alljobs))!=max(alljobs)) warning("Some results are missing")

allrespc10 = NULL
i=1
for (myfile in paste0(resultfolder,allfiles))
{
  if (i %% 100==0) cat(i,'..')
  tmp = data.table::fread(myfile,header = T)
  allrespc10 = rbind(allrespc10,tmp)
  i = i+1
}

idx = order(allrespc10$CHR,allrespc10$BP)
allrespc10 = allrespc10[idx,]
#save(allrespc10,file="../result/allrespc10.RData")
load("../result/allrespc10.RData")
colnames(allrespc10)[8]
sum(allrespc10$P<5e-8) #18
idx1=which(allres$P<5e-8)
idx2=which(allrespc10$P<5e-8)
idx=unique(c(idx1,idx2))
twores=merge(allres,allrespc10,by=colnames(allres)[1:7])
idx1=match(allres$SNP,twores$SNP)
twores=twores[idx1,]
png(file="../result/twopvalues_scatter.png",res=100)
plot(-log10(twores$P.x),-log10(twores$P.y),cex.axis=1.2,cex.lab=1.2,xlab="-log10(P_pc5)",ylab="-log10(P_pc10")
abline(0,1,col="red")
dev.off()
png("../result/qqplot_pc10.png",res=100)
qqplot(allrespc10$P)
dev.off()
#inflation factor
chisq <- qchisq(1-allrespc10$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda
View(twores[idx,])
sum(twores$P.x<5e-8 & twores$P.y<5e-8) #15
tmp=twores[idx,c(1:5,8,12,13,17)]
tmp$P.x=format(tmp$P.x,digits=3)
tmp$P.y=format(tmp$P.y,digits=3)
tmp$OR.x=format(tmp$OR.x,digits=3)
tmp$OR.y=format(tmp$OR.y,digits=3)
# CHR                 SNP       BP EFF REF  OR.x      P.x  OR.y      P.y
# 1:   6   chr6:31297537:T:C 31297537   C   T 1.286 3.90e-08 1.299 1.43e-08
# 2:   6   chr6:31340108:G:A 31340108   A   G 0.772 2.23e-08 0.769 1.69e-08
# 3:   6   chr6:31344281:C:G 31344281   G   C 0.773 2.50e-08 0.770 1.91e-08
# 4:   6   chr6:31347583:G:A 31347583   A   G 1.364 5.23e-10 1.376 2.16e-10
# 5:   6   chr6:31348083:A:G 31348083   G   A 1.286 4.39e-08 1.298 1.83e-08
# 6:   6   chr6:31351614:A:C 31351614   A   C 1.317 7.72e-09 1.331 2.89e-09
# 7:   6   chr6:31351920:A:G 31351920   G   A 1.296 1.57e-08 1.309 5.95e-09
# 8:   6   chr6:32591896:T:G 32591896   G   T 0.703 3.35e-08 0.706 5.32e-08
# 9:   6   chr6:32616916:C:G 32616916   G   C 0.687 1.61e-09 0.691 2.78e-09
# 10:   6   chr6:32662585:C:G 32662585   C   G 0.766 4.43e-09 0.769 9.55e-09
# 11:   6   chr6:32668956:G:A 32668956   G   A 0.774 1.41e-08 0.777 2.90e-08
# 12:   6   chr6:32686725:A:G 32686725   G   A 0.693 3.30e-09 0.697 5.67e-09
# 13:   6   chr6:32698396:A:G 32698396   G   A 0.694 4.23e-09 0.698 7.27e-09
# 14:   6   chr6:32698518:A:G 32698518   G   A 0.694 4.23e-09 0.698 7.27e-09
# 15:   6   chr6:32700634:G:A 32700634   A   G 0.614 1.83e-08 0.616 2.46e-08
# 16:   6  chr6:32700658:A:AG 32700658  AG   A 0.692 3.37e-09 0.696 5.80e-09
# 17:   6 chr6:32704570:AAT:A 32704570   A AAT 0.719 4.38e-08 0.723 7.28e-08
# 18:   6   chr6:31347894:G:A 31347894   A   G 1.285 5.15e-08 1.296 2.14e-08
# 19:   6   chr6:31351977:G:T 31351977   G   T 1.274 9.53e-08 1.286 4.03e-08
# 20:   6   chr6:32707290:T:C 32707290   T   C 0.768 5.92e-08 0.759 2.68e-08
twores[which(grepl(30752534,twores$SNP)),]
# CHR               SNP       BP EFF REF MAF.cases MAF.ctrls      OR.x OR_CI_%95_Low.x OR_CI_%95_High.x       SE.x          P.x
# 1:   6 chr6:30752534:C:T 30752534   T   C 0.3853868 0.4196989 0.7792844       0.7121179        0.8527861 0.04598689 5.866219e-08
# OR.y OR_CI_%95_Low.y OR_CI_%95_High.y       SE.y          P.y
# 1: 0.7796803       0.7117332         0.854114 0.04652168 8.815928e-08

cor(twores$OR.x,twores$OR.y) #0.998
METALdat1=read.table(METALfile,head=T)
mapdat=read.table("../result/locusdata_snpmap.txt")
all(mapdat$V1==METALdat1$MarkerName) #T
idx=match(mapdat$V2,allrespc10$SNP)
METALdat1$Pvalue=allrespc10$P[idx]
METALfile1="../result/sig_snps_cp10_metalfile.txt"
write.table(METALdat1,file=METALfile1,col.names=T,row.names = F,sep="\t",quote=F)
title=""
cmd=paste0("locuszoom"," --metal ",METALfile1," --flank 2000kb"," --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs28771426 --pvalcol Pvalue --plotonly --prefix zoomlocuspc10 --ld-vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz ","ymax=",10," signifLine=",-log10(5*1e-8)," signifLineColor=red ", " axisTextSize=1.4 legendSize=1.0  axisSize=1.4 xlabPos=-2.9")
cmd

library(qqman)
png("../result/manhattan_pc10.png",res=100, width=1200)
manhattan(allrespc10, chr="CHR", bp="BP", snp="SNP", p="P", suggestiveline = -log10(5e-5),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

idx=which(allrespc10$P<1e-7)
plot(allrespc10$BP[idx],-log10(allrespc10$P[idx]))
allrespc10$SNP[idx]

#work on genotyped ChrX results
resultfolder = "/data/DCEGLeiSongData/BinBB/zhiwei/GWAS_and_HBV_DNA_Load/scripts/gwas/out/"
allfiles = list.files(resultfolder,pattern = "REVEAL.*.txt")
alljobs = unlist(strsplit(allfiles,"__"))
alljobs = alljobs[seq(2,length(alljobs),2)]
alljobs = unlist(strsplit(alljobs,".",fixed=T))
alljobs = as.numeric(alljobs[seq(1,length(alljobs),2)])
if (length(unique(alljobs))!=max(alljobs)) warning("Some results are missing")

allres = NULL
i=1
for (myfile in paste0(resultfolder,allfiles))
{
  if (i %% 100==0) cat(i,'..')
  tmp = fread(myfile,header = T)
  allres = rbind(allres,tmp)
  i = i+1
}

idx = order(allres$CHR,allres$BP)
allres = allres[idx,]
chrXres=allres[allres$CHR==23,]
pheno = read.table("../data/ana_final_VL_new.csv",sep=",",header=T)
tmp=fread("/data/DCEGLeiSongData/BinBB/zhiwei/GWAS_and_HBV_DNA_Load/scripts/gwas/splited/REVEAL_TLCN_All.QCed.traw__55.gz")
idx=which(tmp$CHR==23)
tmp=as.data.frame(tmp[idx,])
tmp1=rep(NA,3240)
for (i in 1:3240)
{
  tmp2=unlist(strsplit(colnames(tmp)[i+6],"_"))
  tmp1[i]=paste0(tmp2[1],"_",tmp2[2])
}
colnames(tmp)[7:ncol(tmp)]=tmp1
malesamples=pheno$chip_id[pheno$SEX==1]
femalesamples=pheno$chip_id[pheno$SEX==0]
idx1=match(malesamples,colnames(tmp))
maledat=cbind(tmp[,1:6],tmp[,idx1])
idx1=match(femalesamples,colnames(tmp))
femaledat=cbind(tmp[,1:6],tmp[,idx1])


