#/usr/bin/env Rscript
library(data.table)
plink="/usr/local/apps/plink/1.9/plink"
plink2="/usr/local/apps/plink/2.3-alpha/plink2"

pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
load(pheno.Rdata.file)
hbeagdat=read.csv("../data/GLMdata.csv")
idx=match(rownames(pheno),hbeagdat$id_new)
all(rownames(pheno)==hbeagdat$id_new[idx])
pheno$hbeag=hbeagdat$hbeag[idx]

#form phenotype data
tmp=data.frame(FID=0,IID=rownames(pheno),hbeag=pheno$hbeag)
tmp$hbeag=tmp$hbeag+1
tmp$hbeag[which(is.na(tmp$hbeag))]=-9
# table(tmp$hbeag,useNA="ifany")
# 1    2 
# 2743  497 
write.table(tmp,file="../result/hbeag_pheno.txt",row.names = F,col.names = T,sep=" ",quote=F)


#Run GWAS
prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
cmd=paste0(plink2," --bfile ",prefix,"  --pheno ../result/hbeag_pheno.txt --covar ../result/sero_cov.txt  --logistic hide-covar --ci 0.95 --out ../result/hbeag")
system(cmd)

cmd=paste0(plink2," --bfile ",prefix," --pheno ../result/hbeag_pheno.txt --covar  ../result/sero_cov.txt --logistic no-firth hide-covar --ci 0.95 --out ../result/hbeag")
system(cmd)
logistres=as.data.frame(fread("../result/hbeag.hbeag.glm.logistic"))
logistres=logistres[logistres$TEST=="ADD",]


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
png("../result/hbeag_qqplot.png",res=100)
qqplot(logistres$P) #1.05
dev.off()
library(qqman)
colnames(logistres)[1]="CHR"
png("../result/hbeag_manhattan.png",res=100, width=1200)
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
tmp=logistres[,c(3,16,1,2,17,18,9,10,14)]
write.csv(tmp,file="../result/hbeag_sumstat.csv",row.names = F,quote=F)
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
tmp=data.frame(topsnps[,c(3,16,1,2,17,18,9,10,14,19,20,21)])
tmp=tmp[order(tmp$P),]
write.csv(tmp,file="../result/hbeag_topsnps.csv",row.names = F)
tmp1=as.data.frame(fread("../result/hbeag_topsnps.csv"))
