#!/usr/bin/env Rscript

library("readxl")

#read HLA data
parse.hla.raw.data <- function(raw.data.file = "../data/ana_final_VL_new.csv"){
  # browser()
  raw.data = read.csv(file = raw.data.file, header = TRUE, check.names = FALSE, stringsAsFactors = F)
  samples = raw.data$subject_id
  
  idx0 = grepl("_allele", colnames(raw.data))
  fakesnps = matrix(colnames(raw.data)[idx0], 2, )
  
  fakesnps.data.matrix.list = list()
  
  for(k in 1:ncol(fakesnps)){
    cur.fakesnp.name = unlist(strsplit(fakesnps[1, k], "_"))[1]
    
    cur.fakesnp.raw.data = as.matrix(raw.data[, fakesnps[, k]])
    sub.fakesnps = unique(as.vector(cur.fakesnp.raw.data))
    
    x = matrix(0, length(samples), length(sub.fakesnps))
    colnames(x) = sub.fakesnps
    rownames(x) = samples
    
    for(j in 1:ncol(x)){
      for(i in 1:ncol(cur.fakesnp.raw.data)){
        idx = cur.fakesnp.raw.data[, i] == colnames(x)[j]
        x[idx, j] = x[idx, j] + 1
      }
    }
    
    fakesnps.data.matrix.list[[cur.fakesnp.name]] = x
    cat(k, "\t", cur.fakesnp.name, "\n")
  }
  # browser()

  save(fakesnps.data.matrix.list, file = "../result/fakesnps.data.matrix.list.Rdata")
  pheno=raw.data[,!idx0]
  for (i in 1:length(fakesnps.data.matrix.list))
  {
    pheno=cbind(pheno,fakesnps.data.matrix.list[[i]])
  }

  save(pheno, file = "../result/pheno.Rdata")
}

selsnps=c("rs28780111","rs28771426","rs510205","rs13146927","rs7475986", "rs9561485", "rs1466298","rs4410197")
sumstat=fread("../result/sumstat.cvs")
idx=match(selsnps,sumstat$rsid)
selsnpid=data.frame(snp=sumstat$SNP[idx])

write.table(selsnpid,file="../result/selectedsnps.txt",row.names = F,quote=F)
plink="/usr/local/apps/plink/1.9/plink"

bim=fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed.bim")
idx=match(selsnpid$snp,bim$V2)
cmd=paste0(plink," --bfile /data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed --extract ../result/selectedsnps.txt --recode A-transpose --out ../result/selectedsnps")
system(cmd)
cmd=paste0(plink," --bfile /data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed --extract ../result/selectedsnps.txt --recodeAD --out ../result/selectedsnps")
system(cmd)

# tmp=as.data.frame(data.table::fread("../result/selectedsnps.traw"))
# load("../result/pheno.Rdata")
# tmp1=t(tmp[,7:ncol(tmp)])
# colnames(tmp1)=selsnps[match(tmp$SNP,selsnpid[,1])]
# rownames(tmp1)=gsub("0_","",rownames(tmp1))
# idx=match(pheno$subject_id,rownames(tmp1))
# tmp1=tmp1[idx,]
# pheno0=cbind.data.frame(pheno,tmp1)


load("../result/pheno.Rdata")
tmp=as.data.frame(data.table::fread("../result/selectedsnps.raw"))
idx=match(pheno$subject_id,tmp$IID)
tmp=tmp[idx,]
tmp=tmp[,7:ncol(tmp)]
for (i in 1:length(selsnps))
{
  colnames(tmp)=gsub(selsnpid[i,1],selsnps[i],colnames(tmp))
}

pheno0=cbind.data.frame(pheno,tmp)
#tmp1=read.table("../result/extract_data.txt",header=T)
#all(pheno0[,1:198]==tmp1[,1:198])
write.table(pheno0,file="../result/extract_data.txt",row.names = F,quote=F,sep="\t")

#all snp data
load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/info_res.RData")
genome <- BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(genome) <- "NCBI"
positions <- GPos(seqnames = allres$CHR, pos = allres$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
my_snps = as.data.frame(my_snps)
tmp1=paste0(allres$CHR,"_",allres$BP)
tmp2=paste0(my_snps$seqnames,"_",my_snps$pos)
idx = match(tmp1,tmp2)
allres$rsid=my_snps$RefSNP_id[idx]

#check imputation
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
  tmp = tmp[,1:8]
  tmp$Genotyped[which(tmp$Genotyped=="Imputed")]=0
  tmp$Genotyped[which(tmp$Genotyped=="Genotyped")]=1
  allinfo = rbind(allinfo,tmp)
}
save(allinfo,file="../result/info.RData")
idx=match(allres$SNP,allinfo$SNP)
allres$Genotyped=allinfo$Genotyped[idx]
allres$Beta=log(allres$OR)
allres1=allres[,c("SNP","rsid","CHR","BP","EFF","REF","Beta","SE","P","Genotyped")]
allres1$Genotyped=1-as.numeric(allres1$Genotyped)
fwrite(allres1,file="../result/sumstat.cvs",row.names = F) # HBV_GWAS_sumstat.csv in BOX

#work on GWAS with HBeAg negative
resultfolder = "../result/ordinal_hbeag_neg/"
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

genome <- BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(genome) <- "NCBI"
positions <- GPos(seqnames = allres$CHR, pos = allres$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)
my_snps = as.data.frame(my_snps)
tmp1=paste0(allres$CHR,"_",allres$BP)
tmp2=paste0(my_snps$seqnames,"_",my_snps$pos)
idx = match(tmp1,tmp2)
allres$rsid=my_snps$RefSNP_id[idx]
allres$Beta=log(allres$OR)
allres1=allres[,c("SNP","rsid","CHR","BP","EFF","REF","Beta","SE","P")]
fwrite(allres1,file="../result/HBV_HBeAg_negativeGWAS_sumstat.csv",row.names = F) 
tmp=as.data.frame(fread("../result/sumstat.cvs",sep=","))
all(tmp$SNP==allres1$SNP)
