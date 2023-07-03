#!/usr/bin/env Rscript

library(data.table)
logistres=as.data.frame(fread("../result/HBsAg_seroclearance_sumstat.csv",sep=","))

#Cox res
resultfolder = "../result/cox_sero/"
allfiles = list.files(resultfolder,pattern = "processed.*.txt")
alljobs = unlist(strsplit(allfiles,"__"))
alljobs = alljobs[seq(2,length(alljobs),2)]
alljobs = unlist(strsplit(alljobs,".",fixed=T))
alljobs = as.numeric(alljobs[seq(1,length(alljobs),2)])
if (length(unique(alljobs))!=max(alljobs)) warning("Some results are missing")

coxres = NULL
i=1
for (myfile in paste0(resultfolder,allfiles))
{
  if (i %% 100==0) cat(i,'..')
  tmp = fread(myfile,header = T)
  coxres = rbind(coxres,tmp)
  i = i+1
}

idx = order(coxres$CHR,coxres$BP)
coxres = coxres[idx,]
coxres$rsid=logistres$rsid
all(logistres$ID==coxres$SNP)
table(logistres$a1==coxres$EFF) #T
colnames(coxres)[2]="ID"
colnames(coxres)[7]="beta"
colnames(coxres)[4:5]=c("a1","a2")
idx=which(logistres$a2==coxres$a1)
coxres$beta[idx]=-coxres$beta[idx]
tmp=coxres$a1[idx]
coxres$a1[idx]=coxres$a2[idx]
coxres$a2[idx]=tmp
write.csv(coxres[,c(2,9,1,3:8)],"../result/HBsAg_seroclearance_COX_sumstat.csv")

png(filename="../result/scater_logis_cox_sero_pvalue.png",res=100)
ylim=range(c(-log10(logistres$P),-log10(coxres$P)))
plot(-log10(logistres$P),-log10(coxres$P),xlab="-log10(Logistic p-value)",ylab="-log10(COX pvalue)",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()
png(filename="../result/scater_logis_cox_sero_beta.png",res=100)
ylim=range(c(log(logistres$OR),coxres$beta))
plot(log(logistres$OR),coxres$beta,xlab="Logistic beta",ylab="COX beta",cex.lab=1.2,cex.axis=1.2,xlim=ylim,ylim=ylim)
abline(0,1,lty=2,col="red")
dev.off()

idx1=which(log(logistres$OR)< -0.1 & coxres$beta>0.1)
plot(log(logistres$OR[idx1]),coxres$beta[idx1],xlab="Logistic beta",ylab="COX beta",cex.lab=1.2,cex.axis=1.2)
