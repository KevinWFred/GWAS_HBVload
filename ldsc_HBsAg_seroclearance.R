#!/usr/bin/env Rscript


.libPaths(c("/data/wangx53",.libPaths()))
library(data.table)

HBsAg_serosumdat=as.data.frame(fread("../result/seroclearance.seroclearance.glm.logistic"))
HBsAg_serosumdat=as.data.frame(fread("../result/HBsAg_seroclearance_sumstat.csv"))
HBsAg_serosumdat1=data.frame(snpid=HBsAg_serosumdat$rsid,A1=HBsAg_serosumdat$a1,A2=HBsAg_serosumdat$a2,se=HBsAg_serosumdat$`LOG(OR)_SE`,Z=log(HBsAg_serosumdat$OR)/HBsAg_serosumdat$`LOG(OR)_SE`,MAF=NA,P=HBsAg_serosumdat$P,N=2416+535,info=1)
freq=as.data.frame(fread("../result/sero.afreq",sep="\t"))
idx=match(HBsAg_serosumdat$ID,freq$ID)
HBsAg_serosumdat1$MAF=freq$ALT_FREQS[idx]
write.table(HBsAg_serosumdat1,file="../result/HBsAg_seroclearance_ldsc_sumdat.txt",sep="\t",row.names = F,quote=F)
cmd="ml ldsc; /usr/local/apps/ldsc/1.0.1-20200724/bin/munge_sumstats.py --sumstats ../result/HBsAg_seroclearance_ldsc_sumdat.txt --N 2951 --merge-alleles /data/BB_Bioinformatics/simulated_multi_ances_genotype_600K/snp_infor/hm3rsid.txt --signed-sumstats Z,0 --info-min 0.3  --chunksize 500000 --out ../result/HBsAg_seroclearance_ldsc_sumdat_align"
system(cmd)
cmd="ml ldsc; /usr/local/apps/ldsc/1.0.1-20200724/bin/ldsc.py --h2 ../result/HBsAg_seroclearance_ldsc_sumdat_align.sumstats.gz --ref-ld-chr /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --w-ld-chr /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --out ../result/HBsAg_seroclearance_ldsc --pop-prev 0.01 --samp-prev 0.181 "
system(cmd) #-0.0763 (0.119) # if the true heritability is low, that value plus some sampling error can lead to a negative estimate (especially if the sample size is small and thus there’s lots of sampling variation). 
cmd="ml ldsc; /usr/local/apps/ldsc/1.0.1-20200724/bin/ldsc.py --h2 ../result/HBsAg_seroclearance_ldsc_sumdat_align.sumstats.gz --ref-ld-chr /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --w-ld-chr /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --out ../result/HBsAg_seroclearance_ldsc_obs "
system(cmd)

