#!/usr/bin/env python
#https://github.com/JonJala/mtag/wiki/Tutorial-1:-The-Basics
#python 2.7: source /data/BB_Bioinformatics/Kevin/tools/python2/bin/activate
import pandas as pd
import os
import numpy as np
import subprocess
import sys

plink='/usr/local/apps/plink/1.9.0-beta4.4/plink'
plink2='/usr/local/apps/plink/2.3-alpha/plink2'
os.chdir("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")
hbeagsum=pd.read_csv("../result/hbeag.hbeag.glm.logistic",sep="\t")
#snpid    chr    bpos    a1    a2    freq    z    pval    n
hbeagsum1=pd.DataFrame({"snpid":hbeagsum["ID"],"chr":hbeagsum["#CHROM"],"bpos":hbeagsum["POS"],"a1":hbeagsum["ALT"],"a2":hbeagsum["REF"],
                       "freq":0,"z":hbeagsum["Z_STAT"],"pval":hbeagsum["P"],"n":2743+497})
idx=np.where(hbeagsum["ALT"]!=hbeagsum["A1"])[0]
hbeagsum1.loc[idx,"a1"]=hbeagsum.loc[idx,"REF"]
hbeagsum1.loc[idx,"a2"]=hbeagsum.loc[idx,"ALT"]

tmp=pd.DataFrame({"snp":hbeagsum["ID"],"a1":hbeagsum["A1"]})
tmp.to_csv("../result/hbeag.a1",index=False,sep="\t",header=False)

prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
cmd=[plink2 + " --bfile " + prefix + " --alt1-allele ../result/../result/hbeag.a1 --make-bed --out ../result/hbeag "]
subprocess.call(cmd,shell=True)
cmd=[plink2 + " --bfile ../result/hbeag --keep-allele-order --freq --out ../result/hbeag "]
subprocess.call(cmd,shell=True)
freq=pd.read_csv("../result/hbeag.afreq",sep="\t")
np.all(freq["ID"]==hbeagsum1["snpid"])
hbeagsum1["freq"]=freq["ALT_FREQS"]
rsiddat=pd.read_csv("../result/hbeag_sumstat.csv")
all(rsiddat["ID"]==hbeagsum1["snpid"])
np.sum(rsiddat["ID"].isna())
hbeagsum1["snpid"]=rsiddat["rsid"]
idx=(hbeagsum1["a2"].isin(["A","T","G","C"])) & (hbeagsum1["a1"].isin(["A","T","G","C"]))
idx.value_counts()
hbeagsum1.loc[idx,].to_csv("../result/hbeagsum_mtag.txt",index=False,sep="\t")

#HBVload
hbvsum=pd.read_csv("../result/ordinal_lr_result.txt",sep="\t")
hbvsum1=pd.DataFrame({"snpid":hbvsum["SNP"],"chr":hbvsum["CHR"],"bpos":hbvsum["BP"],"a1":hbvsum["EFF"],"a2":hbvsum["REF"],
                       "freq":0,"z":np.log(hbvsum["OR"])/hbvsum["SE"],"pval":hbvsum["P"],"n":3240})

tmp=pd.DataFrame({"snp":hbvsum["SNP"],"a1":hbvsum["EFF"]})
tmp.to_csv("../result/processed_ordinal.a1",index=False,sep="\t",header=False)
cmd=[plink2 + " --bfile " + prefix + " --a1-allele ../result/processed_ordinal.a1 --make-bed --out ../result/processed_ordinal "]
subprocess.call(cmd,shell=True)
cmd=[plink2 + " --bfile ../result/processed_ordinal --keep-allele-order --freq --out ../result/processed_ordinal "]
subprocess.call(cmd,shell=True)
freq=pd.read_csv("../result/processed_ordinal.afreq",sep="\t")
np.all(freq["ID"]==hbvsum1["snpid"])
hbvsum1["freq"]=freq["ALT_FREQS"]
all(rsiddat["ID"]==hbvsum1["snpid"])
hbvsum1["snpid"]=rsiddat["rsid"]
idx=(hbvsum1["a2"].isin(["A","T","G","C"])) & (hbvsum1["a1"].isin(["A","T","G","C"]))
idx.value_counts()
hbvsum1.loc[idx,].to_csv("../result/hbvsum_mtag.txt",index=False,sep="\t")

cmd=["/data/BB_Bioinformatics/Kevin/tools/mtag/mtag.py --ld_ref_panel /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --sumstats ../result/hbeagsum_mtag.txt,../result/hbvsum_mtag.txt  --force --out ../result/mtag_hbvload_hbeag --n_min 0.0  --stream_stdout"]
subprocess.call(cmd,shell=True)
