#!/usr/bin/env python
#https://github.com/JonJala/mtag/wiki/Tutorial-1:-The-Basics

import pandas as pd
import os
import numpy as np
import subprocess
import sys

plink='/usr/local/apps/plink/1.9/plink'
plink2='/usr/local/apps/plink/2.3-alpha/plink2'
os.chdir("/data/BB_Bioinformatics/Kevin/HBV_GWAS/code")
serosum=pd.read_csv("../result/seroclearance.seroclearance.glm.logistic",sep="\t")
#snpid    chr    bpos    a1    a2    freq    z    pval    n
serosum1=pd.DataFrame({"snpid":serosum["ID"],"chr":serosum["#CHROM"],"bpos":serosum["POS"],"a1":serosum["ALT"],"a2":serosum["REF"],
                       "freq":0,"z":serosum["Z_STAT"],"pval":serosum["P"],"n":2416+535})
idx=np.where(serosum["ALT"]!=serosum["A1"])[0]
serosum1.loc[idx,"a1"]=serosum.loc[idx,"REF"]
serosum1.loc[idx,"a2"]=serosum.loc[idx,"ALT"]

tmp=pd.DataFrame({"snp":serosum["ID"],"a1":serosum["A1"]})
tmp.to_csv("../result/sero.a1",index=False,sep="\t",header=False)

prefix="/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/processed"
cmd=[plink2 + " --bfile " + prefix + " --keep ../result/sero_samples.txt --a1-allele ../result/../result/sero.a1 --make-bed --out ../result/sero "]
subprocess.call(cmd,shell=True)
cmd=[plink2 + " --bfile ../result/sero --keep-allele-order --freq --out ../result/sero "]
subprocess.call(cmd,shell=True)
freq=pd.read_csv("../result/sero.afreq",sep="\t")
np.all(freq["ID"]==serosum1["snpid"])
serosum1["freq"]=freq["ALT_FREQS"]
rsiddat=pd.read_csv("../result/HBsAg_seroclearance_sumstat.csv")
all(rsiddat["ID"]==serosum1["snpid"])
np.sum(rsiddat["ID"].isna())
serosum1["snpid"]=rsiddat["rsid"]
idx=(serosum1["a2"].isin(["A","T","G","C"])) & (serosum1["a1"].isin(["A","T","G","C"]))
idx.value_counts()
serosum1.loc[idx,].to_csv("../result/serosum_mtag.txt",index=False,sep="\t")

#qHBsAg
qHBsAg=pd.read_csv("../result/qHBsAg_sumstat.csv",sep=",")
qHBsAg1=pd.DataFrame({"snpid":qHBsAg["ID"],"chr":qHBsAg["CHR"],"bpos":qHBsAg["POS"],"a1":qHBsAg["a1"],"a2":qHBsAg["a2"],
                       "freq":0,"z":qHBsAg["BETA"]/qHBsAg["SE"],"pval":qHBsAg["P"],"n":3186})

tmp=pd.DataFrame({"snp":qHBsAg["ID"],"a1":qHBsAg["a1"]})
tmp.to_csv("../result/processed_qHBsAg.a1",index=False,sep="\t",header=False)
cmd=[plink2 + " --bfile " + prefix + " --a1-allele ../result/processed_qHBsAg.a1 --make-bed --out ../result/processed_qHBsAg "]
subprocess.call(cmd,shell=True)
cmd=[plink2 + " --bfile ../result/processed_qHBsAg --keep-allele-order --freq --out ../result/processed_qHBsAg "]
subprocess.call(cmd,shell=True)
freq=pd.read_csv("../result/processed_qHBsAg.afreq",sep="\t")
np.all(freq["ID"]==qHBsAg1["snpid"])
np.all(freq["ALT"]==qHBsAg1["a1"])
qHBsAg1["freq"]=freq["ALT_FREQS"]
np.all(rsiddat["ID"]==qHBsAg1["snpid"])
qHBsAg1["snpid"]=rsiddat["rsid"]
idx=(qHBsAg1["a2"].isin(["A","T","G","C"])) & (qHBsAg1["a1"].isin(["A","T","G","C"]))
idx.value_counts()
qHBsAg1.loc[idx,].to_csv("../result/qHBsAg_mtag.txt",index=False,sep="\t")

cmd=["/data/BB_Bioinformatics/Kevin/tools/mtag/mtag.py --sumstats ../result/serosum_mtag.txt,../result/qHBsAg_mtag.txt  --force --out ../result/mtag_qHBsAg_sero --n_min 0.0  --stream_stdout"]
subprocess.call(cmd,shell=True)
