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

cmd=["/data/BB_Bioinformatics/Kevin/tools/mtag/mtag.py --ld_ref_panel /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --sumstats ../result/hbvsum_mtag.txt,../result/hbeagsum_mtag.txt,../result/serosum_mtag.txt  --force --out ../result/mtag_hbvload_hbeag_sero --n_min 0.0  --stream_stdout"]
subprocess.call(cmd,shell=True)
cmd=["/data/BB_Bioinformatics/Kevin/tools/mtag/mtag.py --ld_ref_panel /data/BB_Bioinformatics/Kevin/tools/LDSC/1000G_EAS_Phase3/ --sumstats ../result/hbvsum_mtag.txt,../result/hbeagsum_mtag.txt,../result/qHBsAg_mtag.txt  --force --out ../result/mtag_hbvload_hbeag_qHBsAg --n_min 0.0  --stream_stdout"]
subprocess.call(cmd,shell=True)
