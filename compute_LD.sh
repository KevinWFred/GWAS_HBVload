#!/usr/bin/env bash
module load plink/1.9

infolder=/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/
outfolder="$infolder"processed/

plink --bfile "$outfolder"chr6_filter --extract /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata_snp.txt --make-bed --out /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata
plink --bfile /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata --update-name /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata_snpmap.txt 1 2 --make-bed --out /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata
plink --bfile /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata --update-name /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata_snpmap.txt 1 2 --recode vcf --out /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata

#change filter to PASS
head -7 /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata.vcf > /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf
cat /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata_pass.vcf >> /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf

ml samtools
bgzip -c /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf >/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz
tabix -p vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz

#Error: while calculating LD from VCF file: index SNP position 31315360 does not exist in file.
locuszoom --metal ../result/sig_snps_metalfile.txt --flank 2000kb --build hg19 --pop ASN --source 1000G_March2012 --refsnp rs28771426 --pvalcol Pvalue --plotonly --prefix zoomlocus --ld-vcf /data/DCEGLeiSongData/Kevin/HBVloadGwas/result/locusdata1.vcf.gz ymax=10 signifLine=7.30102999566398 signifLineColor=red  axisTextSize=1.4 legendSize=1.2  axisSize=1.4 xlabPos=-2.9