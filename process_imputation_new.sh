#!/usr/bin/env bash
#code used to process imputed genotype for HBVload GWAS project
#it QCs the genotype (HWE P>=1e-6, MAF>0.01, imputation score >0.3, and missing rate <0.01) and generates the genotype (dosage) data.

module load plink/1.9

#https://si.biostat.washington.edu/sites/default/files/modules/RecommendedReading_Session6.pdf , HRC doesn't include indels

#folder contains imputed data
infolder=/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/
outfolder="$infolder"processed_new/
mkdir -p $outfolder

cd $outfolder
#combine all the reads
readallchr(){
  local chr="$1"
  echo $chr
  prefix="$outfolder"chr${chr}
  #vcf file
  vcffile="$infolder"chr${chr}.dose.vcf.gz
  plink --vcf $vcffile  --make-bed --out $prefix --memory 6400 --threads 1
}
for chr in {1..22}
do
  readallchr $chr &
done

#merge chr1-chr22
rm mergelist.txt
for chr in {1..22}
do
  echo chr${chr} >> mergelist.txt
done

plink --merge-list mergelist.txt --make-bed --out all --memory 128000 --threads 12

#QC
plink --bfile all --maf 0.01 --hwe 0.000001 --geno 0.01 --keep /data/DCEGLeiSongData/Kevin/HBVloadGwas/data/samples.txt --make-bed --out all1 --memory 128000 --threads 12
plink --bfile all --maf 0.05 --hwe 0.000001 --geno 0.01 --keep /data/DCEGLeiSongData/Kevin/HBVloadGwas/data/samples.txt --make-bed --out all1maf05 --memory 128000 --threads 12

zcat ${infolder}chr22.info.gz|head -n 1 > allinfo.txt
for chr in {1..22}
do
  zcat ${infolder}chr${chr}.info.gz |tail -n +2 -q >> allinfo.txt
done

plink --bfile all1 --qual-scores allinfo1.txt 7 1 1 --qual-threshold 0.3 --make-bed --out processed --memory 128000 --threads 12
plink --bfile processed --recode A-transpose --out processed --memory 128000 --threads 12
gzip processed.traw
plink --bfile processed --freq --out processed --memory 100000 --threads 12

wc -l allinfo1.txt #292174913
wc -l all.bim #292174934
# 0 variants removed due to missing genotype data (--geno).
# --hwe: 1808 variants removed due to Hardy-Weinberg exact test.
# 284520511 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 7652615 variants and 3240 people pass filters and QC.
wc -l all1.bim #7652615
wc -l processed.bim #7608863


