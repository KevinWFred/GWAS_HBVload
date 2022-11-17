#!/usr/bin/env bash
module load plink/1.9
#to use snpflip
#ml Python/3.7.4-GCCcore-8.3.0
#install snpflip
#python -m venv /fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1
#source /fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1/bin/activate
#python -m pip install snpflip

#source /fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1/bin/activate
#snpflip=/fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1/bin/snpflip
#reference=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/reference/human_g1k_v37.fasta
#hg38=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/hg38/hg38.fa

#https://si.biostat.washington.edu/sites/default/files/modules/RecommendedReading_Session6.pdf , HRC doesn't include indels

infolder=/data/DCEGLeiSongData/Kevin/TOPMed_Imputation/
outfolder="$infolder"processed/
mkdir -p $outfolder
#run the following for each dataset
#filter

#fileter2 maf=0.01,allvars, they are the same as above. no ins/del included in HRC
# indataset="cambridgewtccc_qc_hrc"
# outdataset="cambridgewtccc_qc_hrc_maf001_var"
# 
# indataset="beacondbgapcontrol_qc_hrc"
# outdataset="beacondbgapcontrol_qc_hrc_maf001_var"

# #check missing
# checkmissing(){
#   local chr="$1"
#   echo $chr
#   vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}.dose.vcf.gz
#   plink --vcf $vcffile  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 2560
# }


cd $outfolder
#filtering

filtering(){
  local chr="$1"
  echo $chr
  prefix="$outfolder"chr${chr}_filter
  vcffile="$infolder"chr${chr}.dose.vcf.gz
  infofile="$infolder"chr${chr}.info.gz
  infofile1="$infolder"chr${chr}.info
  #if [[ ! -f $infofile1 ]];then
  #   gunzip -c $infofile > $infofile1
  #fi
  #plink --vcf $vcffile  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 6400
  #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp.multi-allelic
  #plink --bfile tmp_s1_$chr  --bmerge tmp_s1_$chr --merge-mode 6 --out tmp_$chr   --memory 6400
  #if [[ -f tmp_$chr.missnp ]];then
  #  plink --bfile tmp_s1_$chr  --exclude tmp_$chr.missnp --make-bed --out tmp_s1_$chr --memory 6400
  #fi
  #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. (In fact, some tools, such as BEAGLE 4, require such duplicate variants to be merged/removed during preprocessing.) --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
  #plink --bfile tmp_s1_$chr  --list-duplicate-vars --out tmp_$chr --memory 6400
  #plink --bfile tmp_s1_$chr  --exclude tmp_$chr.dupvar --make-bed --out tmp_s2_$chr --memory 6400
  #plink --bfile tmp_s2_$chr  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --qual-scores $infofile1 7 1 1 --qual-threshold 0.3 --make-bed --out $prefix --memory 6400
  plink --bfile $prefix --make-bed --keep /data/DCEGLeiSongData/Kevin/data/samples.txt --out "$prefix" --memory 6400
  #plink --bfile $prefix --recode A-transpose --out "$prefix" --memory 6400

  #rm $infofile1
  #rm tmp_s1_$chr.*
  #rm tmp_$chr.*
  #rm tmp_s2_$chr.*
  
}
for chr in {1..22}
do
  filtering $chr &
done

rm mergelist.txt
for chr in {1..22}
do
echo chr${chr}_filter>> mergelist.txt
done
plink --merge-list mergelist.txt --make-bed --out processed --memory 6400
plink --bfile processed --recode A-transpose --out processed --memory 6400

cat chr22_filter.traw | parallel --header : --pipe -N10000 'cat | gzip > chr22_filter.traw__{#}.gz'