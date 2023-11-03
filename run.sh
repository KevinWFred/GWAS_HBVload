#!/usr/bin/env bash
#to split genomewied genotype into multiple jobs (each job has 10000 SNPs) and submit

basedir=/data/DCEGLeiSongData/Kevin/HBVloadGwas/
cd $basedir

#split
if [ -d splited ]
then 
  rm -r splited
fi

if [ ! -d splited ]
then
  mkdir -p splited
  cd splited
  cat ${basedir}TOPMed_Imputation/processed/processed.traw | parallel --header : --pipe -N10000 'cat | gzip > processed.traw__{#}.gz'
fi

#create swarm file
if [ -f ${basedir}code/run.swarm ]
then
  rm  ${basedir}code/run.swarm
fi

rscript=" ${basedir}code/ordinal.lr.R"

splited_file_dir="${basedir}splited"
out_dir="${basedir}result"

ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
  ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
  printf "Rscript $rscript $ifile ${out_dir}/$ofile\n" >> ${basedir}code/run.swarm
  echo $ifile
  echo $ofile
done

#submit job
if [ -d ${basedir}logs ]
then
	rm -r ${basedir}logs
fi

mkdir -p ${basedir}logs
cd ${basedir}logs

swarm -f ${basedir}code/run.swarm -g 16 --module R/4.2.0 --time=10:00:00

#For HBeAg negative samples
rscript="/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/ordinal.lr.R"
splited_file_dir="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited"
out_dir="/data/BB_Bioinformatics/Kevin/HBV_GWAS/result/ordinal_hbeag_neg"

ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
  ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
  printf "Rscript $rscript $ifile ${out_dir}/$ofile negative \n" >> run_hbeag_neg.swarm
  echo $ifile
  echo $ofile
done
swarm -f /data/BB_Bioinformatics/Kevin/HBV_GWAS/code/run_hbeag_neg.swarm -g 16 --module R/4.2.0 --time=12:00:00

#for MAF of 0.01
#split
if [ -d splited_maf01 ]
then 
  rm -r splited_maf01
fi

if [ ! -d splited_maf01 ]
then
  mkdir -p splited_maf01
  cd splited_maf01
  zcat ${basedir}TOPMed_Imputation/processed_new/processed.traw.gz | parallel --header : --pipe -N10000 'cat | gzip > processed.traw__{#}.gz'
fi

splited_file_dir="${basedir}splited_maf01"
out_dir="${basedir}result_maf01"
if [ ! -d $out_dir ]
then
  mkdir -p $out_dir
if

rscript="/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/ordinal.lr.R"
rm ${basedir}code/run_maf01.swarm
ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
  ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
  printf "Rscript $rscript $ifile ${out_dir}/$ofile\n" >> ${basedir}code/run_maf01.swarm
  echo $ifile
  echo $ofile
done

#submit job
if [ -d ${basedir}logs_maf01 ]
then
	rm -r ${basedir}logs_maf01
fi

mkdir -p ${basedir}logs_maf01
cd ${basedir}logs_maf01
#8914535,8999929
swarm -f ${basedir}code/run_maf01.swarm -g 16 --module R/4.3.0 --time=10:00:00


#For HBeAg negative samples
rscript="/data/BB_Bioinformatics/Kevin/HBV_GWAS/code/ordinal.lr.R"
splited_file_dir="/data/DCEGLeiSongData/Kevin/HBVloadGwas/splited_maf01"
out_dir="/data/BB_Bioinformatics/Kevin/HBV_GWAS/result/ordinal_hbeag_neg_maf01"

ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
  ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
  printf "Rscript $rscript $ifile ${out_dir}/$ofile negative \n" >> run_hbeag_neg_maf01.swarm
  echo $ifile
  echo $ofile
done
mv run_hbeag_neg_maf01.swarm /data/BB_Bioinformatics/Kevin/HBV_GWAS/code/run_hbeag_neg_maf01.swarm
cd ${basedir}logs_maf01
#9065701
swarm -f /data/BB_Bioinformatics/Kevin/HBV_GWAS/code/run_hbeag_neg_maf01.swarm -g 16 --module R/4.2.0 --time=12:00:00
