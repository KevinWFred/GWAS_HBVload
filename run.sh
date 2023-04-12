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
