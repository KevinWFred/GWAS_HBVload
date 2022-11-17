basedir=/data/DCEGLeiSongData/Kevin/

cd $basedir

if [ ! -d splited ]
then
	mkdir -p splited
	cd splited
	# awk '{if(NR>1) print}' ../REVEAL_TLCN_All.QCed.traw | split -l 10000 -d -a 2 - REVEAL_TLCN_All.QCed.traw__
	cat ${basedir}TOPMed_Imputation/processed/processed.traw | parallel --header : --pipe -N10000 'cat | gzip > processed.traw__{#}.gz'
fi

if [ -f ${basedir}code/run.swarm ]
then
  rm  ${basedir}code/run.swarm
fi

rscript=" ${basedir}code/ordinal.lr.R"

splited_file_dir="${basedir}splited"
out_dir="${basedir}result"


# ordinal.lr <- function(snp.traw.file = "/data/songl5/BinBB/zhiwei/GWAS_and_HBV_DNA_Load/scripts/gwas/splited/REVEAL_TLCN_All.QCed.traw__1.gz",  pheno.Rdata.file = "/data/songl5/BinBB/zhiwei/GWAS_and_HBV_DNA_Load/scripts/pheno.Rdata", out.file = "/data/songl5/BinBB/zhiwei/GWAS_and_HBV_DNA_Load/scripts/gwas/out/REVEAL_TLCN_All.QCed.traw__1.txt"){

ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
  ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
  printf "Rscript $rscript $ifile ${out_dir}/$ofile\n" >> ${basedir}code/run.swarm
  echo $ifile
  echo $ofile
done

if [ -d ${basedir}logs ]
then
	rm -r ${basedir}logs
fi

mkdir ${basedir}logs
cd ${basedir}logs

swarm -f ${basedir}code/run.swarm -g 16 --module R/4.1.0 --time=10:00:00
