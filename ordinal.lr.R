args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop("Two inputs are required: file names of genotype and output")
SNP_TRAW_FILE = args[1]
OUT_FILE = args[2]
ordinal.lr <- function(snp.traw.file = "/data/DCEGLeiSongData/Kevin/splited/processed.traw__1.gz",  pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/result/pheno.RData", out.file = "/data/DCEGLeiSongData/Kevin/result/processed.traw__1.txt"){
	# browser()
	load(pheno.Rdata.file)
	
	snp.traw = read.delim(file = snp.traw.file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
	# CHR     SNP     (C)M    POS     COUNTED ALT     1000_CHB.AxiomGT1_1000_CHB.AxiomGT1
	snp.info = snp.traw[, c("CHR", "SNP", "POS", "COUNTED", "ALT")]
	colnames(snp.info) = c("CHR", "SNP", "BP", "EFF", "REF")
	
	snp.data = snp.traw[, -(1:6)]
	col.names = gsub("^0_","",colnames(snp.data))
	
	colnames(snp.data) = col.names
	idx = match(rownames(pheno), col.names)
	snp.data = t(snp.data[, idx])
	colnames(snp.data) = snp.info$SNP
	
	# if (!require("MASS",character.only = TRUE)){
	# 	install.packages("MASS", repos='http://cran.us.r-project.org')
	# }
	library(MASS)
	# library(lmtest)
	# if (!require("epiDisplay",character.only = TRUE)){
	# 	install.packages("epiDisplay", repos='http://cran.us.r-project.org')
	# }
	library(epiDisplay)
	
	col.names.res = c("MAF.cases", "MAF.ctrls", "OR", "OR_CI_%95_Low", "OR_CI_%95_High", "SE", "P")
	res = matrix(NA, nrow(snp.info), length(col.names.res))
	colnames(res) = col.names.res
	
	for(i in 1:ncol(snp.data)){
	 # for(i in 1:100){
		cur.data = data.frame(SNP = snp.data[, i], pheno, check.names = FALSE, stringsAsFactors = FALSE)
		
		idx.cases = cur.data$HBVDNA.grp != 1
		x = cur.data[idx.cases, "SNP"]
		res[i, "MAF.cases"] = sum(x, na.rm = TRUE)/(2*sum(!is.na(x)))
		
		idx.ctrls = !idx.cases
		x = cur.data[idx.ctrls, "SNP"]
		res[i, "MAF.ctrls"] = sum(x, na.rm = TRUE)/(2*sum(!is.na(x)))


		m1 = tryCatch(
				expr = {
					polr(formula = HBVDNA.grp ~ ., data = cur.data, Hess = TRUE)
				},
				error = function(e){ 
					return(NULL)
				}
			)
			
		if(is.null(m1)){
			cat(i, "\t", colnames(snp.data)[i], "\t", "Error!\n")
			next
		}

		ctable = coef(summary(m1))
		ctable = cbind(ctable, "p_value" = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)
		res[i, "SE"] = ctable["SNP", "Std. Error"]
		res[i, "P"] = ctable["SNP", "p_value"]
		or = exp(cbind(OR = coef(m1), ci = confint.default(m1)))
		res[i, c("OR", "OR_CI_%95_Low", "OR_CI_%95_High")] = or["SNP", ]
		
		#cat(i, "\t", colnames(snp.data)[i], "\t", "Success!\n")
		# browser()
	}
	
	res = data.frame(snp.info, res, check.names = FALSE, stringsAsFactors = FALSE)
	idx = is.na(res$P)
	if(sum(idx) > 0)
		res = res[!idx, ]
	
	idx = which(res[, "MAF.cases"] > 0.5)
	if(length(idx) > 0)
		res[idx, "MAF.cases"] = 1 - res[idx, "MAF.cases"]
	
	idx = which(res[, "MAF.ctrls"] > 0.5)
	if(length(idx) > 0)
		res[idx, "MAF.ctrls"] = 1 - res[idx, "MAF.ctrls"]
	
	write.table(res, file = out.file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

}

ordinal.lr(snp.traw.file = SNP_TRAW_FILE, out.file = OUT_FILE)






