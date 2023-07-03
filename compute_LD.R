
#https://genome.sph.umich.edu/wiki/LocusZoom_Standalone#User-supplied_LD
#https://rdrr.io/cran/genetics/man/LD.html
#https://rdrr.io/cran/genetics/man/genotype.html
METALfile=paste0("../result/sig_snps_metalfile.txt")
locusdata=read.table(METALfile,header = T)
load("../result/sig_res2all.RData")
idx=match(locusdata$MarkerName,sig_res2$rsid)
locusdata$snp=sig_res2$SNP[idx]
write.table(locusdata$snp,file="../result/locusdata_snp.txt",row.names = F,col.names = F,quote=F)
write.table(locusdata[,c(1,5)],file="../result/locusdata_snpmap.txt",row.names = F,col.names = F,quote=F)
chr6data=data.table::fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/chr6_filter.bim")

#change locusdata from hg38 to hg19
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("../../tools/hg38ToHg19.over.chain")
locusdata$pos=NA
tmp=unlist(strsplit(locusdata$snp,":"))
locusdata$pos=as.numeric(tmp[seq(2,length(tmp),4)])
gr_dat=GRanges(seqnames = "chr6",ranges=IRanges(start=locusdata$pos,width=1))
tmp=as.data.frame(liftOver(gr_dat,chain))
write.table(tmp$start,file="../result/locusdata_hg19.txt",row.names = F,col.names = F,quote=F)

tmp1=read.table("../result/locusdata.bim",header=F)
tmp1$V4=tmp$start
write.table(tmp1,file="../result/locusdata.bim",row.names = F,col.names = F,quote=F)

tmp2=read.table("../result/locusdata.vcf")
tmp2$V7="PASS"
write.table(tmp2,sep="\t",file="../result/locusdata_pass.vcf",row.names = F,col.names = F,quote=F)

#check LD on chr6 peak snps
load("../result/info_res.RData")
idx = which(allres$P<5.9e-8)
sig_res=allres[idx,]
plot(sig_res$BP,-log10(sig_res$P),xlab="Position",ylab="-log10(P)",cex.axis=1.2,cex.lab=1.2)
chr6genotype=data.table::fread("/data/DCEGLeiSongData/Kevin/HBVloadGwas/TOPMed_Imputation/processed/chr6_filter.traw")
idx=match(sig_res$SNP,chr6data$V2)
locusgenotype=chr6genotype[idx,7:ncol(chr6genotype)]
tmp=cor(t(locusgenotype))
tmp=tmp^2
heatmap(tmp)
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
clusters=c("cluster1",rep("cluster2",8),rep("cluster3",10))
rownames(tmp)=colnames(tmp)=sig_res$SNP
drawheatmap=function(clusterrows=F,clustercolumns=F,roworder=NULL)
{
  ha = HeatmapAnnotation(cluster =clusters, 
                         col = list(cluster = c("cluster1" = "green", "cluster2"="darkblue","cluster3"="skyblue")),
                         show_annotation_name = T,
                         #gp = gpar(col = "black"),
                         annotation_name_gp=gpar(fontface = "bold",fontsize = 12),
                         # annotation_name_offset = unit(2, "cm"),
                         # annotation_name_rot = c(0, 0),
                         annotation_name_side = "right")
  
  ht_list = Heatmap(tmp, #name = "R^2",
                    heatmap_legend_param = list(title = expression("R"^2),title_gp = gpar(fontsize = 12, fontface = "bold")),
                    row_order = roworder,
                    col = colorRamp2(c(0, 1), c("white", "red")), 
                    #column_dend_height = unit(4, "cm"),
                    cluster_rows = clusterrows,
                    cluster_columns=clustercolumns,
                    column_dend_reorder=F,
                    row_names_gp=gpar(fontsize = 8, fontface = "bold"),
                    column_names_gp=gpar(fontsize = 8, fontface = "bold"),
                    row_dend_width=unit(2.5,"cm"),
                    top_annotation = c(ha),
                    show_row_dend=F,
                    show_column_names = T, show_row_names=T) 
  ht_list = draw(ht_list, heatmap_legend_side = "right",legend_title_gp = gpar(fontsize = 12, fontface = "bold"))
  return(ht_list)
}
pdf(file="../result/chr6_peaksnps_R2.pdf",width=8)
drawheatmap()
dev.off()

#conditional analysis
pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno.RData"
#for 10 pc
#pheno.Rdata.file = "/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/pheno10pc.RData"
load(pheno.Rdata.file)
tmp=unlist(strsplit(colnames(locusgenotype),"_"))
colnames(locusgenotype)=tmp[seq(2,length(tmp),2)]
rownames(locusgenotype)=sig_res$SNP
all(rownames(pheno)%in%colnames(locusgenotype))
idx=match(colnames(locusgenotype),rownames(pheno))
pheno=pheno[idx,]
library(MASS)
dat1=data.frame(snp=t(locusgenotype[which(rownames(locusgenotype)=="chr6:30752534:C:T"),]),pheno,
                snp1=t(locusgenotype[which(rownames(locusgenotype)=="chr6:31347583:G:A"),]),
                snp2=t(locusgenotype[which(rownames(locusgenotype)=="chr6:32616916:C:G"),]))
#dat1=data.frame(snp=t(locusgenotype[which(rownames(locusgenotype)=="chr6:30752534:C:T"),]),pheno)
fm = polr(formula = HBVDNA.grp ~ ., data = dat1, Hess = TRUE)
or = exp(cbind(OR = coef(fm), ci = confint.default(fm)))
ctable = coef(summary(fm))
ctable = cbind(ctable, "p_value" = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2)
ctable
# Value Std. Error     t value      p_value
# snp            -0.1885529 0.04668262  -4.0390393 5.367057e-05
# PC1            -0.3257649 1.97858391  -0.1646455 8.692230e-01
# PC2             4.4638750 1.92510441   2.3187703 2.040749e-02
# PC3            -1.2234477 1.75660200  -0.6964854 4.861249e-01
# PC4            -5.4242927 1.92472462  -2.8182175 4.829108e-03
# PC5            -3.1667139 2.11879477  -1.4945826 1.350234e-01
# SEX1            0.2261426 0.06504520   3.4766992 5.076269e-04
# AGE.grp(39,49] -0.1274991 0.08011884  -1.5913742 1.115254e-01
# AGE.grp(49,59] -0.2827505 0.08052812  -3.5112021 4.460850e-04
# AGE.grp(59,79] -0.6988602 0.11922471  -5.8617061 4.581353e-09
# snp1            0.2863769 0.05051293   5.6693785 1.433165e-08
# snp2           -0.3543211 0.06252970  -5.6664450 1.457906e-08
# 1|2            -1.2934212 0.09158345 -14.1228705 2.745750e-45
# 2|3             0.1177008 0.08884821   1.3247407 1.852572e-01
# 3|4             0.9346521 0.09082963  10.2901673 7.804065e-25
# 4|5             1.5199137 0.09398459  16.1719457 7.954593e-59

library(genetics)
x=unlist(locusgenotype[which(rownames(locusgenotype)=="chr6:30752534:C:T"),])
g1<- as.genotype.allele.count(x, alleles=c("C","T") )
x=unlist(locusgenotype[which(rownames(locusgenotype)=="chr6:31347583:G:A"),])
g2<- as.genotype.allele.count(x, alleles=c("G","A") )
x=unlist(locusgenotype[which(rownames(locusgenotype)=="chr6:32616916:C:G"),])
g3<- as.genotype.allele.count(x, alleles=c("C","G") )
data <- makeGenotypes(data.frame(g1,g2,g3))
LD(data)
