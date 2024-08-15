
library(gwasglue) 
library(dplyr)
library(gassocplot) 
library(coloc)
library(locuscomparer)
eqt1ID="eqtl-a-ENSG00000163464"
outcomeID="ebi-a-GCST004133"		
geneChr=2
geneStart=219027564	
geneEnd=219031685 	

chrpos <- paste0(geneChr,":",geneStart - 100000,"-",geneEnd +100000)
out <- ieugwasr_to_coloc(id1=eqt1ID,id2=outcomeID,chrompos=chrpos)
result <- coloc::coloc.abf(out[[1]],out[[2]])
k=out[[1]]
k=data.frame(k[[1]],k[[7]])
k=k[,c(2,1)]
colnames(k)=c("rsid","pval")

p=out[[2]]
p=data.frame(p[[1]],p[[7]])
p=p[,c(2,1)]
colnames(p)=c("rsid","pval")
write.table(k,file = "eqtl.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(p,file = "gwas.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
locuscompare(in_fn1 = "D:/中性粒细胞孟德尔/共定位/eqtl.tsv",in_fn2 ="D:/中性粒细胞孟德尔/共定位/gwas.tsv" )


