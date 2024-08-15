

library(gridExtra)
library(limma) 
library(ggpubr)
library(pROC)# 用于ROC曲线分析
# 设置工作目录和文件路径
expFile = "geneMatrix.txt" 
conFile = "s1.txt" 
treatFile = "s2.txt" 
geneFile = "intersect.txt"

# 读取基因表达数据并进行预处理
rt = read.table(expFile,header = TRUE, sep = "\t", check.names =FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1] 
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp),colnames(exp))
data = matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)
rt =avereps(data)
data = normalizeBetweenArrays(rt)



con =read.table(conFile, header = FALSE, sep = "\t", check.names=FALSE)
treat = read.table(treatFile,header =FALSE, sep = "\t", check.names =FALSE)	
conData = data[,as.vector(con[,1])]
treatData = data[,as.vector(treat[,1])]
data = cbind(conData,treatData) 
conNum = ncol(conData)
treatNum =ncol(treatData)

Type = c(rep("con",conNum),rep("treat", treatNum))
outData = rbind(id = paste0(colnames(data),"_",Type),data)
write.table(outData, file = "test.normalize.txt", sep = "\t", quote =FALSE, col.names =FALSE)
nrow(data)
ncol(data)


geneRT = read.table (geneFile,header =FALSE,sep = "\t",check.names=FALSE)	
data = data[as.vector(geneRT[,1]),,drop =FALSE]
Type = c(rep("Con",conNum), rep("Treat",treatNum)) 
my_comparisons = list()
my_comparisons[[1]] = levels(factor (Type) ) 
newGeneLists = c() 
outTab = data.frame ()
for (i in row.names(data)) {
  rtl = data.frame (expression = data[i,], Type =Type)
  boxplot = ggboxplot (rtl,x = "Type",y = "expression", color = "Type",
                     xlab = "",
                     ylab = paste(i,"expression"), 
                     legend.title = "",
                     #palette = c("blue", "orange"),
                     #palette = c("#C388FE", "#77C034"),
                     palette = c("#0f86a9", "#FC8452"),
                     add = "jitter") +
  stat_compare_means(comparisons = my_comparisons)
pdf (file = paste0("boxplot.", i, ".pdf"), width = 3, height =4.5) 
print (boxplot) 
dev.off()
}



exp <- mRNA_exp[c("ENSG00000180871","ENSG00000197561"),]
exp=t(exp)
exp=data.frame(exp)
ggscatter(exp, x = colnames(exp)[1], y = colnames(exp)[2],
                            size = 1,
                            add = "reg.line",  # 添加回归线
                             add.params = list(color = "#0AA1FF", fill = "#a5dff9", size = 1),  # 自定义回归线的颜色
                           conf.int = TRUE  # 添加置信区间
                    ) +
    stat_cor(method = "spearman",  label.x =10, label.y = 7,label.sep = "\n") +
   xlab("CXCR2") +
  ylab("")###NET genes



library(pheatmap)
logFCfilter=1        
adjPfilter=0.05      
gene="ENSG00000163735"      
expFile="GSE66407.normalize.txt"    
setwd("")     
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
rt=rt[,Type=="Treat"]
low=rt[gene,]<=median(rt[gene,])
high=rt[gene,]>median(rt[gene,])
lowRT=rt[,low]
highRT=rt[,high]
conNum=ncol(lowRT)        
treatNum=ncol(highRT)    
rt=cbind(lowRT,highRT)



Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="all.txt",sep="\t",quote=F)


diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adjPfilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.txt",sep="\t",quote=F,col.names=F)




geneNum=20     
diffUp=diffSig[diffSig$logFC>0,]
diffDown=diffSig[diffSig$logFC<0,]
geneUp=row.names(diffUp)
geneDown=row.names(diffDown)
if(nrow(diffUp)>geneNum){geneUp=row.names(diffUp)[1:geneNum]}
if(nrow(diffDown)>geneNum){geneDown=row.names(diffDown)[1:geneNum]}
hmExp=rt[c(geneUp,geneDown),]

Type=c(rep("Low",conNum),rep("High",treatNum))
Type=factor(Type, levels=c("Low", "High"))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
colnames(Type)[1]=gene

pdf(file="heatmap1.pdf", width=7.5, height=5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


hmExp=rt[unique(c(gene,geneUp,geneDown)),]
hmExpOut=rbind(id=colnames(hmExp),hmExp)
write.table(hmExpOut,file="corInput.txt",sep="\t",quote=F,col.names=F)


pdf(file="vol.pdf",width=5,height=5)
xMax=4      
yMax=10    
plot(allDiff$logFC, -log10(allDiff$adj.P.Val), ylab="-log10(adj.P.Val)",xlab="logFC", col="grey",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1)
diffSub=subset(allDiff, adj.P.Val<adjPfilter & logFC>logFCfilter)
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="red",cex=1.2)
diffSub=subset(allDiff, adj.P.Val<adjPfilter & logFC<(-logFCfilter))
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="green",cex=1.2)
abline(v=0,lty=2,lwd=3)
dev.off()

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

gene="CXCR1"     
expFile="GSE66407.normalize.txt"         
gmtFile="c2.cp.kegg.Hs.symbols.gmt"    
setwd("")             




rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,Type=="Treat",drop=F]
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]   
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=meanH-meanL
logFC=sort(logFC, decreasing=T)
genes=names(logFC)
gmt=read.gmt(gmtFile)
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)


termNum=5    
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in high expression group")
  pdf(file="GSEA.highExp.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}


termNum=5    
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in low expression group")
  pdf(file="GSEA.lowExp.pdf", width=6.5, height=5.5)
  print(gseaplot)
  dev.off()
}




