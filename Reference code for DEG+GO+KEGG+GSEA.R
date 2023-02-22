#### GEO数据库的使用 ####
# GEO网站：https://www.ncbi.nlm.nih.gov/geo/

#### GSE84402 ####
setwd("GSE84402")

###加载R包
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)

###下载数据，如果文件夹中有会直接读入
chooseBioCmirror()
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)

###提取子集
gset[[1]]

#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))
##2.2 通过exprs函数获取表达矩阵并校正
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###数据校正
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
#使用R包转换id
index = gset[[1]]@annotation
if(!require("hgu133plus2.db"))
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(29:30)]
write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####GEO手动注释####
####GSE31056####
setwd("GSE31056")
###加载R包
library(tidyverse)
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
chooseBioCmirror()
gset = getGEO('GSE31056', destdir=".", AnnotGPL = F, getGPL = F)
#有时会报错  Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
class(gset)
###提取子集
gset[[1]]
#读取表达谱
exp <- exprs(gset[[1]])
#把表达谱转为数据框格式
exp <- as.data.frame(exp)
##转换id
#读取GPL文件
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]
exp1 <- cbind(GPL,exp)
exp1 <- exp1[!duplicated(exp1$SYMBOL),]
rownames(exp1) <- exp1$SYMBOL
exp1 <- exp1[,-(1:5)]
write.table(exp1, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####中秋快乐####
setwd("GSE84402")
###加载R包
library(tidyverse)
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))

##读取上节课整理好的表达数据exp##
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#差异分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

##热图##
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####GEO三大富集分析
setwd("GSE84402")
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.table("deg_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
deg <- deg %>% filter(change != "stable")

DEG <- deg
DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result

#KEGG
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result

#GSEA
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.0.entrez.gmt"    #c2.all.v7.0.entrez.gmt 或 c5.all.v7.0.entrez.gmt
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,2]
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)

set.seed(1)
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
KEGG_result_df <- as.data.frame(KEGG)
write.table(KEGG_result_df,file="GSEA_MSigDb_C5_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_SPP1.rda")
