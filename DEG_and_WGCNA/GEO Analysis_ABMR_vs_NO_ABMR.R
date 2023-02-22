rm(list = ls())

# 在调用as.data.frame的时，将stringsAsFactors设置为FALSE可以避免character类型自动转化为factor类型
options(stringsAsFactors = F) 

library(tidyverse)

load('GSE147089_eSet.rda')  ## 载入数据

dat <- gset[[1]] %>% exprs() # gset[[1]]现在是一个对象，取这个对象通过看说明书知道要用exprs这个函数
dim(dat) # 看一下dat这个矩阵的维度

# 预览dat这个矩阵的1至4行和1至4列，逗号前为行，逗号后为列
dat[1:4,1:4] 

# 箱型图查看表达量数据是否已进行标准化
boxplot(dat,las=2) 
dev.off()

# 通过查看说明书知道取对象gset[[1]]里的临床信息用pData()
pd <- gset[[1]] %>% pData() 

# 挑选一些感兴趣的临床表型。
# group_list <- str_split(pd$title,' ',simplify = T)[,4]
group_list <- as.character(pd$`phenotype:ch1`)
table(group_list)
dat[1:4, 1:4] 
pd[1:4, 1:4]
# ABMR与NON_ABMR分组
group_list1 <- ifelse(str_detect(dat$group_list, "No "), "NO_ABMR" , "ABMR")
table(group_list1)

# GPL570
f <- 'gpl570_probe2gene.rda'
if (!file.exists(f)) {
  # BiocManager::install('hgu133plus2.db')
  library(hgu133plus2.db)
  ls("package:hgu133plus2.db")
  probe2gene <- toTable(hgu133plus2SYMBOL)
  head(probe2gene)
  save(probe2gene, file = f)
}

load(file='gpl570_probe2gene.rda')
ids <- probe2gene 
head(ids)
dim(ids)
ids <- ids[ids$symbol != '',] # 去除ids中symbol为空的行（探针没有对应的注释基因）
ids <- ids[ids$probe_id %in%  rownames(dat),] # 提取出在dat中出现过的探针所在的行

dat[1:4,1:4]
ids[1:4,] # 可以看到经过上一步提取操作后ids与dat的行名顺序对齐了

# 取出dat中探针有注释的行，过滤掉不表达的基因？
dat <- dat[ids$probe_id,] 
# ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids$median <- apply(dat,1,median) 
# 对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids <- ids[order(ids$median, decreasing = T),]

nrow(ids)
# 将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果
ids <- ids[!duplicated(ids$symbol),]
nrow(ids)

# 新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
dat <- dat[ids$probe_id,] 
dim(dat)

# 把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat) <- ids$symbol 
dat[1:4,1:4]


save(dat, group_list, group_list1, 
     file = 'step1-output_ABMR_vs_NO_ABMR.rda')


### step2 check ----

rm(list = ls())
options(stringsAsFactors = F)
load(file = 'step1-output_ABMR_vs_NO_ABMR.rda')
table(group_list)
table(group_list1)
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat <- t(dat) # 画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
dat <- as.data.frame(dat) # 将matrix转换为data.frame
# dat <- cbind(dat,group_list) # cbind横向追加，即将分组信息追加到最后一列
# 追加一列ABMR与NON_ABMR的分类
dat <- cbind(dat,group_list1) 

library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 

# The variable group_list (index = 20858) needs to be removed before PCA analysis.
# 现在dat最后一列是group_list，需要重新赋值给一个dat.pca, 这个矩阵是不含有分组信息的
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list1, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA_grouped_by_ABMR.pdf')


rm(list = ls())
load(file = 'step1-output_ABMR_vs_NO_ABMR.rda') 
dat[1:4,1:4] 

# apply按行（'1'是按行取，'2'是按列取）取每一行的方差，从小到大排序，取最大的1000个
cg <- names(tail(sort(apply(dat,1,sd)),200))
library(pheatmap)
# 将提取出来的1000个基因所在的每一行取出，组合起来为一个新的表达矩阵，绘热图
dat[cg,] %>% pheatmap(show_colnames =F,show_rownames = F) # 表现不好
n <- dat[cg,] %>% t() %>% scale() %>% t() # 'scale()'对log-ratio数值进行归一化
n[n > 2] = 2
n[n < -2] = -2
n[1:4,1:4]
pheatmap(n, show_colnames =F, show_rownames = F)
ac <- data.frame(g = group_list1)
rownames(ac) <- colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息

pheatmap(n, show_colnames =F, show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top200_sd.pdf')
dev.off()

### step3 DEG ----
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
load(file = 'step1-output_ABMR_vs_NO_ABMR.rda')
# 每次都要检测数据
dat[1:4,1:4] 
table(group_list)
table(group_list1)

#通过为每个数据集绘制箱形图，比较数据集中的数据分布
boxplot(dat[1,]~group_list) #按照group_list分组画箱线图
boxplot(dat[1,]~group_list1) #按照group_list分组画箱线图
dev.off()

#定义一个函数g，函数为{}里的内容
bp <- function(gene, group){
  library(ggpubr)
  df <- data.frame(gene = gene, 
                   stage= group)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}

bp(dat[1,], group_list) ## 调用上面定义好的函数，避免同样的绘图代码重复多次敲。
bp(dat[1,], group_list1)
dim(dat)

library(limma)
design <- model.matrix(~factor( group_list1 ))
fit <- lmFit(dat, design)
fit <- eBayes(fit)

## 上面是limma包用法的一种方式 
options(digits = 4) #设置全局的数字有效位数为 4
#topTable(fit,coef=2,adjust='BH') 
topTable(fit, coef=2, adjust='BH') 
## 但是上面的用法做不到随心所欲的指定任意两组进行比较

design <- model.matrix(~0+factor(group_list1))
colnames(design) <- levels(factor(group_list1))
head(design)
exprSet <- dat
rownames(design) <- colnames(exprSet)
design

contrast.matrix <- makeContrasts("ABMR-NO_ABMR", 
                                 levels = design)
contrast.matrix ##这个矩阵声明，我们要把 ABMR 跟 NO_ABMR 进行差异分析比较

deg <- function(exprSet, design, contrast.matrix){
  ##step1
  fit <- lmFit(exprSet, design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

deg <- deg(exprSet, design, contrast.matrix)

head(deg)

save(deg, file = 'deg_ABMR_vs_NO_ABMR.rda')
# > rownames(deg)[254]
# [1] "STAT1"
# > rownames(deg)[81]
# [1] "IRF1"


## for volcano（火山图） 
if(T){
  nrDEG <- deg
  head(nrDEG)
  attach(nrDEG)
  plot(logFC,-log10(P.Value))
  library(ggpubr)
  df <- nrDEG
  df$v <- -log10(P.Value) # df新增加一列'v',值为-log10(P.Value)
  ggscatter(df, x = "logFC", y = "v", size=0.5)
  
  df$g <- ifelse(df$P.Value>0.01,'stable', # if 判断：如果这一基因的P.Value>0.01，则为stable基因
              ifelse( df$logFC >2,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再 if 判断：如果logFC >1.5,则为up（上调）基因
                      ifelse( df$logFC < -2,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
  )
  table(df$g)
  df$name <- rownames(df)
  head(df)
  ggscatter(df, x = "logFC", y = "v", size=0.5, color = 'g')
  ggscatter(df, x = "logFC", y = "v", color = "g",size = 0.5,
            label = "name", repel = T,
            #label.select = rownames(df)[df$g != 'stable'] ,
            label.select = c(head(rownames(deg),10), 'IRF1'), # 挑选一些基因在图中显示出来
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
  ggsave('volcano_ABMR_vs_NO_ABMR.pdf')
  
  ggscatter(df, x = "AveExpr", y = "logFC",size = 0.2)
  df$p_c = ifelse(df$P.Value<0.001,'p<0.001',
                  ifelse(df$P.Value<0.01,'0.001<p<0.01','p>0.01'))
  table(df$p_c )
  ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
            palette = c("green", "red", "black") )
  ggsave('MA_ABMR_vs_NO_ABMR.pdf')
  
  
}

## for heatmap 
if(T){ 
  load(file = 'step1-output_ABMR_vs_NO_ABMR.rda')
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list)
  x <- deg$logFC # deg取logFC这列并将其重新赋值给x
  names(x) <- rownames(deg) # deg取probe_id这列，并将其作为名字给x
  # 对x进行从小到大排列，取前100及后100，并取其对应的探针名，作为向量赋值给cg
  cg <- c(names(head(sort(x),100)), 
          names(tail(sort(x),100)))
  library(pheatmap)
  library(tidyverse)
  # 对dat按照cg取行，所得到的矩阵来画热图
  dat[cg,] %>% pheatmap(show_colnames =F, show_rownames = F) 
  
  ## 通过“scale”对log-ratio数值进行归一化，现在的dat是行名为探针，列名为样本名，
  ## 由于scale这个函数应用在不同组数据间存在差异时，需要行名为样本，因此需要用t(dat[cg,])来转换，最后再转换回来
  n <- dat[cg,] %>% t() %>% scale() %>% t() 
  
  n[n > 2] = 2
  n[n < -2] = -2
  n[1:4, 1:4]
  pheatmap(n, show_colnames =F, show_rownames = F)
  ac <- data.frame(g=group_list1)
  rownames(ac) <- colnames(n) #将ac的行名也就分组信息 给到n的列名，即热图中位于上方的分组信息
  pheatmap(n,
           show_colnames = F,
           show_rownames = F,
           cluster_cols = F, 
           annotation_col = ac, 
           filename = 'heatmap_top200_DEG.pdf') #列名注释信息为ac即分组信息
  
}
dev.off()

write.csv(deg, file = 'deg_ABMR_vs_NO_ABMR.csv')


### step4 anno-go-kegg ----


rm(list = ls())  ## 魔幻操作，一键清空~
load(file = 'deg_ABMR_vs_NO_ABMR.rda')
head(deg)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t <- 1.0
deg$g <- ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
deg$symbol <- rownames(deg)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), 
           fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG <- deg
head(DEG)

DEG <- merge(DEG, df, by.y='SYMBOL', by.x='symbol')
head(DEG)
save(DEG,
     file = 'anno_DEG_ABMR_vs_NO_ABMR.rda')

gene_up <- DEG[DEG$g == 'UP','ENTREZID'] 
gene_down <- DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff <- c(gene_up, gene_down)
gene_all <- as.character(DEG[ ,'ENTREZID'] )
data(geneList, package = "DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)
dev.off()

geneList <- DEG$logFC
names(geneList) <- DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)


## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  dotplot(kk.up );ggsave('kk.up.dotplot_ABMR_vs_NO_ABMR.pdf')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down );ggsave('kk.down.dotplot_ABMR_vs_NO_ABMR.pdf')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  dotplot(kk.diff );ggsave('kk.diff.dotplot_ABMR_vs_NO_ABMR.pdf')
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  
  down_kegg <- kegg_down_dt[kegg_down_dt$pvalue<0.05,] ; down_kegg$group = -1
  up_kegg <- kegg_up_dt[kegg_up_dt$pvalue<0.05,] ; up_kegg$group = 1
  source('functions.R')
  g_kegg <- kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg, filename = 'kegg_up_down_ABMR_vs_NO_ABMR.pdf')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,] ; down_kegg$group = -1
  up_kegg <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,] ; up_kegg$group = 1
  
  g_kegg <- kegg_plot(up_kegg, down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea_ABMR_vs_NO_ABMR.pdf')
  
  
}

### GO database analysis 
### 做GO数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
{
  
  g_list <- list(gene_up = gene_up, 
                 gene_down = gene_down, 
                 gene_diff = gene_diff)
  
  if(F){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP', 'MF', 'CC') , function(ont) {
        cat(paste('Now process ', ont ))
        ego <- enrichGO(gene          = gene,
                        universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file = 'go_enrich_results_ABMR_vs_NO_ABMR.rda')
    
  }
  
  
  load(file = 'go_enrich_results_ABMR_vs_NO_ABMR.rda')
  
  n1 <- c('gene_up','gene_down','gene_diff')
  n2 <- c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn = paste0('dotplot_',n1[i],'_',n2[j],'_ABMR_vs_NO_ABMR','.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
}


### step5 anno-GSEA ----


library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

### 对 MigDB中的全部基因集 做GSEA分析。
# http://www.bio-info-trainee.com/2105.html
# http://www.bio-info-trainee.com/2102.html 
{
  load(file = 'anno_DEG_ABMR_vs_NO_ABMR.rda')
  geneList <- DEG$logFC
  names(geneList) <- DEG$symbol
  geneList <- sort(geneList, decreasing = T)
  #选择gmt文件（MigDB中的全部基因集）
  d <- '../MsigDB/symbols'
  gmts <- list.files(d, pattern = 'all')
  gmts
  #GSEA分析
  library(GSEABase) # BiocManager::install('GSEABase')
  ## 下面使用lapply循环读取每个gmt文件，并且进行GSEA分析
  ## 如果存在之前分析后保存的结果文件，就不需要重复进行GSEA分析。
  f <- 'gsea_results_ABMR_vs_NO_ABMR.rda'
  if(!file.exists(f)){
    gsea_results <- lapply(gmts, function(gmtfile){
      # gmtfile=gmts[2]
      geneset <- read.gmt(file.path(d,gmtfile)) 
      print(paste0('Now process the ',gmtfile))
      egmt <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
      head(egmt)
      # gseaplot(egmt, geneSetID = rownames(egmt[1,]))
      
      return(egmt)
    })
    # 上面的代码耗时，所以保存结果到本地文件
    save(gsea_results,file = f)
  }
  load(file = f)
  #提取gsea结果，熟悉这个对象
  gsea_results_list <- lapply(gsea_results, function(x){
    cat(paste(dim(x@result)),'\n')
    x@result
  })
  gsea_results_df <- do.call(rbind, gsea_results_list)
  gseaplot(gsea_results[[2]],'KEGG_CELL_CYCLE') 
  gseaplot(gsea_results[[2]],'FARMER_BREAST_CANCER_CLUSTER_6') 
  gseaplot(gsea_results[[5]],'GO_CONDENSED_CHROMOSOME_OUTER_KINETOCHORE') 
  
}


### step6 anno-GSVA ----

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

### 对 MigDB中的全部基因集 做GSVA分析。
## 还有ssGSEA, PGSEA
{
  load(file = 'step1-output_ABMR_vs_NO_ABMR.rda')
  # 每次都要检测数据
  dat[1:4,1:4]  
  
  X <- dat
  table(group_list1)
  ## Molecular Signatures Database (MSigDb) 
  d='../MSigDB/symbols/'
  gmts <- list.files(d,pattern = 'all')
  gmts
  library(GSVA) # BiocManager::install('GSVA')
  library(GSEABase)
  
  if(T){
    timestamp()
    es_max <- lapply(gmts, function(gmtfile){ 
      #gmtfile=gmts[8];gmtfile
      geneset <- getGmt(file.path(d, gmtfile))  
      es.max <- gsva(X, geneset, 
                     mx.diff = FALSE, verbose = FALSE, 
                     parallel.sz = 1)
      return(es.max)
    })
    timestamp()
    ## 上一步太久了，保存一下结果
    save(es_max, file = 'anno_gsva_es_max.rda')
    
    adjPvalueCutoff <- 0.001
    logFCcutoff <- log2(2)

    es_deg <- lapply(es_max, function(es.max){
      table(group_list1)
      dim(es.max)
      design <- model.matrix(~0+factor(group_list1))
      colnames(design) <- levels(factor(group_list1))
      rownames(design) <- colnames(es.max)
      design
      library(limma)
      # contrast.matrix <- makeContrasts(paste0(unique(group_list1), collapse = "-"),
      #                                  levels = design)
      contrast.matrix <- makeContrasts("ABMR-NO_ABMR",
                                       levels = design)
      
      contrast.matrix ##这个矩阵声明，我们要把ABMR跟NO_ABMR进行差异分析比较
      
      deg <- function(es.max,design,contrast.matrix){
        ##step1
        fit <- lmFit(es.max, design)
        ##step2
        fit2 <- contrasts.fit(fit, contrast.matrix) 
        ##这一步很重要，大家可以自行看看效果
        
        fit2 <- eBayes(fit2)  ## default no trend !!!
        ##eBayes() with trend=TRUE
        ##step3
        res <- decideTests(fit2, p.value = adjPvalueCutoff)
        summary(res)
        tempOutput <- topTable(fit2, coef=1, n=Inf)
        nrDEG <- na.omit(tempOutput) 
        #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
        head(nrDEG)
        return(nrDEG)
      }
      
      re <- deg(es.max, design, contrast.matrix)
      nrDEG <- re
      head(nrDEG) 
      return(nrDEG)
    })

    } 
  
  gmts
  
  save(es_max, es_deg, 
       file = 'gsva_msigdb_ABMR_vs_NO_ABMR.rda')
  
  load(file = "gsva_msigdb_ABMR_vs_NO_ABMR.rda")
  
  library(pheatmap)
  lapply(1:length(es_deg), function(i){
    # i=2
    print(i)
    dat <- es_max[[i]]
    df <- es_deg[[i]]
    df <- df[df$P.Value<0.01 & abs(df$logFC) > 0.3,]
    print(dim(df))
    if(nrow(df)>5){
      n <- rownames(df)
      dat <- dat[match(n,rownames(dat)),]
      ac <- data.frame(g=group_list1)
      rownames(ac) <- colnames(dat)
      rownames(dat) <- substring(rownames(dat),1,50)
      pheatmap::pheatmap(dat, 
                         fontsize_row = 8,height = 11,
                         annotation_col = ac,show_colnames = F,
                         filename = paste0('gsva_',strsplit(gmts[i],'[.]')[[1]][1],'_ABMR_vs_NO_ABMR','.pdf'))
      
    }
  })
  dev.off()
  
  adjPvalueCutoff <- 0.001
  logFCcutoff <- log2(2)
  df <- do.call(rbind ,es_deg)
  es_matrix <- do.call(rbind ,es_max)
  df <- df[df$P.Value<0.01 & abs(df$logFC) > 0.5,]
  write.csv(df, file = 'GSVA_DEG_ABMR_vs_NO_ABMR.csv')
}
