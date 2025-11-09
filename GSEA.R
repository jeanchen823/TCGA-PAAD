####GSEA####
#安装加载包
#GO KEGG 时已经装过 直接library即可
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db) #org.Hs.eg.db包主要注释基因:用于不同数据库ID间的转化
library(clusterProfiler)
DEG <- as.data.frame(res)%>% arrange(padj) 

DEG <- DEG %>% rownames_to_column("Gene")

geneList = DEG[,3]
names(geneList) = as.character(DEG[,'Gene'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

#GSEA基因集：https://zhuanlan.zhihu.com/p/504101161

msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "h.all.v7.0.symbols.gmt"    
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(1) #设置种子
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "GSEA_CXCL10_h.all.rda")
#绘图
#安装enrichplot
library(enrichplot)
#单个结果绘制
gseaplot2(gsea,1,color="red",pvalue_table = T)
#多个结果绘制
#A
gseaplot2(gsea, geneSetID = c(1,2,3,4,5,6,8,10), subplots = 1:3)
#B
gseaplot2(gsea, geneSetID = c(7,9,11,13,14,16,17), subplots = 1:3)
gseaplot2(gsea, geneSetID = 1:10, subplots = 1:3)
dev.off()

####换C7跑####
msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c7.all.v7.0.symbols.gmt"    
#读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

set.seed(1) #设置种子
gsea <-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
#转换成数据框
gsea_result_df <- as.data.frame(gsea)
save(gsea,gsea_result_df,file = "GSEA_CXCL10_c7.rda")
#绘图
#C
gseaplot2(gsea, geneSetID = 1:7, subplots = 1:3)
#D
gseaplot2(gsea,782,color="red",pvalue_table = T)
gseaplot2(gsea, geneSetID = 782, subplots = 1:3)
dev.off()