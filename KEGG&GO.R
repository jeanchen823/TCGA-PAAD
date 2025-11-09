####用交集差异基因做富集分析####
setwd("TCGA-PAAD")
setwd("FUJI_Immune_Stromal_DEG")
library(tidyverse)
library("BiocManager")
#安装加载包
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
#org.Hs.eg.db包主要注释人类基因:用于不同数据库ID间的转化
library(clusterProfiler)
#导入DEG_final.txt
#导入immune或stromal差异分析结果 均可
DEG <- as.data.frame(res)
DEG <- DEG[DEG_final$SYMBOL,]
DEG <- rownames_to_column(DEG,"SYMBOL")
genelist <- bitr(DEG$SYMBOL, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL")

####GO####
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_DEG_final.Rda")

####KEGG####
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
save(kk,kk_res,file = "KEGG_DEG_final.Rda")

#网络图
library(ggnewscale)
#install.packages("ggnewscale")
List = DEG$log2FoldChange
names(List)= DEG$ENTREZID
head(List)
List = sort(List,decreasing = T)
#GO
cnetplot(ego, foldChange = List, circular = TRUE, colorEdge = TRUE)
#KEGG
cnetplot(kk, foldChange = List, circular = TRUE, colorEdge = TRUE)
