####差异分析####
####免疫评分####
setwd("TCGA-PAAD")
setwd("Immune_DEG") 
library(BiocManager)
#BiocManager::install('DESeq2')
library(DESeq2)
library(tidyverse)
#TCGA差异分析用counts来做 因为是把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#整理分组信息
x <- "ImmuneScore"
med <- as.numeric(median(estimate[,x]))
estimate <- as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))

conditions=data.frame(sample=colnames(counts_01A),
                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

#差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

#开始差异分析
dds <- DESeq(dds)
#这句很重要
resultsNames(dds)
#提取结果
res <- results(dds)
save(res,file="DEG_ImmuneScore.Rda")

####热图绘制####
DEG <- as.data.frame(res)
#读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#下载pheatmap包 
#install.packages("pheatmap")
library(pheatmap)
#提取差异基因表达谱
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
#设置分组信息
annotation_col <- conditions
#对exp_diff 列的顺序进行处理
a <- filter(annotation_col,group == 'high')
b <- filter(annotation_col,group == 'low')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)
#开始画图
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
#保存图片 调整大小
dev.off()#关闭画板

####基质评分####
setwd("TCGA-PAAD")
setwd("Stromal_DEG")
library(BiocManager)
library(DESeq2)
library(tidyverse)
#TCGA差异分析用counts来做 因为是把01A患者分组做差异分析所以读取01A
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#因为是用免疫评分分组所以读取ESTIMATE_result
estimate <- read.table("ESTIMATE_result.txt", sep = "\t",row.names = 1,check.names = F,header = T)
#整理分组信息
x <- "StromalScore"
med <- as.numeric(median(estimate[,x]))
estimate <- as.data.frame(t(estimate))
identical(colnames(counts_01A),colnames(estimate))

conditions=data.frame(sample=colnames(counts_01A),
                      group=factor(ifelse(estimate[x,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")
#差异分析准备工作
dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

#开始差异分析
dds <- DESeq(dds)
resultsNames(dds)
#提取结果
res <- results(dds)
save(res,file="DEG_StromalScore.Rda")

####热图绘制####
DEG <- as.data.frame(res)
#读取表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

library(pheatmap)
#提取差异基因表达谱
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
c <- rbind(a,b)
d <- rownames(c)
exp_diff <- exp[d,]
#设置分组信息
annotation_col <- conditions
#对exp_diff 列的顺序进行处理
a <- filter(annotation_col,group == 'high')
b <- filter(annotation_col,group == 'low')
exp_diff_high <- exp_diff[,rownames(a)]
exp_diff_low <- exp_diff[,rownames(b)]
exp_diff <- cbind(exp_diff_high,exp_diff_low)
#开始画图
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows = T,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
#保存图片 调整大小
dev.off()#关闭画板

####将两次差异分析的差异基因取交集####
setwd("TCGA-PAAD")
setwd("Immune_Stromal_DEG")
#打开DEG_ImmuneScore.rda
DEG <- as.data.frame(res)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#提取上下调基因
library(tidyverse)
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
write.csv(a, file = "Immune_up.csv")
write.csv(b, file = "Immune_down.csv")

#打开DEG_StromalScore.rda
DEG <- as.data.frame(res)
#添加上下调信息
logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
#提取上下调基因
a <- filter(DEG,change == 'UP')
b <- filter(DEG,change == 'DOWN')
write.csv(a, file = "Stromal_up.csv")
write.csv(b, file = "Stromal_down.csv")
#仙桃画图韦恩图
