####ESTIMATE####
#计算患者免疫评分与肿瘤纯度# 基质组分 免疫组分 肿瘤组分 肿瘤纯度
setwd("TCGA-PAAD")
setwd("ESTIMATE")  #设置工作目录
#安装包
library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#计算免疫评分
filterCommonGenes(input.f = "tpms01A_log2.txt",   #输入文件名
                  output.f = "tpms01A_log2.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("tpms01A_log2.gct",   #刚才的输出文件名
              "tpms01A_log2_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#提取结果并整理
ESTIMATE_result <- read.table("tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]   
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
#保存结果
write.table(ESTIMATE_result, file = "ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F) 