####整理TCGA临床信息####
setwd("TCGA-PAAD")
setwd("clinical")
library(tidyverse)
load("paad.gdc_2022.rda")
#提取临床信息
clinical <- as.data.frame(expquery2@colData) %>%   
  .[!duplicated(.$sample),]
#提取需要的临床信息数据
clinical <-clinical[,c("gender","age_at_index","ajcc_pathologic_stage",
                       "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m")]

class(clinical$gender)
class(clinical$age_at_index)
class(clinical$ajcc_pathologic_stage)
class(clinical$ajcc_pathologic_t)
class(clinical$ajcc_pathologic_n)
class(clinical$ajcc_pathologic_m)

table(clinical$gender)
table(clinical$age_at_index)
table(clinical$ajcc_pathologic_stage)
table(clinical$ajcc_pathologic_t)
table(clinical$ajcc_pathologic_n)
table(clinical$ajcc_pathologic_m)

clinical$ajcc_pathologic_stage <- gsub("A","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub("B","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_t <- gsub("a","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub("b","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_m <- gsub("a","",clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_m <- gsub("b","",clinical$ajcc_pathologic_m)

#提取01A临床数据
rownames(clinical) <- substring(rownames(clinical),1,16)

##将基因表达谱和临床数据合并并保存
exp01A <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

clinical01A <- clinical[colnames(exp01A),]   
exp01A <- exp01A %>% t() %>% as.data.frame()

identical(rownames(clinical01A),rownames(exp01A))

clinical.expr01A <- cbind(clinical01A,exp01A)

write.table(clinical.expr01A,"clinical.expr01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

##将ESTIMATE_result和临床数据合并并保存
ESTIMATE_result <- read.table("ESTIMATE_result.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

identical(rownames(clinical01A),rownames(ESTIMATE_result))

clinical.ESTIMATE_result01A <- cbind(clinical01A,ESTIMATE_result)

write.csv(clinical.ESTIMATE_result01A,file = "clinical.ESTIMATE_result01A.csv")
#解螺旋：https://www.helixlife.cn/class