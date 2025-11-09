####TCGA-PAAD 数据下载####
setwd('TCGA-PAAD')
setwd("TCGAdata")
library(tidyverse)
library(BiocManager)
library(TCGAbiolinks)
chooseBioCmirror()
cancer_type = "TCGA-PAAD"
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file = "paad.gdc_2022.rda") 

#提取counts 
load("paad.gdc_2022.rda")#导入文件，rda格式文件也可直接从文件夹双击导入
load("gene_annotation_2022.rda")#导入gene注释文件
table(gene_annotation_2022$type)
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES
counts <- counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames("symbol") 
table(counts$type)#（注：可通过table(counts$type)查看基因类型）#lncRNA
counts <- counts[counts$type == "protein_coding",]
#counts <- counts[counts$type == "lncRNA",]
counts <- counts[,-c(1,ncol(counts))]
colnames(counts) <- substring(colnames(counts),1,16)
counts <- counts[,!duplicated(colnames(counts))]
table(substring(colnames(counts),14,16))
# 保留01A  （注：可通过table(substring(colnames(counts),14,16))查看样本类型）
counts01A <- counts[,substring(colnames(counts),14,16) == c("01A")]
# 保留11A
counts11A <- counts[,substring(colnames(counts),14,16) == c("11A")]
table(substring(colnames(counts01A),14,16))
table(substring(colnames(counts11A),14,16))

####tpms####
#和counts基本一模一样
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames("symbol") 
# 保留mRNA （注：可通过table(tpms$type)查看基因类型）
tpms <- tpms[tpms$type == "protein_coding",]
tpms <- tpms[,-c(1,ncol(tpms))]
# 把TCGA barcode切割为16位字符,并去除重复样本
colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
# 保留01A  （注：可通过table(substring(colnames(tpms),14,16))查看样本类型）
tpms01A <- tpms[,substring(colnames(tpms),14,16) == c("01A")]
# 保留11A
tpms11A <- tpms[,substring(colnames(tpms),14,16) == c("11A")]

#判断counts和tpms的行列名是否一致
identical(rownames(counts01A),rownames(counts11A))
identical(rownames(tpms01A),rownames(tpms11A))
identical(rownames(counts01A),rownames(tpms01A))
identical(colnames(counts01A),colnames(tpms01A))
identical(colnames(counts11A),colnames(tpms11A))
#保存counts和tpms数据
write.table(counts01A,"counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts11A,"counts11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A,"tpms01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A,"tpms11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#cbind和rbind 合并 col row
#cbind之前需要确认两个数据框的行名
counts <- cbind(counts01A,counts11A)
tpms <- cbind(tpms01A,tpms11A)
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms,"tpms.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####tpms_log2####
range(tpms)#查看数据范围
range(tpms01A)
range(tpms11A)
tpms_log2 <- log2(tpms+1)#log2转换 为什么要加1
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)
#保存log2转换后的数据
write.table(tpms_log2,"tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#表达谱整理完毕
