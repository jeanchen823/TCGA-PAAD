####BTK在肿瘤样本与正常样本中的表达差异####
####柱状图####
setwd("TCGA-PAAD")
setwd("CXCL10")
library(tidyverse)
tpms01A_log2 <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms11A_log2 <- read.table("tpms11A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "CXCL10"#以后修改这里即可
a <- tpms01A_log2[gene,]
b <- tpms11A_log2[gene,]
#t转换
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "CXCL10_01A.csv")
write.csv(b, file = "CXCL10_11A.csv")

####配对图绘制####
tpms01A_log2 <- tpms01A_log2 %>% t() %>% as.data.frame()
tpms11A_log2 <- tpms11A_log2 %>% t() %>% as.data.frame()
rownames(tpms01A_log2) <- substring(rownames(tpms01A_log2),1,12)
rownames(tpms11A_log2) <- substring(rownames(tpms11A_log2),1,12)
a <- intersect(rownames(tpms01A_log2),rownames(tpms11A_log2))
tpms01A_log2 <- tpms01A_log2[a,]
tpms11A_log2 <- tpms11A_log2[a,]
peidui <- cbind(tpms11A_log2[,gene],tpms01A_log2[,gene])#11A放在前面
peidui <- as.data.frame(peidui)
write.csv(peidui,file = "peidui.csv")


####根据CXCL10高低组做生存分析####
setwd("TCGA-PAAD")
setwd("survival")
surv <- read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv$OS.time <- surv$OS.time/365

#median 中位数
#CXCL10
surv$CXCL10
surv$group <- ifelse(surv$CXCL10 > median(surv$CXCL10),"High","Low")
class(surv$group)
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 拟合生存曲线
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, # 显示置信区间
           risk.table = TRUE, # 显示风险表
           risk.table.col = "strata",
           palette = "jco", # 配色采用jco
           legend.labs = c("Low", "High"), # 图例
           size = 1,
           xlim = c(0,20), # x轴长度
           break.time.by = 5, # x轴步长为5
           legend.title = "CXCL10",
           surv.median.line = "hv", # 限制垂直和水平的中位生存
           ylab = "Survival probability (%)", # 修改y轴标签
           xlab = "Time (Years)", # 修改x轴标签
           ncensor.plot = TRUE, # 显示删失图块
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()

####不同分期BTK的表达####
setwd("TCGA-PAAD")
setwd("CXCL10")
library(tidyverse)
clinical.expr01A = read.table("clinical.expr01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- "CXCL10"
clinical_CXCL10 <- cbind(clinical.expr01A[,1:6],clinical.expr01A[,gene])
write.csv(clinical_CXCL10, file = "clinical_CXCL10.csv")


####CXCL10差异分析####
setwd("TCGA-PAAD")
setwd("CXCL10_DEG")
library(DESeq2)
library(tidyverse)
counts_01A <- read.table("counts01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("tpms01A_log2.txt", sep = "\t",row.names = 1,check.names = F,header = T)
identical(colnames(counts_01A),colnames(exp))#习惯性判断以防万一
gene <- "CXCL10"#每次运行只改这个基因名
med=median(as.numeric(exp[gene,]))

conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A,
  colData = conditions,
  design = ~ group)

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)
save(res,file="DEG_CXCL10.Rda")