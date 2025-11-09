####COX####
setwd("TCGA-PAAD")
setwd("COX")
#安装加载R包
#install.packages("survival")
#install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#读取文件
exp_surv_01A = read.table("exp_surv_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#手动读取DEG_final.txt 
#提取DEG_final
surv.expr <- cbind(exp_surv_01A[,1:2],exp_surv_01A[,DEG_final$SYMBOL])
#a <- exp_surv_01A[,1:2]
#b <- exp_surv_01A[,DEG_final$SYMBOL]
#Cox分析
#如何修改特定列的列名
#colnames(surv.expr)[ ] <- ""  #[]内填特定列数字 ""内填写修改的名字
Coxoutput <- NULL 

for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}

Coxoutput <- arrange(Coxoutput,pvalue)
###筛选top基因
gene_sig <- Coxoutput[Coxoutput$pvalue < 0.005,] # 取出p值小于0.05的基因
write.csv(gene_sig, file = "gene_sig.csv")
topgene <- gene_sig #为了下面不改topgene
#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "53" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
#保存图片 大小30*30
dev.off()
##cox&ppi韦恩图取交集：CXCL10 CXCL9