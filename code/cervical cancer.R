# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO
gset <- getGEO("GSE63514", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("111111111111111111111111XXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "0000000000000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("tumor","normal"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=25000000)

tT1 <- subset(tT, select=c("ID","Gene.ID","adj.P.Val","P.Value","t","B","logFC"
                           ,"Gene.symbol","Gene.title"))

tT1 <- tT1[!duplicated(tT1$Gene.ID),]


tT1 <- tT1[-1,]


genes <- c()
symbols <- c()
FC <- c()
pvalue <- c()
probe <- c()
for(i in 1:length(tT1$Gene.ID)){
  m <- strsplit(tT1[i,]$Gene.ID,'///')[[1]]
  m1 <- strsplit(tT1[i,]$Gene.symbol,'///')[[1]]
  
  
  FC1 <- rep(tT1[i,]$logFC,length(m))
  pvalue1 <- rep(tT1[i,]$adj.P.Val, length(m))
  m2 <- rep(tT1$ID[i],length(m))
  
  genes <- c(genes, m)
  symbols <- c(symbols, m1)
  FC <- c(FC, FC1)
  pvalue <- c(pvalue, pvalue1)
  probe <- c(probe, m2)
}


tg_all <- data.frame(probe=probe,genes = genes,symbols=symbols, FC = FC, pvalue = pvalue)

tg_all1 <- tg_all[!duplicated(tg_all$genes),]

deg <- tg_all1[tg_all1$pvalue<0.05,]
#deg1 <- deg[abs(deg$FC)>1,]



disf <- c("SLC7A11", "SLC3A2", "RPN1", "NCKAP1", "NUBPL", "NDUFA11", "LRPPRC", "OXSM", "NDUFS1", 
          " GYS1", "MYH9", "MYH10", "MYL6", "PRDX1", 
          "RAC1", "DSTN", "IQGAP1", "CD2AP", "ACTN4",
          "PDLIM1", "FLNB", "ACTB", "INF2", "WASF2", "CYFIP1", "ABI2", "BRK1")

deg1 <- deg

m <- intersect(deg1$symbols, disf)


deg_F <- deg1[match(m,deg1$symbols),]
deg2 <- deg1[match(m,deg1$symbols),]
ex1 <- exprs(gset)
ex_f <- ex1[tg_all1$probe,]
rownames(ex_f) <- tg_all1$symbols
ex_f <- ex1[deg2$probe,]
rownames(ex_f) <- deg2$symbols
ex1 <- as.data.frame(ex1)




#volcano plot of degs

#tg_all1 <- tT1

library(dplyr)
library(ggplot2)
dat  = tg_all1[match(disf, tg_all1$symbols),]

dat <- tg_all1

dat <- dat[!is.na(dat$symbols),]


dat$change=ifelse(dat$pvalue>0.05,'stable', #if 判断：如果这一基因的P.Value>0.01，则为stable基因
                  ifelse( dat$FC >0,'up', #接上句else 否则：接下来开始判断那些P.Value<0.01的基因，再if 判断：如果logFC >1.5,则为up（上调）基因
                          ifelse( dat$FC < -0,'down','stable') )#接上句else 否则：接下来开始判断那些logFC <1.5 的基因，再if 判断：如果logFC <1.5，则为down（下调）基因，否则为stable基因
)

p <- ggplot(data = dat,
            aes(x = FC,y = -log10(pvalue))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-0,0),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  theme_bw()



#heatmap
n=ex_f
library(pheatmap)
annotation_col=data.frame(group=gs1)
rownames(annotation_col) <- colnames(n)
n2 = t(scale(t(n)))
n=t(scale(t(ex_f))) # 'scale'可以对log-ratio数值进行归一???
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(group=gs1)
rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息

pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac)

ex_f1 <- t(ex_f)
write.csv(ex_f1,'ex_f1_all.csv')









library(bseqsc) #携带大量CIBERSORT的依赖
library(CIBERSORT)

library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

data(LM22) # 自带LM22文件
data(LM22)
##       B cells naive B cells memory Plasma cells T cells CD8 T cells CD4 naive
## ABCB4        555.71          10.74        7.226       4.311             4.606
## ABCB9         15.60          22.09      653.392      24.224            35.672
## ACAP1        215.31         321.62       38.617    1055.613          1790.097
## ACHE          15.12          16.65       22.124      13.428            27.188
## ACP5         605.90        1935.20     1120.105     306.313           744.657

# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr

ex1[tT1$ID,]
results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)

# 理解一下results的结果
# 你可以理解为返回了一个列名为细胞类型、行名为样本名的细胞浸润程度（占比）的矩阵
# 此外result中还会多出三列：
# P-value: 用来展示去卷积的结果在所有细胞类群中是否具有差异
# Correlation:参考矩阵与输入矩阵的特征基因相关性
# RMSE: Root mean squared error，参考矩阵与输入矩阵的特征基因标准差

# heatmap
# 按行（样本内部）标准化可以看出在各类样本内部，M2浸润程度（占比）最高
rowscale <- results[,1:ncol(LM22)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=F,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度



# perm置换次数 = 1000，QN分位数归一化 = TRUE 
FPKM <- ex1[tg_all1$probe,]
rownames(FPKM) <- tg_all1$symbols

FPKM1 <- as.matrix(FPKM)

results <- cibersort(sig_matrix = LM22, mixture_file = FPKM1)

results[1:6,1:5]


library(reshape)

# 构造注释文件
group_list <- ifelse(as.numeric(substring(rownames(results), 14, 15)) < 10, "Tumor","Normal") %>% factor(., levels = c("Normal", "Tumor"))

# 取前22列细胞亚群组
LUAD_data <- as.data.frame(results[, 1:22])
group_list <- gset$group %>% factor(., levels = c("normal", "tumor"))
LUAD_data$group <- group_list
LUAD_data$sample <- row.names(results)

# 融合数据
LUAD_New = melt(LUAD_data)
colnames(LUAD_New) = c("Group", "Sample", "Celltype", "Composition")  #设置行名
head(LUAD_New)
##   Group                       Sample      Celltype Composition
##  1 Tumor TCGA-35-5375-01A-01R-1628-07 B cells naive  0.00000000
##  2 Tumor TCGA-55-A4DF-01A-11R-A24H-07 B cells naive  0.06919370
##  3 Tumor TCGA-95-8039-01A-11R-2241-07 B cells naive  0.06389702
##  4 Tumor TCGA-MP-A4T4-01A-11R-A262-07 B cells naive  0.03455365
##  5 Tumor TCGA-62-A471-01A-12R-A24H-07 B cells naive  0.13974877
##  6 Tumor TCGA-L9-A5IP-01A-21R-A39D-07 B cells naive  0.11928189

library(ggpubr)

ggplot(LUAD_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell Proportion", x =  NULL, title = "Cervel Cancer Cell Proportion")+  
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#096EA9", "#B33D27")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)



genes <- m

genes_expr <- as.data.frame(t(FPKM1[rownames(FPKM1) %in% genes,]))
genes_expr <- genes_expr[match(rownames(LUAD_data),rownames(genes_expr)),]
identical(rownames(LUAD_data),rownames(genes_expr))



library(psych)
x <- genes_expr
y <- immnune_data
d <- corr.test(x,y,use="complete",method = 'spearman')

r <- d$r
p <- d$p

library(ggcorrplot)
ggcorrplot(t(d$r), show.legend = T, 
           p.mat = t(d$p.adj), digits = 2,  sig.level = 0.05,insig = 'blank',lab = T)

library(linkET)
immnune_data <- LUAD_data[,-24]
immnune_data <- immnune_data[,-23]

cor_res <- correlate(genes_expr, immnune_data,method = "spearman")
qcorrplot(cor_res) +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

library(tidyr)
library(tibble)
# 先整理下数据
df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "correlation")

df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`

library(ggplot2)

ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())


ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='blue', high='red',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())



#Lasso

library(glmnet)
df <- ex_f1

X <- as.matrix(df)

Y <- as.matrix(sml)

lasso_model <- glmnet(X, Y, family = 'binomial', alpha = 1)

seleced_gene <- which(coef(lasso_model)!=0)

plot(lasso_model)

set.seed(123)
cvfit<-cv.glmnet(X,Y, family = "binomial",alpha=1,type.measure = 'deviance')
print(cvfit)
# Call:  cv.glmnet(x = x, y = y, type.measure = "deviance", family = "binomial",      alpha = 1) 
# 
# Measure: Binomial Deviance 
# 
#      Lambda Index Measure      SE Nonzero
# min 0.02141    27  0.8091 0.08923      19
# 1se 0.04945    18  0.8944 0.06460      13
plot(cvfit)

cvfit$lambda.min
#[1] 0.02140756
cvfit$lambda.1se
#[1] 0.04945423
coef<-coef(cvfit, s = "lambda.1se")
# 31 x 1 sparse Matrix of class "dgCMatrix"
#                      s1
# (Intercept)  0.20679470
# V1           .         
# V2           0.27873628
# V3          -0.16579666
# V4          -0.69275750
# V5          -0.07549238
# V6          -0.42352625
# ...

coef<-coef@Dimnames[[1]][which(!coef==0)]
coef[-1]
#[1] "V2"  "V3"  "V4"  "V5"  "V6"  "V8"  "V9"  "V10" "V22" "V23" "V25" "V26" "V29"
length(coef[-1])



#SVM

library(tidyverse)
library(caret)
library(ggplot2)
library(cowplot)
library(ggplotify)

Y<- as.numeric(sml)

set.seed(21) # 设置种子
control <- rfeControl(functions = caretFuncs, method = "cv", number = 10) # cv 交叉验证次数10
# 执行SVM-RFE算法
num <- 11
results <- rfe(x = X, # 除去最后一列，其余列均为预测变量（也就是hubgene的表达量）
               y = Y, # 分组信息
               sizes = c(1:num), 
               rfeControl = control,
               method = "svmRadial"
)

## 结果分析
svmrfe_result <- data.frame(symbol = predictors(results)) ## 7个基因

# SVM-RFE结果简单可视化
p1 <- plot(results, type=c("o"),
           xgap.axis = 1)
p1 <- as.ggplot(plot_grid(p1))+
  labs(title="SVM_RFE_analyse", x="", y = "",size=25) +
  # theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=25),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p1


#randomforest

library(randomForest)

randomForest_select = function(gene_exp,sample_info,tree_num,rate){
  #data = t(gene_exp)
  data = bind_cols(sample_info,gene_exp)
  
  rf_model <- randomForest(data[,-1], data[,1], ntree = tree_num)
  # 获取特征的重要性分数
  importance_scores <- importance(rf_model)
  
  # 将特征名称和对应的重要性分数合并
  feature_scores <- data.frame(Features = row.names(importance_scores), Importance = importance_scores[, 1])
  
  # 按重要性分数降序排序
  sorted_scores <- feature_scores[order(-feature_scores$Importance), ]
  
  number = round(nrow(sorted_scores)*rate,digits = 0)
  
  result_rf = sorted_scores[1:number,]
  
  return(result_rf)
}

library(ggsci)

sample_info <- as.data.frame(Y)
randomforest_result = randomForest_select(X,sample_info,5000,0.6)

ggplot(randomforest_result,aes(x = Features %>% sort(), y = Importance))+
  geom_col(aes(fill = Importance),#注意条形图的颜色使用fill
           width = 0.7)+ 
  scale_fill_gsea()+  #注意这里要用fill而不是color
  theme_bw() + theme(axis.text = element_text(hjust = 1,
                                              angle = 60), legend.position = c(0.9,
                                                                               
                                                                               
                                                                               0.9), legend.direction = "horizontal") +labs(x = NULL)
data = bind_cols(sample_info,X)

rf_model <- randomForest(data[,-1], data[,1], ntree = 5000)
varImpPlot(rf_model)



#

library(rms)

rt <- data[,c('Y',"NDUFA11", "NDUFS1",  "RAC1" , "BRK1")]



rt$NDUFA11 <- ifelse(rt$NDUFA11 >= mean(rt$NDUFA11),1,0)
rt$NDUFA11 <- factor(rt$NDUFA11,labels = c(0,1))


rt$NDUFS1 <- ifelse(rt$NDUFS1 >= mean(rt$NDUFS1),1,0)
rt$NDUFS1 <- factor(rt$NDUFS1,labels = c(0,1))

rt$RAC1 <- ifelse(rt$RAC1 >= mean(rt$RAC1),1,0)
rt$RAC1 <- factor(rt$RAC1,labels = c(0,1))

rt$BRK1 <- ifelse(rt$BRK1 >= mean(rt$BRK1),1,0)
rt$BRK1 <- factor(rt$BRK1,labels = c(0,1))


ddist<-datadist(rt)
options(datadist="ddist")

fit <- lrm(Y~ NDUFA11+ NDUFS1 + RAC1 + BRK1, data=rt, x= T, y = T)

nom<- nomogram(fit, 
               lp.at=seq(-20,40,by=2),
               fun=plogis, 
               fun.at=c(0.1,0.3,0.5,0.7,0.9),
               lp=F,
               funlabel="Risk"
               
)   #构建Nomogram图5  
plot(nom)   #输出Nomogram图

plot(nom,
     #lplabel="linear Predictor", 
     fun.side=c(3,1,3,1,3), 
     label.every=2, 
     col.grid=gray(c(0.85,0.95))) 

nom<-nomogram(fit,
              lp.at=seq(-2,40,by=0.5), 
              fun=function(x)1/(1+exp(-x)), 
              funlabel="Risk of Death",
              fun.at=c(0.05,seq(0.1,0.9,by=0.1),0.95), 
              conf.int=c(0.1,0.3))



nom <- nomogram(fit)
library(nomogramFormula)##加载nomogramFormula包
results<-formula_rd(nomogram=nom)
rt$points<-points_cal(formula = results$formula,rd=rt)##生成每个个体分数
pre<-rt$points

library(pROC)##加载pROC包
plot.roc(rt$Y, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, of="thresholds",
         thresholds="best",
         print.thres="best")##构建roc曲线
rocplot1 <- roc(rt$Y,pre)
ci.auc(rocplot1)##计算ROC下面积AUC区间

plot(rocplot1, 
     main = paste0("GSE63514 "), # 设置主标题
     col = "red", # 设置曲线颜色
     lwd = 2,  # 线条宽度
     print.auc = T, # 打印AUC的值
     print.auc.x = 0.4, # AUC值得位置（x轴）
     print.auc.y = 0.5, # AUC值得位置（y轴）
     print.auc.pattern = 'AUC of Nomogram = %.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     grid.col = "gray", # 网格线颜色
     xlab = "False Positive Rate (1 - Specificity)",  # x轴标签  
     ylab = "True Positive Rate (Sensitivity)",  # y轴标签 
     font.lab = 2,  # 轴标签字体样式  
     cex.axis = 1.2,  # 轴标签字体大小  
     cex.main = 1.5  # 主标题字体大小  
)




## 绘制校正曲线

## 重新调整模型函数f2，也即添加x=T, y=T
f2 <- lrm(Y~ NDUFA11+ NDUFS1 + RAC1 + BRK1, data=rt, x=T, y=T) 
## 构建校正曲线



cal2 <- calibrate(f2, B=1000)
plot(cal2,lwd=2,lty=1,
     conf.int=T,# 是否显示置信区间
     errbar.col="blue",#直线曲线bar颜色
     col="red", # 曲线颜色
     xlim=c(0,1),ylim=c(0,1),
     xlab="Predicted Probability",
     ylab="Observed Probability",
     subtitles = F)




#单基因分类
library(glmnet)
set.seed(123)
data <- genes_expr
data$Label <- LUAD_data$group


library(pROC)

roc <- roc(response = data$Label, 
           predictor = data[,1],
           levels = c('normal', 'tumor'))


roc2 <- roc(response = data$Label, 
            predictor = data[,2],
            levels = c('normal', 'tumor'))

roc3 <- roc(response = data$Label, 
            predictor = data[,8],
            levels = c('normal', 'tumor'))

roc4 <- roc(response = data$Label, 
            predictor = data[,9],
            levels = c('normal', 'tumor'))

plot(roc, 
     main = paste0("ROC Curve for Disease Prediction of ", colnames(data)[1]), # 设置主标题
     col = "#DD7123", # 设置曲线颜色
     lwd = 2,  # 线条宽度
     print.auc = T, # 打印AUC的值
     print.auc.x = 0.4, # AUC值得位置（x轴）
     print.auc.y = 0.5, # AUC值得位置（y轴）
     print.auc.pattern = 'NDUFS1 AUC=%.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     grid.col = "gray", # 网格线颜色
     xlab = "False Positive Rate (1 - Specificity)",  # x轴标签  
     ylab = "True Positive Rate (Sensitivity)",  # y轴标签 
     font.lab = 2,  # 轴标签字体样式  
     cex.axis = 1.2,  # 轴标签字体大小  
     cex.main = 1.5  # 主标题字体大小  
)

plot(roc2,add=TRUE, 
     print.auc = T,
     print.auc.x = 0.4, # AUC值得位置（x轴）
     print.auc.y = 0.45, # AUC值得位置（y轴）
     print.auc.pattern = 'NDUFA11 AUC=%.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     col="blue")

plot(roc3,add=TRUE, 
     print.auc = T,
     print.auc.x = 0.4, # AUC值得位置（x轴）
     print.auc.y = 0.4, # AUC值得位置（y轴）
     print.auc.pattern = 'BRK1 AUC=%.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     col="black")

plot(roc4,add=TRUE, 
     print.auc = T,
     print.auc.x = 0.4, # AUC值得位置（x轴）
     print.auc.y = 0.35, # AUC值得位置（y轴）
     print.auc.pattern = 'RAC1 AUC=%.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     col="green")




data <- t(ex_f)
data <- as.data.frame(data)
rt <- data[,c("NDUFA11", "NDUFS1",  "RAC1" , "BRK1")]
rt$Y <- gs


ddist<-datadist(rt)
options(datadist="ddist")

fit <- lrm(Y~ NDUFA11+ NDUFS1 + RAC1 + BRK1, data=rt, x= T, y = T)


nom <- nomogram(fit)
library(nomogramFormula)##加载nomogramFormula包
results<-formula_rd(nomogram=nom)
rt$points<-points_cal(formula = results$formula,rd=rt)##生成每个个体分数
pre<-rt$points

library(rms)
library(pROC)##加载pROC包
plot.roc(rt$Y, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, of="thresholds",
         thresholds="best",
         print.thres="best")##构建roc曲线
rocplot1 <- roc(rt$Y,pre)
ci.auc(rocplot1)##计算ROC下面积AUC区间


plot(rocplot1, 
     main = paste0("GSE67522 "), # 设置主标题
     col = "red", # 设置曲线颜色
     lwd = 2,  # 线条宽度
     print.auc = T, # 打印AUC的值
     print.auc.x = 0.6, # AUC值得位置（x轴）
     print.auc.y = 0.2, # AUC值得位置（y轴）
     print.auc.pattern = 'AUC of Nomogram = %.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     grid.col = "gray", # 网格线颜色
     xlab = "False Positive Rate (1 - Specificity)",  # x轴标签  
     ylab = "True Positive Rate (Sensitivity)",  # y轴标签 
     font.lab = 2,  # 轴标签字体样式  
     cex.axis = 1.2,  # 轴标签字体大小  
     cex.main = 1.5  # 主标题字体大小  
)



#
gene <- tT[match(c("NDUFA11", "NDUFS1",  "RAC1" , "BRK1"),tT$Gene.symbol),]
df <- ex1[match(gene$ID,rownames(ex1)),]
rownames(df) <- gene$Gene.symbol
rt <- t(df)
rt <- as.data.frame(rt)
rt$Y <- gs


ddist<-datadist(rt)
options(datadist="ddist")

fit <- lrm(Y~ NDUFA11+ NDUFS1 + RAC1 + BRK1, data=rt, x= T, y = T)


nom <- nomogram(fit)
library(nomogramFormula)##加载nomogramFormula包
results<-formula_rd(nomogram=nom)
rt$points<-points_cal(formula = results$formula,rd=rt)##生成每个个体分数
pre<-rt$points

library(rms)
library(pROC)##加载pROC包
plot.roc(rt$Y, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, of="thresholds",
         thresholds="best",
         print.thres="best")##构建roc曲线
rocplot1 <- roc(rt$Y,pre)
ci.auc(rocplot1)##计算ROC下面积AUC区间


plot(rocplot1, 
     main = paste0("GSE52903"), # 设置主标题
     col = "red", # 设置曲线颜色
     lwd = 2,  # 线条宽度
     print.auc = T, # 打印AUC的值
     print.auc.x = 0.6, # AUC值得位置（x轴）
     print.auc.y = 0.2, # AUC值得位置（y轴）
     print.auc.pattern = 'AUC of Nomogram = %.3f', # AUC值得格式，表明保留三位小数
     print.auc.cex = 1.2,  # AUC值字体大小
     grid = c(0.5,0.2), # 网格线设置
     grid.col = "gray", # 网格线颜色
     xlab = "False Positive Rate (1 - Specificity)",  # x轴标签  
     ylab = "True Positive Rate (Sensitivity)",  # y轴标签 
     font.lab = 2,  # 轴标签字体样式  
     cex.axis = 1.2,  # 轴标签字体大小  
     cex.main = 1.5  # 主标题字体大小  
)



