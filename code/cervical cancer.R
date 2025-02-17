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
#differentially expressed genes of CC
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

#27 disulfidptosis-related genes
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










