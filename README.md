# cervical-cancer

# Getting started

## Step 1. Pre run the method for installation

You should ensure that you have the necessary system dependencies configured.

For Windows (8.1 / 10 / 11): Rtools should be installed to the system path.

The latest base R is recommended. The compatibility of the earlier version (v4.0.x) is under evaluation.
We use R version is [64-bit] d:\Program Files\R\R-4.2.2

## Step 2. Install the package
The dependency `GEOquery`, `limma`, `SPIA` and `ggplot2` are unavailable on the CRAN but available on [BioConductor](https://www.bioconductor.org/). So we need to install the BiocManager manually. 

``` r
if (!"BiocManager" %in% as.data.frame(installed.packages())$Package)
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma","SPIA","SPIA","ggplot2"))

## Examples

 Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
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
