# download necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("IHW")


# clear environment
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "C:/Users/sophi/Documents/McGill/1_Academics/prof_li_research/DE_genes"
setwd(PARENTDIR)

library("DESeq2")
library("pasilla")
library("IHW")

# Read in matrix and Sample annotation
mat <- read.csv("MH.SARDs.RNASeq.COUNT_PrePost.csv",row.names=1)
SampleTableCorrected <- read.csv("CovariableModifiedRNASeq_PrePost.csv",row.names=1)

#only keep counts for TCells and samples that are for TCells for RA patients pre-treatment
mat <- mat[,grepl("_TC_",names(mat))]
SampleTableCorrected <- SampleTableCorrected[SampleTableCorrected$cell_type=="TCELL",]
SampleTableCorrected <- SampleTableCorrected[SampleTableCorrected$Disease=="RA",]
SampleTableCorrected <- SampleTableCorrected[SampleTableCorrected$Timing=="Pre-treatment",]

# Find the overlapping samples, 
inter <- intersect(colnames(mat),rownames(SampleTableCorrected))
SampTabC <- SampleTableCorrected[inter,]
mat <- mat[,match(rownames(SampleTableCorrected),colnames(mat))]
#check if the samples match up
all(rownames(SampTabC) == colnames(mat))

# Create object
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = SampTabC,
                              design = ~ Sex + R_vs_NR_comb)

# prefiltering, keep rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#drop levels which do not have samples in the current dds
# (unnesessary for this dataset)
dds$R_vs_NR_comb <- droplevels(dds$R_vs_NR_comb)

# Reset the levels of colData factor variables
colData(dds) <- DataFrame(droplevels(data.frame(colData(dds))))

dds$R_vs_NR_comb <- relevel(dds$R_vs_NR_comb, ref = "Yes")

dds <- DESeq(dds)
res <- results(dds, filterFun=ihw, alpha=.01)

#plot heatmap of count matrix
library("pheatmap")
ntd <- normTransform(dds) # normalize, log2(n+1) i think
vsd <- vst(dds, blind=TRUE) # variance stabilizing transformation
rld <- rlog(dds, blind=TRUE) # regularized log transformation


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[c('CENPK', 'HLA-B', 'KLRD1', 'IFNG', 'CCL4')]
df <- as.data.frame(colData(dds)[,c("R_vs_NR_comb", "Sex")])
pheatmap(assay(vsd)[c('CENPK', 'HLA-B', 'KLRD1', 'IFNG', 'CCL4'),], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

exp_info <- assay(vsd)[c('CENPK', 'HLA-B', 'KLRD1', 'IFNG', 'CCL4'),]

# plotPCA(vsd, intgroup=c("R_vs_NR_comb", "Sex"))

resOrderedpPval <- res[order(res$pvalue),]
resOrderedPadj <- res[order(res$padj),]

a <- rowMeans(counts(dds,normalized=TRUE))
b <- order(a[c('CENPK', 'HLA-B', 'KLRD1', 'IFNG', 'CCL4')])

summary(res)

#How many p-values were less than 0.01?
sum(res$pvalue < 0.01, na.rm=TRUE)

#How many adjusted p-values were less than 0.01?
sum(res$padj < 0.01, na.rm=TRUE)

write.table(resOrderedpPval,'IHW_results/DESEQ2_TCELL_RvsNR_pvalue.csv',sep=",", quote=F, col.names=NA)
write.table(resOrderedPadj,'IHW_results/DESEQ2_TCELL_RvsNR_padj.csv',sep=",", quote=F, col.names=NA)
