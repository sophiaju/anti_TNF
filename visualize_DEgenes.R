# clear environment
rm(list=ls())
# CHANGE This to directory downloaded to
PARENTDIR <- "C:/Users/sophi/Documents/McGill/1_Academics/prof_li_research/DE_genes"
setwd(PARENTDIR)

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

mat[c('CENPK', 'HLA-B', 'KLRD1', 'IFNG', 'CCL4'), ]
