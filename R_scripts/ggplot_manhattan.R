require(stringr)
require(ggplot2)
require(dplyr)
library(ggrepel)

rm(list=ls())

PARENTDIR <- "/Users/sophi/Documents/McGill/1_Academics/prof_li_research"
setwd(PARENTDIR)

gwas <- read.table('manhattan/allchr.qassoc', header=TRUE)
#gwas <- read.table('manhattan/new_deseq/all_gene_loc_pvalue.txt', header=TRUE)
#gwas <- read.table('manhattan/new_deseq/DESEQ2_TCELL_RvsNR_padj.csv', header=TRUE, sep=',')


#osnps <- str_split(readLines('McGill/1_Academics/prof_li_research/manhattan/GWAS_snps_uniq_Mb.txt'), pattern=' ')
#osnps_c <- as.character(osnps)

#deseq_overlap <- str_split(readLines('McGill/1_Academics/prof_li_research/manhattan/DESEQ_uniq_overlap.txt'), pattern=' ')
#deseq_overlap <- as.list(str_split(readLines('DE_genes/Tao_ADA_DEgenes.txt', warn=FALSE), "[, ]+")[[1]])
#osnps_c <- as.character(deseq_overlap)

#for pvalue graph
snp_toDEgenes_overlap <- read.table('DE_genes/rsid_to_DEgenes.csv', sep = ',')
snpsOfInterest <- as.list(snp_toDEgenes_overlap[[1]])
temp <- snp_toDEgenes_overlap[,1]
rownames(snp_toDEgenes_overlap) <- temp

# for Tao ADA overlaps
# snp_toDEgenes_overlap <- read.table('manhattan/multiomics/SNP_ada_overlap.csv', sep = ',')
# snpsOfInterest <- as.list(snp_toDEgenes_overlap[[1]])
# temp <- snp_toDEgenes_overlap[,1]
# rownames(snp_toDEgenes_overlap) <- temp

#for padj graph
# snp_toDEgenes_overlap_padj <- read.table('manhattan/new_deseq/overlaps_padj01_Mb.csv', sep=',', header=TRUE)
# snpsOfInterest <- as.list(snp_toDEgenes_overlap_padj[["SNP_id"]])
# temp <- snp_toDEgenes_overlap_padj[["SNP_id"]]
# rownames(snp_toDEgenes_overlap_padj) <- temp


#snpsOfInterest <- osnps_c
gwas <- na.omit(gwas) 
gwas %>% filter(-log10(P)>1)

don <- gwas %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric((chr_len)))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate(is_annotate=ifelse(SNP %in% snpsOfInterest, "yes", "no"))%>%
  # add DE gene label information for labeling
  #pvalue
  mutate(DEgenes=ifelse(SNP %in% snpsOfInterest, snp_toDEgenes_overlap[SNP,2], 'no'))
  #padj
  #mutate(DEgenes=ifelse(SNP %in% snpsOfInterest, snp_toDEgenes_overlap_padj[SNP,3], 'no'))

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  xlab("Chromosome") +
  
  #set y axis
  #ylim(0,10) +
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=DEgenes), size=2, max.overlaps=30) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print('done')
