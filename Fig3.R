# Fig3 PNS microglia-like cells share epigenetic profiles with CNS microglia.
library(DESeq2)
library(ComplexHeatmap)

setwd('/Volumes/Nexus/Projects/Project_Mag/bam/')
dat <- read.table('H3K27ac_merged.counts', header = T, row.names = 1)
head(dat)

coldata <- data.frame(row.names = colnames(dat),
                      sample = rep(c('Brain', 'DRG', 'Pia'), c(4,2,2)))
coldata

dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = coldata,
                              design= ~ sample)
dds <- DESeq(dds)
vsd <- vst(dds)
res <- results(dds)
res <- na.omit(res)
res

## Fig4F 
# --------------------------------------------Monocle3 cell order----------------------------------------------------------
