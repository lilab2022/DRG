# Fig3 PNS microglia-like cells share epigenetic profiles with CNS microglia.
library(DESeq2)
library(ComplexHeatmap)
setwd('E:/R/DRG/Fig3')



# --------------------- DEseq2 ------------------------------------
dat <- read.table('E:/R/DRG/Fig4/Atac_merged.counts', header = T, row.names = 1)
head(dat)
coldata <- read.csv('E:/R/DRG/Fig4/atac_info.csv', row.names = 1)
coldata <- coldata[, c(rownames(coldata))]
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = coldata,
                              design= ~ sample)
dds <- DESeq(dds)
t <- 'DRG11'
c <- 'Brain11'
pcw <- '11'
res <- results(dds, contrast = c('condition', t, c))
res_order <- res[order(res$padj), ]
summary(res)
res_filter <- as.data.frame(res)
res_filter$chr_info <- rownames(res_filter)
#res_filter <- filter(res_filter, res_filter$padj < 0.05)
res_filter <- na.omit(res_filter)
FC <- 1.5
FDR <- 0.05
res_filter$Significant = 'Normal'

up <- intersect(which(res_filter$log2FoldChange > log2(FC)),
                which(res_filter$padj < FDR) )

down <- intersect(which(res_filter$log2FoldChange< (-log2(FC))),
                  which(res_filter$padj < FDR))
res_filter$Significant[up] <- paste(t, 'Up', sep = ' ')
res_filter$Significant[down] <- paste(c, 'Up', sep = ' ')
table(res_filter$Significant)

write.csv(res_filter, paste('deseq2',pcw, t, c, '.csv', sep = '_'))


## Fig4B 
# --------------------------------------------cor plot----------------------------------------------------------
res_cor <- results(dds, contrast = c('condition', t, c))
res2 = as.data.frame(res_cor)
res2$abs_fc <- abs(res2$log2FoldChange)
res2 <- res2[order(-res2[, 'abs_fc']), ]
now2 <- rownames(res2)[1:1000]
dat_now2 <- assay(vsd)[now2,]
dat_now_centered2 <- t(apply(dat_now2, 1, function(x) x-mean(x)))
r = cor(dat_now_centered2, method = 'pearson')
r
pdf('./cor_atac.pdf', width = 6, height = 6)
pheatmap(r, 
         show_colnames = TRUE,  
         show_rownames=TRUE,    
         fontsize=10,             
         color = colorRampPalette(c('navyblue','#ffffff',  'red'))(101),
         annotation_legend=TRUE, 
         border_color=NA,        
         scale="none",      
         breaks=unique(c(seq(-0.8,.8, length=101))),
         main = 'atac')
dev.off()



                             
## Fig4C
# --------------------------------------------- Heatmap -------------------------------------------------------------------
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  repelled.y <- function(d, d.select, k = repel.degree){
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  grid.newpage()
  grid.draw(heatmap)
  invisible(heatmap)
}

chr_label <- c('10_122461037_122462265', #HTRA1
               '16_51151515_51152835',  #SALL1
               '12_108597927_108598397', #TMEM119
               '19_6693636_6694638', '19_6706562_6706965', '19_6718760_6719380',#C3,
               '3_151387641_151388022',#P2RY12
               '10_17808928_17809837',#MRC1
               '5_39424741_39425482', #DAB2
               '5_140632539_140633511',#CD14
               '11_60280361_60280873')


bp <- read.csv('deseq2_11_Brain11_Pia11.csv')
bp <- filter(bp, abs(bp$log2FoldChange) >= 0.5 & bp$padj < 0.05)
bp$group <- 'bp'
dp <- read.csv('deseq2_11_DRG11_Pia11.csv')
dp <- filter(dp, abs(dp$log2FoldChange) >= 0.5 & dp$padj < 0.05)
dp$group <- 'dp'
dat <- rbind(bp, dp)
table(dat$Significant)
chrs <- unique(dat$X)


col_sums <- colSums(dat) # mer_cou是1中读入的counts矩阵
expnor <- sweep(dat*10000, 2, col_sums, `/`) 
expnor <- log2(expnor + 1)
exp0 <- as.data.frame(t(expnor[chrs, ]))
exp_t <- as.data.frame(t(exp0))
ex_H <- exp_t[grep(c("^H"), rownames(exp_t)), ] 
ex_K <- exp_t[grep(c("^K"), rownames(exp_t)), ] 
ex_G <- exp_t[grep(c("^G"), rownames(exp_t)), ] 
ex_MT <- exp_t[grep(c("^MT"), rownames(exp_t)), ] 
ex <- rbind(ex_H, ex_K, ex_G, ex_MT)
save_chr <- setdiff(colnames(exp0), rownames(ex))
expmean_final <- exp0[, save_chr]



                             
pdf('./atac_heatmap.pdf', width = 15, height = 40)
p1 <- pheatmap(as.matrix(t(expmean_final)),
               clustering_method = 'ward.D2',
               clustering_distance_cols = 'maximum',
               scale = 'row',
               color =colorRampPalette((c('navyblue','white',  '#A00021')))(300),
               fontsize_row = 8,
               cellwidth = 25, 
               cellheight = 2,
               show_colnames = T,
               cluster_rows = T,
               cluster_cols = T,
               #gaps_col = c(11,20),
               angle_col = 45,
               #annotation_row = annotation_cols,
               breaks=unique(c(seq(-1,1, length=300))))
add.flag(p1,
         kept.labels = chr_label,
         repel.degree = 0.2)
dev.off()
