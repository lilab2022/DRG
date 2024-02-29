#Fig1 PNS microglia are present in human, macaque and pig but absent from rodent models.
library(SingleR)
library(celldex)
library(pheatmap)
library(ggplot2)
library(Seurat)
library(dplyr)
library(reshape2)
library(symphony)
library(Matrix)
library(cowplot)
library(foreach)
library(metacell)
library(tgconfig)
library(chameleon)
library(homologene)
library(SeuratWrappers)
library(biomaRt)
library(stringr)
library(corrplot)
library(DropletUtils)
library(PerformanceAnalytics)
library(tinyarray) 
setwd('E:/R/DRG/Fig1')

## Fig1b-----------------------2D-projection of adult immune cells of different species----------------------------------

 # After Metacell
hmeta <- read.csv('E:/R/DRG/Fig1/h83.csv', row.names = 1)
drg_theme <- function(){
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = 'black'),
        plot.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, color = "black"), 
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 22, color = "black"),
        axis.text = element_text(size = 15, color = "black"), 
        axis.title = element_text(size = 15, color = "black"))
        }

p1 <- ggplot(hmeta, aes(x = sc_x, y = sc_y, color = anno)) + 
     geom_point(shape = 19, size = 2, stroke = 0.1) +
     scale_color_manual(values = c('#B5B5B6', 'royalblue', 'hotpink1')) +
     drg_theme()
p21 <- ggplot(hmeta)+
       geom_point(aes(x = sc_x, y = sc_y),shape = 19,size = 1, color = "grey88") +
       geom_point(data = hmeta[hmeta$Tissue_group %in% c('CNS'),  ], aes(x=sc_x, y=sc_y), 
       color ="#269ED1", shape = 19, size = 1)+
       labs(title = 'CNS') +
       drg_theme()
p22 <- ggplot(hmeta)+
       geom_point(aes(x = sc_x, y = sc_y),shape = 19,size = 1, color = "grey88") +
       geom_point(data = hmeta[hmeta$Tissue_group %in% c('PNS'),  ], aes(x=sc_x, y=sc_y), 
       color ="#D5221E", shape = 19, size = 1)+
       labs(title = 'PNS') +
       drg_theme()       
p <- plot_grid(p1, p21, p22)
p
ggsave("./fig1_2d_human.pdf", p, width = 6, height = 6)




 ## FigS1b--------------------------------------- Dot plot---------------------------------------------------------------
markers <- c('C3', 'P2RY12', 'TMEM119', 'SALL1', 'MRC1', 'LYVE1', 'DAB2',
             'S100A8', 'S100A9', 'LYZ', 'CLEC10A', 'CD1C', 'CD1E', 'CPA3', 'LMO4', 'TPSAB1', 
             'GNLY', 'FGFBP2', 'NCAM1', 'TRAC', 'CD3G', 'CD3D', 'CD79A', 'MS4A1', 'IGHM', 'IGLC2', 'IGKC')
order_ct <- rev(c('Microglia/Microglia-like cell', 'Mac',
                 'Mono.', 'DC', 'Gr.', 'NK', 'T', 'B'))
# using metacell to get the annotation and the mcid of each single cell
mc <- scdb_mc(matid) 
mcmc <- as.data.frame(mc@mc)
colnames(mcmc) <- 'mcid'
mcmc$cellid <- rownames(mcmc)
a <- as.data.frame(mc@annots)
colnames(a) <- 'anno'
a$mcid <- rownames(a)
mcmc <- merge(mcmc, a, by = 'mcid')
head(mcmc)   

# using Seurat to get the normalized matrix of the marker genes
mat <- scdb_mat(matid)
drg <- CreateSeuratObject(counts = mat@mat[, mcmc$cellid])
huamn_seu <- NormalizeData(huamn_seu)
dat0 <- huamn_seu@assays$RNA@data
dat1 <- dat0[markers, ]
gene_exp <- as.data.frame(t(as.matrix(dat1)))
gene_exp$cellid <- rownames(gene_exp)

# the exp
tmp <- mcmc[, c("cellid", "anno")]
colnames(tmp) <- c("cellid", "subtype")
gene_exp <- merge(gene_exp, tmp, by = "cellid", all.x = T)
gene_exp <- gene_exp[,-1]
gene_exp <- melt(gene_exp)
gene_exp2 <- gene_exp[, c("subtype","variable")]
gene_exp2 <- gene_exp2[!duplicated(gene_exp2), ]
gene_exp2$mean <- 100
gene_exp2$ratio <- 100
gene_exp2 <- na.omit(gene_exp2)
gene_exp[is.na(gene_exp)] <- 0
gene_exp2 <-  gene_exp2[gene_exp2$variable %in% markers, ]
for (i in unique(gene_exp2$subtype)){
  print(i)
  for (j in unique(gene_exp2$variable[gene_exp2$subtype == i])) {
    gene_exp2$mean[gene_exp2$subtype == i & gene_exp2$variable == j] = mean(gene_exp$value[gene_exp$subtype == i & gene_exp$variable == j])
    gene_exp2$ratio[gene_exp2$subtype == i & gene_exp2$variable == j] = nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j & gene_exp$value > 0,])/nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j,])
  }
}

# the scaled exp and the ratio
exp2 <- spread(data = gene_exp2[, c('subtype', 'variable', 'mean')],
              key = subtype,
              value = mean)
exp3 <- as.data.frame(t(exp2))
colnames(exp3) <- exp3[1, ]
exp3 <- exp3[-1, ]
names <- rownames(exp3)
exp3 <- as.data.frame(lapply(exp3, as.numeric))
rownames(exp3) <- names
exp3 <- as.data.frame(scale(exp3, center = TRUE, scale = TRUE))
exp3$subtype <- rownames(exp3)
exp4 <- melt(exp3)
exp4$ratio <- 1
for (i in unique(exp4$variable)){
  for (j in unique(exp4$subtype)) {
    exp4$ratio[exp4$subtype == j & exp4$variable == i] = gene_exp2$ratio[gene_exp2$variable == i & gene_exp2$subtype == j]
  }
}

# plot
gene_exp3 <- exp4
gene_exp3$subtype <- factor(gene_exp3$subtype, levels = rev(c(order_ct)))
gene_exp3$variable <- factor(gene_exp3$variable, levels = (c(markers)))
gene_exp3$value[gene_exp3$value >3] = 3
gene_exp3$value[gene_exp3$value < 0] = 0
p <- ggplot(gene_exp3,aes(x = variable, y= subtype, fill = value, size =ratio))+
    geom_point(shape = 21,color = "grey80" )+
    scale_fill_gradient(low = "white",high = 'darkred')+
    xlab("")+ylab("")+
    theme_bw()+
    theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,color = "black"),
          axis.text.y = element_text(color = "black"))
p
ggsave("marker_dot_fig1_human.pdf", p, width = 8, height =4)



#---------------------------------------------pie plot----------------------------------------------------------
tis <- unique(hmeta$Tissue)
for (i in tis) {
  t <- filter(hmeta, hmeta$Tissue == tis)
  t$anno <- ifelse(t$anno != 'Microglia', 'Other cells', t$anno)
  write.csv(table(t$anno), 'microglia_tissue.csv')
  ttp <- read.csv('microglia_tissue.csv', row.names = 1)
  percent <- round(100*(ttp$Freq)/sum(ttp$Freq), 2)
  percent <-paste(percent, "%", sep = "")
  pdf(paste0(i, '.pdf'), width=5, height=5)
  pie(ttp$Freq, 
      labels = percent,
      radius = 0.9 ,
      col = c('hotpink1','lightpink'),
      clockwise = TRUE)
  legend("topright", ttp$Var1, cex=0.8, 
         fill = c('hotpink1','lightpink'))
  dev.off()
}



## Fig1c ----------------------------------------------violon plot---------------------------------------------
mt <- hmeta
mt <- filter(mt, mt$anno %in% c('Mac', 'Microglia'))
mt$at <- paste(mt$anno, mt$Tissue, sep = '_')
mt2 <- filter(mt, mt$at %in% c("Microglia_Brain", "Microglia_DRG",       
                               "Mac_DRG", "Mac_ST ganglia", "Microglia_ST ganglia"))
features <- c('SALL1', 'P2RY12', 'TMEM119', 'MRC1', 'LYVE1', 'AIF1')
data <- as.data.frame(t(huamn_seu@assays$RNA@data[features, mt2$cellid]))

data$cellid <- rownames(data)
data2 <- merge(data, mt2[, c('cellid', 'at')], by = 'cellid')
rownames(data2) <- data2[, 1]
data2 <- data2[, -1]
data2 <- melt(data2)
colnames(data2) <- c('Celltype', 'Gene', 'Expr')
data2$Gene <- factor(data2$Gene, levels = features)
order_ct <- c('Microglia_DRG', 'Mac_DRG','Microglia_ST ganglia','Mac_ST ganglia', 'Microglia_Brain')
data2$Celltype <- factor(data2$Celltype, levels = rev(order_ct))

data3 <- filter(data2, data2$Gene != 'SALL1')
sall1 <- filter(data2, data2$Gene == 'SALL1' & data2$Expr >= 2)
data3 <- rbind(data3, sall1)

p1 <- ggplot(data = data3,
             aes(x = Expr, y = Celltype, fill = Celltype)) +
      geom_violin(scale = 'width',
                  draw_quantiles = c(0.25, 0.5, 0.75),
                  color = 'black',
                  linewidth = 0.3, 
                  alpha = 0.8) +
      facet_grid(cols = vars(Gene), scales = 'free_x')
p1

p2 <- p1 +
  scale_fill_manual(values = rev(c('hotpink1', 'royalblue', 'hotpink1', 'royalblue', 'hotpink1'))) +
  theme_bw() +
  labs(x = 'Log Normalized Expression')
p2 
ggsave(p2, filename = './fig1_violon.pdf',width = 6, height = 5, units = "in", device = 'pdf', dpi = 300)


##------------------------------------------------Intergrate------------------------------------------------
hg2 <- read.csv('ensemble_homology.csv', row.names = 1)

## Take monkey data as an example
monkeygene <- unique(hg2$Monkey_gene)
monmeta <- read.csv('E:/R/DRG/Fig1/Macaca/macaca_meta.csv', row.names = 1)
rownames(monmeta) <- monmeta[, 1]
monmeta <- filter(monmeta, monmeta$anno %in% c('Microglia', 'Mac'))
monkey <- Read10X('E:/R/DRG/Fig1/Macaca/Matrix')[, monmeta$cellid]
monkey <- as.data.frame(monkey[intersect(rownames(monkey), monkeygene), ])
monkey$Monkey_gene <- rownames(monkey)
mon2h <- hg2[, c('Human_gene', 'Monkey_gene')]
mon2h <- mon2h[mon2h$Monkey_gene %in% intersect(rownames(monkey), monkeygene), ]
mon2h <- mon2h[!duplicated(mon2h), ]
monkey2 <- merge(monkey, mon2h, by = 'Monkey_gene')
monkey3 <- monkey2[, -1]
monkey3 <- monkey3 %>%
  group_by(Human_gene) %>%
  summarise(across(everything(), sum)) %>%
  data.frame()
rownames(monkey3) <- monkey3[, 1]
monkey3 <- monkey3[, -1]
write10xCounts("monkey2human_matrix/", as(as.matrix(monkey3), "dgCMatrix"), version = "3")

## All species data
meta_in <- read.csv('E:/R/DRG/Fig1/all_spe_meta.csv',row.names = 1)
mhmrp <- Read10X('mm_all_species/')
nms <- rownames(mhmrp)
pre_nr_term <- c("^ERCC-","^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBB","^MTATP","^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR","^MCM", '^ERCC')
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("MALAT1", "XIST", "XIST_intron")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
print("pre_bad_genes")
print(length(pre_bad_genes))

mhmrp2 <- mhmrp[setdiff(nms, pre_bad_genes),]

spe_inter <- CreateSeuratObject(counts = mhmrp[, intersect(meta_in$cellid, colnames(mhmrp))])
spe_inter
me <- spe_inter@meta.data
me$cellid <- rownames(me)
spe_inter@meta.data <- merge(me, meta_in, by='cellid')
rownames(spe_inter@meta.data) <- spe_inter@meta.data$cellid
head(spe_inter@meta.data)
ucolkey <-  unique(meta_in[, c('anno', 'cellid')])
col2grp <- ucolkey$meta_in
names(col2grp) <- ucolkey$cellid
col2grp <- factor(col2grp)
spe_inter@active.ident <- col2grp
levels(spe_inter@active.ident)
spe_inter <- NormalizeData(spe_inter)
spe_inter <- FindVariableFeatures(spe_inter, selection.method = "vst")
spe_inter <- RunFastMNN(object.list = SplitObject(spe_inter, split.by = "Datast"))
spe_inter <- RunUMAP(spe_inter, reduction = "mnn", dims = 1:10)
spe_inter <- FindNeighbors(spe_inter, reduction = "mnn", dims = 1:35)
spe_inter <- FindClusters(spe_inter)
DimPlot(spe_inter, group.by = c("Datast", "Tissue", 'anno'), ncol = 3)

spe_inter.list <- SplitObject(spe_inter, split.by = "Datast")

# normalize and identify variable features for each dataset independently
spe_inter.list <- lapply(X = spe_inter.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
 hvgs_per_dataset <- lapply(spe_inter.list, function(x) {
   x@assays$RNA@var.features
 })
features <- SelectIntegrationFeatures(object.list = spe_inter.list)
spe_inter.anchors <- FindIntegrationAnchors(object.list = spe_inter.list)

spe_inter.integrated <- IntegrateData(anchorset = spe_inter.anchors, dims = 1:25)

DefaultAssay(spe_inter.integrated) <- "integrated"

features = row.names(spe_inter.integrated)

spe_inter.integrated <- ScaleData(spe_inter.integrated, verbose = FALSE)

spe_inter.integrated <- RunPCA(spe_inter.integrated, npcs = 20, verbose = FALSE)
ElbowPlot(spe_inter.integrated)
spe_inter.integrated <- RunUMAP(spe_inter.integrated, reduction = "pca", dims = 1:5)

spe_inter.integrated@active.ident = col2grp
levels(spe_inter.integrated@active.ident)
DimPlot(spe_inter.integrated, reduction = "umap", group.by = c("Datast", "Tissue", 'anno'), ncol = 3,
        label = TRUE, repel = TRUE) + NoLegend()
umap <- as.data.frame(spe_inter.integrated@reductions[["umap"]]@cell.embeddings)
colnames(umap) <- c('umap1', 'umap2')
umap$cellid <- rownames(umap)
meta <- merge(umap, meta_in, by = 'cellid')
rownames(meta) <- meta$cellid
write.csv(meta, 'E:/R/DRG/Fig1/all_spe_intergrated.csv')



## Fig1f ------------------------------------------2D plot-------------------------------------------------
meta <- read.csv('E:/R/DRG/Fig1/all_spe_intergrated.csv',row.names = 1)

p <- ggplot(meta, aes(x = mnn1, y = mnn2, color = anno)) + 
    geom_point(shape = 19, size = 1.5,  stroke = 0.1) +
    scale_color_manual(values=c('royalblue', "hotpink1")) +
    drg_theme()+
    theme(legend.position ='none')
p
ggsave(p,filename = './all.pdf',width = 4.5, height = 4, units = "in", device = 'pdf',dpi = 300) 


sp <- 'hg38'
p <-  ggplot(meta)+
     geom_point(aes(x = mnn1, y = mnn2),shape = 19,size = 1.5,color='grey88')+
     geom_point(data = meta[meta$Datast == sp & meta$Tissue == 'DRG', ],
                aes(x = mnn1, y = mnn2,color = anno),
                shape = 19,size = 1.5, stroke = 0.1)+
      scale_color_manual(values = c('royalblue1','hotpink1'))+
      labs(title = sp)+
      theme(legend.position ='none')+
      drg_theme() 
p
ggsave(paste(sp, 'drg.pdf'),p,width = 4,height = 4)


## Fig1g ------------------------------------------cor plot-------------------------------------------------
corr <- read.csv('E:/R/DRG/Fig1/cor_plot.csv', row.names = 1)
col <- rev(COL2('RdBu', 10))
pdf('./cor_plot.pdf',width = 12, height = 10)
corrplot(corr, 
         order = 'hclust', 
         addrect = 2, 
         method = "square",
         col = rev(COL2('RdBu', 10)))
dev.off()


## Fig1i -------------------------------------------DEGs--------------------------------------------------------
#绘图
allmarkers <- read.csv('E:/R/DRG/Fig1/all_spe_dge.csv')
allup <- filter(allmarkers, allmarkers$Significant == 'Up')
alldown <-filter(allmarkers, allmarkers$Significant == 'Down')
allud <- read.csv('E:/R/DRG/Fig1/all_label.csv')
p <- ggplot(allmarkers,
           aes(x = pct.2 - pct.1,y = avg_log2FC)) +
  geom_point(color = 'grey80', size = 0.3) +
  geom_hline(yintercept = c(-0.25,0.25), lty = 'dashed', size = 1, color = 'grey50') +

  geom_text_repel(data = allup,
                  aes(x = pct.2 - pct.1, y = avg_log2FC,
                      label = gene, color = Significant),
                  show.legend = F, direction = 'y',
                  hjust = 0, 
                  nudge_y = 0.25, 
                  force = 5, 
                  nudge_x =1- (allup$pct.1 - allup$pct.2)) +

  geom_text_repel(data = alldown,
                  aes(x = pct.2 - pct.1,y = avg_log2FC,
                      label = gene,color = Significant),
                  show.legend = F,direction = 'y',
                  hjust = 1, 
                  force = 2.5, 
                  nudge_x = -1- (alldown$pct.1 - alldown$pct.2)) +
 
  geom_point(data = allud, show.legend = F,
             aes(x = pct.2 - pct.1, y = avg_log2FC, color = Significant),
             size = 0.3) +
  scale_color_manual(values = c('royalblue', 'hotpink1')) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_rect(color = NA,fill = 'grey90')) +
  xlab(expression(Delta~'Percentage Diffrence')) +
  ylab('Log2-Fold Change') +
  facet_wrap(~Group, nrow = 1,scales = 'fixed')
p
ggsave(p, filename = './deg_all_species.pdf',width = 18, height = 15, units = "in", device='pdf',dpi=300)




