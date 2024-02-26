# Fig3 PNS microglia differentiate in parallel with CNS microglia during human prenatal development.
library(SingleR)
library(celldex)
library(pheatmap)
library(ggplot2)
library(Seurat)
library(dplyr)
library(reshape2)
library(CytoTRACE)
library(Monocle3)
library(ggalluvial)
library(ImageGP)
library(ClusterGVis)
setwd('E:/R/DRG/Fig3')


## Fig3F 
# --------------------------------------------Monocle3 cell order----------------------------------------------------------
fdg <- read.csv('E:/R/DRG/Fig3/fdg_meta.csv')
counts <- Read10X(data.dir = "E:/R/DRG/Fig3/dsspv2/")
counts <- as.matrix(counts[, fdg$cellid])
gene_ann <- as.data.frame(rownames(counts), row.names = rownames(counts))
colnames(gene_ann) <- c('gene_short_name')
rownames(fdg) <- fdg$cellid

cds <- new_cell_data_set(as.matrix(counts),
                         cell_metadata = fdg,
                         gene_metadata = gene_ann)
str(cds)
cds <- preprocess_cds(cds, num_dim = 20)
cds <- align_cds(cds)
cds <- reduce_dimension(cds)

xy <- fdg[,c('cellid', 'fdg1', 'fdg2')]
rownames(xy) <- xy[,1]
xy <- xy[,-1]
colnames(xy) <- c('UMAP_1','UMAP_2')
xy <- as.matrix(xy)
    # drg_seu, seurat object
drg_seu@reductions[["umap"]]@cell.embeddings = xy 
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- as.matrix(xy)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)

embed <- data.frame(Embeddings(drg_seu, reduction = "umap"))
embed <- subset(embed,  UMAP_1 < -10000 &  UMAP_2 < -17000)
root.cell <- rownames(embed)  
cds <- order_cells(cds, root_cells = root.cell)

    # order: the single cell order
order <- as.data.frame(cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]])
colnames(order) <- 'timeorder'
order$cellid <- rownames(order)
    # orderm: the MC order value
order <- merge(order, fdg[, c('cellid', 'mcid')], by = 'cellid')
order <- order[, -1]
orderm <- order %>%
  group_by(mcid) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  data.frame()



# ---------------------------------------Actual time----------------------------------------------------------------
dp = fdg
dp$day <- NA
dp$day <- ifelse(dp$Time %in% c('CS10', 'CS11', 'CS12', 'CS13'), 28, dp$day)
dp$day <- ifelse(dp$Time %in% c('CS14', 'CS15'), 35, dp$day)
dp$day <- ifelse(dp$Time %in% c('CS16', 'CS17'), 42, dp$day)
dp$day <- ifelse(dp$Time %in% c('CS18', 'CS19'), 49, dp$day)
dp$day <- ifelse(dp$Time %in% c('CS21', 'CS22'), 56, dp$day)
dp$day <- ifelse(dp$Time %in% c('CS23'), 59.5, dp$day)
dp$day <- ifelse(dp$Time %in% c('9W'), 63, dp$day)
dp$day <- ifelse(dp$Time %in% c('10W'), 70, dp$day)
dp$day <- ifelse(dp$Time %in% c('12W'), 84, dp$day)
dp$day <- ifelse(dp$Time %in% c('15W'), 105, dp$day)
dp$day <- ifelse(dp$Time %in% c('16W'), 112, dp$day)
dp$day <- ifelse(dp$Time %in% c('17W'), 119, dp$day)
dp$day <- ifelse(dp$Time %in% c('20W'), 140, dp$day)
dp$day <- ifelse(dp$Time %in% c('21W'), 147, dp$day)
dp$day <- ifelse(dp$Time %in% c('24W'), 168, dp$day)

    # dpm: the actual time of each MC
dp <- dp[, c('mcid', 'day')]
dpm <- dp %>%
  group_by(mcid) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  data.frame()
dpm <- merge(dpm, unique(dp[, c('mcid', 'anno')]), by = 'mcid')


# ------------------------------------Cell rank--------------------------------------------------------------------------
phe <- fdg$anno
phe <- as.character(phe)
names(phe) <- rownames(fdg)
mat_3k <- as.matrix(pbmc@assays$RNA@counts)
results <- CytoTRACE(mat = mat_3k)
plotCytoTRACE(results, phenotype = phe)
    
    # cyto: Differentiation capacity of single cell
cyto <- as.data.frame(results$CytoTRACE)
colnames(cyto) <- 'Rank'
cyto$cellid <- rownames(cyto)

cyto <- merge(cyto, fdg[, c('cellid', 'mcid')], by = 'cellid')
cyto <- cyto[, -1]
    # cytom: Differentiation capacity of MC
cytom <- cyto %>%
  group_by(mcid) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  data.frame()
cytom <- merge(cytom, unique(dp[, c('mcid', 'anno')]), by = 'mcid')

# -----------------------------------Module score------------------------------------------------------------------------
data <- as.data.frame(as.matrix(pbmc@assays$RNA@data))
microgenes <- c("OLFML3",	"C3",	"GAL3ST4",	"PLXDC2",	"ENTPD1",	"HTRA1",	
            "P2RY12",	"TMEM119",	"SIGLEC8",	"LHFPL2",	"A2M",	"SORL1",	"APBB1IP",	"PDPN",
            "LILRB4",	"FCGR1A",	"NBL1",	"DNASE2",	"TREM2",	"FCGR1B")
micro <- as.data.frame(t(data[microgenes, ]))
micro$cellid <- rownames(micro)
micro <- merge(micro, fdg[, c('cellid', 'mcid')], by = 'cellid')
rownames(micro) <- micro[, 1]
micro <- micro[, -1]
microm <- micro %>%
  group_by(mcid) %>%
  summarise(across(everything(), mean)) %>%
  data.frame()
rownames(microm) <- microm[, 1]
microm <- microm[, -1]

macgenes <- c("LYVE1",	"DAB2",	"CD14",	"F13A1",	"MRC1",	"FCGRT",	"CD36",
          "NCF2",	"FOLR2",	"IGFBP4",	"RGL1",	"MS4A4A",	"ITSN1",	"CD163",
          "PMP22",	"STAB1",	"WWP1",	"SCN9A",	"NRP1")
mac <- as.data.frame(t(data[macgenes, ]))
mac$cellid <- rownames(mac)
mac <- merge(mac, fdg[, c('cellid', 'mcid')], by = 'cellid')
rownames(mac) <- mac[, 1]
mac <- mac[, -1]
macm <- mac %>%
  group_by(mcid) %>%
  dplyr::summarise(across(everything(), mean)) %>%
  data.frame()
rownames(macm) <- macm[, 1]
macm <- macm[, -1]

ms <- merge(microm, macm, by = 'mcid')

# -----------------------------------merge the above objects and plot figs--------------------------------------------------
cor_merge <- merge(dpm, orderm, by = 'mcid')
cor_merge <- merge(cor_merge, cytom[, c('mcid', 'Rank')], by = 'mcid')
cor_merge <- merge(cor_merge, ms, by = 'mcid')
order_ct <- c('Macrophage progenitor', 
             'Microglia precursor', 'PNS microglia precursor',
             'Microglia', 'PNS microglia')
cor_merge$anno <- factor(cor_merge$anno, levels = order_ct)

cor_merge$Group <- NA
cor_merge$Group <- ifelse(cor_merge$anno == 'Macrophage progenitor', 'a', too$Group)
cor_merge$Group <- ifelse(cor_merge$anno %in% c('Microglia precursor', 'Microglia'), 'b', too$Group)
cor_merge$Group <- ifelse(cor_merge$anno %in% c('PNS microglia precursor', 'PNS microglia'), 'c', too$Group)
write.csv(cor_merge, './cor_figs.csv')

    # plot
p = ggplot(cor_merge, aes(x = timeorder, y = day)) +
    geom_point(aes(x = timeorder, y = day , fill = anno), color = 'grey40', shape = 21, size = 8) +
    scale_fill_manual(values=rev(c('#921D39', '#9DB81F', '#D25690', "#E8B885", '#3AB69C'))) +
    geom_smooth(data = subset(cor_merge, Group %in% c("a", 'c')), method = lm, fill = "lightpink", color = "red") +
    geom_smooth(data = subset(cor_merge, Group %in% c("a", 'b')), method = lm, fill = "#B2DF8A", color = "#33A02C") +
    stat_cor(method = "pearson",
             data = subset(cor_merge, Group %in% c("a", 'c')), color = 'red') +
    stat_cor(method = "pearson", 
             data = subset(cor_merge, Group %in% c("a", 'b')), color = '#33A02C') +
    drg_theme()
p 
ggsave(p, filename='./monocle3_day.pdf',width = 6.5, height = 4.8, dpi=300, units = "in", device='pdf')





