rm(list = ls())

setwd("/data/hoan/re_analyze/multimodal")
library(Seurat)
library(ggplot2)
library(mclust)
library(plyr)
library(dplyr) 

# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative controls
# for the protein measurements. For this reason, the gene expression matrix has HUMAN_ or MOUSE_
# appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "../seurat/data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
# export to file 
# write.table(t(as.matrix(cbmc.rna)), file = paste0("data/cbmc_raw_rna_count.csv"), row.names = F, col.names = F, sep = ',' )

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "../seurat/data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))

# When adding multimodal data to Seurat, it's okay to have duplicate feature names. Each set of
# modal data (eg. RNA, ADT, etc.) is stored in its own Assay object.  One of these Assay objects
# is called the 'default assay', meaning it's used for all analyses and visualization.  To pull
# data from an assay that isn't the default, you can specify a key that's linked to an assay for
# feature pulling.  To see all keys for all objects, use the Key function.  Lastly, we observed
# poor enrichments for CCR5, CCR7, and CD10 - and therefore remove them from the matrix
# (optional)
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]

#-------------------------------------------------------------------------------------------------
#################################### Analyse RNA-seq ######################################
#-------------------------------------------------------------------------------------------------

cbmc <- CreateSeuratObject(counts = cbmc.rna)
# standard log-normalization
cbmc <- NormalizeData(cbmc)

cbmc <- FindVariableFeatures(cbmc)
# standard scaling (no regression)
cbmc <- ScaleData(cbmc)
# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
cbmc <- RunPCA(cbmc, verbose = FALSE)
ElbowPlot(cbmc, ndims = 50)
cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")
# Find the markers that define each cluster, and use these to annotate the clusters, we use
# max.cells.per.ident to speed up the process
cbmc.rna.markers <- FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)

new.cluster.ids <- c("CD4 T", "CD14+ Mono", "CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs")

names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
cbmc.rna.markers <- FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)

# plot
DimPlot(cbmc, label = TRUE) + NoLegend()

#-------------------------------------------------------------------------------------------------
################ Add the protein expression levels to the Seurat object ####################
#-------------------------------------------------------------------------------------------------
# We will define an ADT assay, and store raw counts for it

# If you are interested in how these data are internally stored, you can check out the Assay
# class, which is defined in objects.R; note that all single-cell expression data, including RNA
# data, are still stored in Assay objects, and can also be accessed using GetAssayData
cbmc[["ADT"]] <- CreateAssayObject(counts = cbmc.adt)

# Now we can repeat the preprocessing (normalization and scaling) steps that we typically run
# with RNA, but modifying the 'assay' argument.  For CITE-seq data, we do not recommend typical
# LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# independently for each feature.  This is a slightly improved procedure from the original
# publication, and we will release more advanced versions of CITE-seq normalizations soon.
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR")
cbmc <- ScaleData(cbmc, assay = "ADT")

# Visualize protein levels on RNA clusters
# in this plot, protein (ADT) levels are on top, and RNA levels are on the bottom
FeaturePlot(cbmc, features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16", "CD3E", "ITGAX", "CD8A", 
                               "FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

# Let's plot CD4 vs CD8 levels in T cells
tcells <- subset(cbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")

cbmc <- subset(cbmc, idents = c("Multiplets", "Mouse"), invert = TRUE)

## cluster on protein
# Because we're going to be working with the ADT data extensively, we're going to switch the
# default assay to the 'CITE' assay.  This will cause all functions to use ADT data by default,
# rather than requiring us to specify it each time
DefaultAssay(cbmc) <- "ADT"
cbmc <- RunPCA(cbmc, features = rownames(cbmc), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
DimPlot(cbmc, reduction = "pca_adt")

# Since we only have 10 markers, instead of doing PCA, we'll just use a standard euclidean
# distance matrix here.  Also, this provides a good opportunity to demonstrate how to do
# visualization and clustering using a custom distance matrix in Seurat
adt.data <- GetAssayData(cbmc, slot = "data")
adt.dist <- dist(t(adt.data))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
cbmc[["rnaClusterID"]] <- Idents(cbmc)

# Now, we rerun tSNE using our distance matrix defined only on ADT (protein) levels.
cbmc[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
cbmc[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
cbmc <- FindClusters(cbmc, resolution = 0.2, graph.name = "adt_snn")

# We can compare the RNA and protein clustering, and use this to annotate the protein clustering
# (we could also of course use FindMarkers)
clustering.table <- table(Idents(cbmc), cbmc$rnaClusterID)
clustering.table

new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "CD34+", "T/Mono doublets", 
                     "CD16+ Mono", "pDCs", "B")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)
cbmc[["adtClusterID"]] <- Idents(cbmc)

tsne_rnaClusters <- DimPlot(cbmc, reduction = "tsne_adt", group.by = "rnaClusterID") + NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on scRNA-seq") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 4)

tsne_adtClusters <- DimPlot(cbmc, reduction = "tsne_adt", pt.size = 0.5) + NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "ident", size = 4)
# Note: for this comparison, both the RNA and protein clustering are visualized on a tSNE
# generated using the ADT distance matrix.
CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

#-------------------------------------------------------------------------------------------------
################ Analyse Specter clustering ####################
# Cell types: 
#-------------------------------------------------------------------------------------------------

sce_labels <- read.csv("/data/hoan/Specter/output/cbmc_specter_adtK_15_rnaK_16_labels_v4.csv", header = F) #the best
new.cluster.ids <- c("NK","NK","pDCs","Eryth","CD8 T","CD8 T","CD14+ Mono","NK","NK","NK","CD14+ Mono","Eryth","CD4 T","DC","CD34+","CD16+ Mono", "MK", "B")

sce_labels <- as.factor(sce_labels$V1)
names(sce_labels) <- names(cbmc$rnaClusterID)
cbmc_sce <- cbmc 
sce_labels <- mapvalues(sce_labels, from = levels(sce_labels), to = new.cluster.ids)

cbmc_sce[["sceClusterID"]] <- sce_labels
DefaultAssay(cbmc_sce) <- "RNA"
Idents(cbmc_sce) <- sce_labels

specter.rna.markers <- FindAllMarkers(cbmc_sce, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
dc_markers <- specter.rna.markers[specter.rna.markers['cluster'] == 'DC', ]
mk_markers <- specter.rna.markers[specter.rna.markers['cluster'] == 'MK', ]


## identify NK subpopulations: CD56 dim vs CD56 bright
# RidgePlot(cbmc_sce, assay = "ADT", features = c("adt_CD56", "adt_CD16"), ncol = 1)
### Assign the labels
#Plot # tsne_RNA
tsne_rnaClusters <- DimPlot(cbmc_sce, reduction = "tsne", group.by = "sceClusterID") + NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on sceClusterID") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "sceClusterID", size = 4)
tsne_adtClusters <- DimPlot(cbmc_sce, reduction = "tsne", group.by = "rnaClusterID", pt.size = 0.5) + NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on RNA signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "rnaClusterID", size = 4)
CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

# tsne_adt
tsne_rnaClusters <- DimPlot(cbmc_sce, reduction = "tsne_adt", group.by = "sceClusterID") + NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering based on sceClusterID") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "sceClusterID", size = 4)
tsne_adtClusters <- DimPlot(cbmc_sce, reduction = "tsne_adt", group.by = "adtClusterID", pt.size = 0.5) + NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering based on ADT signal") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "adtClusterID", size = 4)
CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

########### Summary for paper #########################
#Plot clusters on RNA and ADT t-SNE. 
tsne_rnaClusters <- DimPlot(cbmc_sce, reduction = "tsne", group.by = "sceClusterID") + NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Citefuse's clusters on RNA") + theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "sceClusterID", size = 4)
tsne_adtClusters <- DimPlot(cbmc_sce, reduction = "tsne_adt", group.by = "sceClusterID", pt.size = 0.5) + NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Citefuse's clusters on ADT") + theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "sceClusterID", size = 4)
CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

############################################################################################################
##-------------------------------------- tSNE PLot for paper ---------------------------------------------

#---------------------------------- Multimodal Specter on ADT-based tSNE ----------------------------------------------
sce_labels <- read.csv("/data/hoan/Specter/output/cbmc_specter_adtK_15_rnaK_16_labels_v4.csv", header = F) #the best
new.cluster.ids <- c("NK","NK","pDCs","Eryth","CD8 T","CD8 T","CD14+ Mono","NK","NK","NK","CD14+ Mono","Eryth","CD4 T","DC","CD34+","CD16+ Mono", "MK", "B")
sce_labels <- as.factor(sce_labels$V1)
names(sce_labels) <- names(cbmc$rnaClusterID)
sce_labels <- mapvalues(sce_labels, from = levels(sce_labels), to = new.cluster.ids)
cbmc_sce <- cbmc
tsne_data <- cbmc_sce@reductions[["tsne_adt"]]@cell.embeddings
data <- data.frame(adtTSNE_1 = tsne_data[,1], adtTSNE_2 = tsne_data[,2], group = sce_labels)
data_nk <- subset(data, group=='MK')
data_dc <- subset(data, group=='DC')

p <- ggplot(data = data) +
  # blue plot
  geom_point(data=data, aes(x=adtTSNE_1, y=adtTSNE_2, col=group), size = 0.2) +
  geom_point(data=data_nk, aes(x=adtTSNE_1, y=adtTSNE_2), size = 2.0, color = "red") +
  geom_point(data=data_dc, aes(x=adtTSNE_1, y=adtTSNE_2), size = 2.0, color = "gold")

p <- p+  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p <- p + theme(legend.position="none")  #remove legend
p <- p + ggtitle("Joint clusters (CBMC)") + theme(plot.title = element_text(size = 14, hjust = 0.5)) #, face = 'bold'
#p <- LabelClusters(plot = p, id = "group", size = 5, color = "black", fontface = 'bold')
p
# ggsave("results/cbmc_specter_on_adt.pdf")
# ggsave(filename = "results/cbmc_specter_on_adt_K18.pdf", width = 5, height = 4, dpi = 300, units = "in", device='pdf')
ggsave(filename = "results/cbmc_specter_on_adt_K18_nocap.pdf", width = 5, height = 4, dpi = 300, units = "in", device='pdf')


#---------------------------------- Multimodal Specter on RNA-based tSNE ----------------------------------------------
sce_labels <- read.csv("/data/hoan/Specter/output/cbmc_specter_adtK_15_rnaK_16_labels_v4.csv", header = F) #the best
new.cluster.ids <- c("NK","NK","pDCs","Eryth","CD8 T","CD8 T","CD14+ Mono","NK","NK","NK","CD14+ Mono","Eryth","CD4 T","DC","CD34+","CD16+ Mono", "MK", "B")
sce_labels <- as.factor(sce_labels$V1)
names(sce_labels) <- names(cbmc$rnaClusterID)
sce_labels <- mapvalues(sce_labels, from = levels(sce_labels), to = new.cluster.ids)
tsne_data <- cbmc_sce@reductions[["tsne"]]@cell.embeddings
data <- data.frame(TSNE_1 = tsne_data[,1], TSNE_2 = tsne_data[,2], group = sce_labels)
data_nk <- subset(data, group=='MK')
data_dc <- subset(data, group=='DC')

p <- ggplot(data = data) +
  # blue plot
  geom_point(data=data, aes(x=TSNE_1, y=TSNE_2, col=group), size = 0.4) +
  geom_point(data=data_nk, aes(x=TSNE_1, y=TSNE_2), size = 0.4, color = "red") +
  geom_point(data=data_dc, aes(x=TSNE_1, y=TSNE_2), size = 0.4, color = "gold")

p <- p+  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p <- p + theme(legend.position="none")  #remove legend
p <- p + ggtitle("Specter") + theme(plot.title = element_text(size = 16, hjust = 0.5))
#p <- LabelClusters(plot = p, id = "group", size = 4, color = "black", fontface = 'bold')
p
#ggsave("results/cbmc_specter_on_rna.pdf")
ggsave(filename = "results/cbmc_specter_on_rna_K18_nocaption.pdf", width = 6, height = 5, dpi = 300, units = "in", device='pdf')


#---------------------------------- Multimodal Citefuse on RNA-based tSNE ----------------------------------------------
# sce_labels <- read.csv("/data/hoan/re_analyze/multimodal/input/cbmc_citefuse_15labels.csv", header = F)
# new.cluster.ids <- c("CD4 T","CD14+ Mono","B","NK","CD34+","CD16+ Mono","Eryth","pDCs",
#                      "CD8 T","NK","T/Mono doublets","NK","DC","B","CD4 T")
sce_labels <- read.csv("/data/hoan/re_analyze/multimodal/input/cbmc_citefuse_18labels.csv", header = F)
# new.cluster.ids <- c("CD4 T","CD14+ Mono","B","NK","NK","CD34+","B","pDCs","NK","T/Mono doublets","NK","CD8 T","CD16+ Mono","NK","Eryth","CD4 T", "DC", "CD4 T")
new.cluster.ids <- c("CD4 T","CD14+ Mono","B","NK","NK","CD34+","B","pDCs","NK","CD14+ Mono","NK","CD8 T","CD16+ Mono","NK","Eryth","CD4 T", "DC", "CD4 T")
sce_labels <- as.factor(sce_labels$V1)
names(sce_labels) <- names(cbmc$rnaClusterID)
sce_labels <- mapvalues(sce_labels, from = levels(sce_labels), to = new.cluster.ids)
tsne_data <- cbmc_sce@reductions[["tsne"]]@cell.embeddings
data <- data.frame(TSNE_1 = tsne_data[,1], TSNE_2 = tsne_data[,2], group = sce_labels)
# data_nk <- subset(data, group=='MK')
data_dc <- subset(data, group=='DC')

p <- ggplot(data = data) +
  # blue plot
  geom_point(data=data, aes(x=TSNE_1, y=TSNE_2, col=group), size = 0.4) +
  # + geom_point(data=data_nk, aes(x=TSNE_1, y=TSNE_2), size = 1, color = "red") +
  geom_point(data=data_dc, aes(x=TSNE_1, y=TSNE_2), size = 0.4, color = "gold")

p <- p+  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p <- p + theme(legend.position="none")  #remove legend
p <- p + ggtitle("CiteFuse") + theme(plot.title = element_text(size = 16, hjust = 0.5))
#p <- LabelClusters(plot = p, id = "group", size = 4, color = "black", fontface = 'bold')
p
# ggsave("results/cbmc_citefuse_on_rna.pdf")
ggsave(filename = "results/cbmc_citefuse_on_rna_K18_nocaption.pdf", width = 6, height = 5, dpi = 300, units = "in", device='pdf')


