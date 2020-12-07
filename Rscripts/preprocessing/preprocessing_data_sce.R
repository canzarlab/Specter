library(SingleCellExperiment)
library(scater)
library(plyr)
library(stats)

process_data_HG <- function(data_dir, dataname){
  # Parameters:
  ## data_dir: input directory of data (SCE object)
  ## dataname: name of the dataset
  
  sce <- readRDS(data_dir)
  x <- as.factor(sce$cell_type1)
  levels(x) <- 1:length(levels(x))
  labels <- as.numeric(x)
  if(!("logcounts") %in% names(assays(sce)) && ("counts") %in% names(assays(sce)))
  {
    counts = assays(sce)$counts
    new_counts <- scale(counts, center = FALSE,
                        scale = colSums(counts))     # normalized by the feature, |Samples| = 1
    assays(sce)$logcounts <- log2(new_counts*(10^4)+1)
  }
  
  dat <- logcounts(sce)
  n_HV <- 2000
  if (dim(dat)[1] > n_HV)
  {
    sel = order(apply(dat, 1, var), decreasing=TRUE)[1: n_HV]
    dat <- dat[sel, ]
  }

  # PCA
  pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
  # print(dim(dat))
  N <- min(100, dim(dat)[2]-1) #handle data with smaller than 100 cells
  Data1 <- pca$x[, 1:N, drop = FALSE]
  write.table(Data1,file = paste0(dataname, "_HG2k_pca100.csv"),
              sep = ",", row.names = F, col.names = F)
  
  ## save labels for simulated data (uncomment below)
  # x <- as.factor(sim1$Group)
  # levels(x) <- 1:length(levels(x))
  # labels <- as.numeric(x)
  # write.table(labels,file = paste0(dataname, "_labels.csv"),
  #             sep = ",", row.names = F, col.names = F)
}

# An example
## Input data dir object
data_dir <- "../splatter_simulation/splat_data/DE1GpbmcN1k.rds"
## input data name
dataname <- "DE1GpbmcN1k"
## Run preprocessing
process_data_HG(data_dir, dataname)
