library(mclust)
library(SingleCellExperiment)
library(scater)
library(DuoClustering2018)
library(Seurat)

## download dataset from DuoClustering2018
sim1 <- readRDS("data/sce_full/sce_full_Zhengmix4eq.rds")

cell_type1 <- "regulatory.t" 
cell_type2 <- "naive.cytotoxic" 

sim1 <- sim1[,(sim1$phenoid==cell_type1) | (sim1$phenoid==cell_type2)]
n_cells <- dim(sim1)[2]
n_reps <- 10
n_groups <- 2

for (rare_pct in c(1, 2.5, 5, 10))
{
  for (n_run in 1:10)
  {
    set.seed(100*rare_pct + n_run) # avoid repetitive subsample
    n_rare_pop <- round(n_cells * rare_pct / (2*(100-rare_pct)))
    group1 <- which(sim1$phenoid %in% cell_type1) 
    group2 <- which(sim1$phenoid %in% cell_type2) 
    subgroup1 <- sample(group1, n_rare_pop)
    subgroup2 <- sample(group2, n_rare_pop)
    print(subgroup1)
    print(subgroup2)
    
    if (n_run <= n_reps/2)
    {
      subsim1 <- sim1[, c(group1, subgroup2)]
    }  else {
      subsim1 <- sim1[, c(subgroup1, group2)]
    }
    
    ## save subsample data
    simdatName <- paste0("Sim_rarecellExp3_pct_", toString(rare_pct), "_run_", toString(n_run))
        
    sce <- subsim1
    x <- as.factor(sce$phenoid)
    levels(x) <- 1:length(levels(x))
    labels <- as.numeric(x)
    dat <- logcounts(sce)
    n_HV <- 2000
    if (dim(dat)[1] > n_HV)
    {
        sel = order(apply(dat, 1, var), decreasing=TRUE)[1: n_HV]
        dat <- dat[sel, ]
    }
    
    # PCA data
    pca <- prcomp(t(dat), center = TRUE, scale. = FALSE)
    # print(dim(dat))
    N <- min(100, dim(dat)[2]-1) #handle data with smaller than 100 cells
    Data1 <- pca$x[, 1:N, drop = FALSE]
    write.table(Data1,file = paste0("rare_cell_data/", dataname, "pca.csv"),
                sep = ",", row.names = F, col.names = F)
    write.table(labels,file = paste0("rare_cell_data/", dataname, "_labels.csv"),
                sep = ",", row.names = F, col.names = F)

  }
}

