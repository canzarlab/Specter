library(splatter)
library(scater)
library(Rtsne)
library(mclust)
library(Seurat)

######################
n_cells <- 10000
n_genes <- 1000
params.groups <- newSplatParams(batchCells = n_cells, nGenes = n_genes)
n_groups <- 2
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.01, n_groups)

group.prob <- c(90, 10)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.65,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

n_reps <- 20

for (rare_pct in c(1, 2.5, 5))
{
  for (n_run in 1:10)
  {
    set.seed(100*rare_pct + n_run) # avoid repetitive subsample
    n_rare_pop <- round(n_cells * rare_pct / (2*(100-rare_pct)))
    group1 <- which(sim1$Group %in% "Group1") #sim1$Group=="Group1"
    group2 <- which(sim1$Group %in% "Group2") # sim1$Group=="Group2"
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
    
    dat <- counts(subsim1)
    pbmc <- CreateSeuratObject(counts = dat, project = "exp2", min.cells = 0, min.features = 0) #pbmc.data
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    # HVG
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    # Scaling the data
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    Data1 <- pbmc@reductions[["pca"]]@cell.embeddings
    x <- as.factor(subsim1$Group)
    levels(x) <- 1:length(levels(x))
    labels <- as.numeric(x)
    
    ## save subsample data
    simdatName <- paste0("Sim_rarecellExp2_pct_", toString(rare_pct), "_run_", toString(n_run))
    write.table(Data1,file = paste0("rare_cell_data/", simdatName, "_pca.csv"),
                sep = ",", row.names = F, col.names = F)
    write.table(labels,file = paste0("rare_cell_data/", simdatName, "_labels.csv"),
                sep = ",", row.names = F, col.names = F)

  }
}

