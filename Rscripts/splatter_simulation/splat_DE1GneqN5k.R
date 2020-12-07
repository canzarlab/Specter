rm(list = ls())
library(splatter)
library(scater)
library(Rtsne)
library(mclust)

######################
n_genes <- 1000
params.groups <- newSplatParams(batchCells = 5000, nGenes = n_genes)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.01, n_groups)

group.prob <- c(1, 5, 14, 30, 50)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.4,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

## save data
simdatName <- "DE1GneqN5k"
if (n_genes==10000){
  simdatName <- "DE1GneqN5kD10k" 
}
colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/", simdatName, ".rds"))
