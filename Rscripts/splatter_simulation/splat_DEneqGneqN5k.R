rm(list = ls())
library(splatter)
library(scater)

######################
n_genes <- 10000
params.groups <- newSplatParams(batchCells = 5000, nGenes = n_genes)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- c(0.01, 0.01, 0.02, 0.02, 0.05)

group.prob <- c(1, 5, 14, 30, 50)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.65,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

## save data
if (n_genes==10000){
  simdatName <- "DEneqGneqN5kD10k"
} else {
  simdatName <- "DEneqGneqN5k"
}

colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/", simdatName, ".rds"))




