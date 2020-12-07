rm(list = ls())
library(splatter)
library(scater)

######################
n_genes <- 1000
n_cells <- 100000 ## input number of cells here
params.groups <- newSplatParams(batchCells = n_cells, nGenes = n_genes)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.01, n_groups)

group.prob <- c(1, 5, 14, 30, 50)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.7,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)

colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/big_data_100k.rds"))

