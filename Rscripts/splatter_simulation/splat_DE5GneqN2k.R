rm(list = ls())
library(splatter)
library(scater)
library(Rtsne)

######################
params.groups <- newSplatParams(batchCells = 2000, nGenes = 1000)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.05, n_groups)

group.prob <- c(1, 5, 14, 30, 50)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.2,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

## save data
simdatName <- "DE5GneqN2k"
colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/", simdatName, ".rds"))
