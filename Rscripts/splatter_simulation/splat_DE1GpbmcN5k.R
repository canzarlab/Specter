rm(list = ls())
library(splatter)
library(scater)

######################
params.groups <- newSplatParams(batchCells = 5000, nGenes = 1000)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.01, n_groups)

# DC: 2, NK: 20, B: 10, Mono: 8, T: 60
group.prob <- c(2, 8, 10, 20, 60)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 0.3,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

## save data
simdatName <- "DE1GpbmcN5k"
colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/", simdatName, ".rds"))
