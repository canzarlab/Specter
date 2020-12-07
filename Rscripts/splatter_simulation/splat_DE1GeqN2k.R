rm(list = ls())
library(splatter)
library(scater)

######################
params.groups <- newSplatParams(batchCells = 2000, nGenes = 1000)
n_groups <- 5
# de.prob - DE probability This parameter controls the probability that a gene will be selected to be differentially expressed.
# The higher the more separation between groups
deprob <- rep(0.01, n_groups)
group.prob <- c(20, 20, 20, 20, 20)/100

sim1 <- splatSimulateGroups(params.groups, group.prob = group.prob, de.prob = deprob, de.facLoc = 1.1,
                            verbose = FALSE)
sim1 <- logNormCounts(sim1)
sim1 <- runPCA(sim1)

## save data
simdatName <- "DE1GeqN2k"
colData(sim1)$phenoid <- sim1$Group
saveRDS(sim1, paste0("splat_data/", simdatName, ".rds"))


