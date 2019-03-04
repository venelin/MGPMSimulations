# Calculates the mean vector mu and the variance covariance matrix Sigma
# for the (N x k) - variate Gaussian distribution at the tips of each tree, expected
# under the true (simulated) model of evolution for each one of the 512 parameter
# sets in testData_t5. Stores the resulting mean vectors and covariance matrices
# in a data.table called trueModelDistributions_t5 and stores this data.table
# in a file trueModelDistributions_t5.RData
library(PCMBase)
library(MGPMSimulations)
library(data.table)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

trueModelDistributions_t5 <-
  testData_t5[, {
    m <- model[[1]]
    tr <- treeWithRegimes[[1]]
    metaI <- PCMInfo(X = NULL, tree = tr, model = m)
    list(
      mu = list({
        as.vector(PCMMean(tree = tr, model = m, metaI = metaI))
      }),
      Sigma = list({
        PCMVar(tree = tr, model = m, metaI = metaI)
      })
    )}, keyby = list(IdTree, IdClusteringForTree, IdMappingForClustering, IdParamForMapping)]

save(trueModelDistributions_t5, file = "trueModelDistributions_t5.RData")
