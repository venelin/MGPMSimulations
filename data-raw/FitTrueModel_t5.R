# This script can be run locally or on a cluster using a command like:
#
# bsub -n 1 -W 2:00 sh R --vanilla --slave -f ../../FitTrueModel_t5.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package MGPMSimulations The blank XXX stays for the type of model to be inferred, e.g. ScalarOU or MixedGaussian
#
library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(MGPMSimulations)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
} else {
  id <- 1548
}

prefixFiles = paste0("FitTrueModel_t5_2_id_", id)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-7)

print(PCMOptions())

tree <- testData_t5$treeWithRegimes[[id]]
values <- testData_t5$X[[id]][, 1:PCMTreeNumTips(tree)]
model <- testData_t5$model[[id]]
vecParams <- PCMParamGetShortVector(model, k = PCMNumTraits(model), R = PCMNumRegimes(model))

matParInit1 <- matrix(vecParams, nrow = 1L)
matParInitGuess <- GuessInitVecParams(
  model, n = 10000, X = values, tree = tree)
matParInitGuessVaryTheta <- GuessInitVecParams(
  model, n = 10000, X = values, tree = tree, varyTheta = TRUE)
matParInitJitter <- jitterModelParams(
  model,
  numJitterRootRegimeFit = 10000,  sdJitterAllRegimeFits = 0.05,
  numJitterAllRegimeFits = 10000, sdJitterRootRegimeFit = 0.05)


fit <- PCMFit(

  X = values, tree = tree, model = model, positiveValueGuard = 10000,

  metaI = PCMInfoCpp,

  matParInit = rbind(matParInit1,
                     matParInitGuess,
                     matParInitGuessVaryTheta,
                     matParInitJitter),

  numCallsOptim = 60L,
  numRunifInitVecParams = 0L,
  numGuessInitVecParams = 0L,

  verbose = FALSE)

save(fit, file = paste0("FinalResult_", prefixFiles, ".RData"))
