# This script can be run locally or on a cluster using a command like:
#
# bsub -n 1 -W 2:00 sh R --vanilla --slave -f ../../FitTrueModel_NoHint_t5.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in
# testData_t5 data.table from the package MGPMSimulations
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

prefixFiles = paste0("FitTrueModel_NoHint_t5_id_", id)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-7)

print(PCMOptions())

tree <- testData_t5$treeWithRegimes[[id]]
values <- testData_t5$X[[id]][, seq_len(PCMTreeNumTips(tree))]

model <- MixedGaussian(
  k = 2,
  modelTypes = MGPMDefaultModelTypes(),
  mapping = testData_t5$mapping[[id]],
  Sigmae_x = Args_MixedGaussian_MGPMDefaultModelTypes()$Sigmae_x)

fit <- PCMFit(

  X = values, tree = tree, model = model, positiveValueGuard = 10000,

  metaI = PCMInfoCpp,

  numCallsOptim = 60L,
  numRunifInitVecParams = 20000L,
  numGuessInitVecParams = 20000L,

  verbose = FALSE)

save(fit, file = paste0("FinalResult_", prefixFiles, ".RData"))
