# This script can be run locally or on a cluster using a command like:
#
# bsub -n 1 -W 2:00 sh R --vanilla --slave -f ../../FitModelFWithTrueShiftsOnPRCData_t5.R --args 96
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

prefixFiles = paste0("FitModelFWithTrueShiftsOnPRCData_t5_id_", id)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-7)

print(PCMOptions())

tree <- testData_t5$treeWithRegimes[[id]]
values <- testData_t5$X[[id]][, seq_len(PCMTreeNumTips(tree))]
rownames(values) <- paste0("V", seq_len(nrow(values)))
valuesPRC <- t(prcomp(
  formula = ~V1+V2,
  data = as.data.table(t(values)))$x)

modelF <- MixedGaussian(
  k = 2,
  modelTypes = MGPMDefaultModelTypes(),
  mapping = rep(6, length(testData_t5$mapping[[id]])),
  Sigmae_x = Args_MixedGaussian_MGPMDefaultModelTypes()$Sigmae_x)

modelTrue <- testData_t5$model[[id]]
PCMParamSetByName(modelF, modelTrue, inplace = TRUE, deepCopySubPCMs = TRUE)

cat("==============================\n\n")
cat("logLik(trueModel): ", toString(PCMLik(values, tree, modelTrue)), ";\n")
cat("logLik(modelF): ", toString(PCMLik(values, tree, modelF)), ";\n")

options(MGPMSimulations.LowerLimitTheta = apply(values, 1, min))
options(MGPMSimulations.UpperLimitTheta = apply(values, 1, max))
options(MGPMSimulations.LowerSigma_xUppTri = -1.0)

randVecParams <- PCMParamRandomVecParams(
  modelF, k = PCMNumTraits(modelF), R = PCMNumRegimes(modelF))

PCMParamLoadOrStore(
  modelF, randVecParams, 0, k = PCMNumTraits(modelF), R = PCMNumRegimes(modelF),
  load=TRUE)

fit <- PCMFit(

  X = values, tree = tree, model = modelF, positiveValueGuard = 10000,

  metaI = PCMInfoCpp,

  numCallsOptim = 200L,
  numRunifInitVecParams = 100L,
  numGuessInitVecParams = 100000L,

  verbose = FALSE)

cat("logLik(fit): ", toString(logLik(fit)), ";\n")


options(MGPMSimulations.LowerLimitTheta = apply(valuesPRC, 1, min))
options(MGPMSimulations.UpperLimitTheta = apply(valuesPRC, 1, max))
options(MGPMSimulations.LowerSigma_xUppTri = -1.0)

fitPRC <- PCMFit(

  X = valuesPRC, tree = tree, model = modelF, positiveValueGuard = 10000,

  metaI = PCMInfoCpp,

  numCallsOptim = 200L,
  numRunifInitVecParams = 100L,
  numGuessInitVecParams = 100000L,

  verbose = FALSE)

save(fit, fitPRC, file = paste0("FinalResult_", prefixFiles, ".RData"))

cat("logLik(fitPRC): ", toString(logLik(fitPRC)), "\n")
cat("\n==============================\n")
