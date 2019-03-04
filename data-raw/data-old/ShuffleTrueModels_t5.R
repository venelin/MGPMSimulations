library(PCMFit)
library(PCMBase)
library(PCMBaseCpp)
library(MGPMSimulations)
library(data.table)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

fits_ShuffledTrueModels_t5 <-
  testData_t5[
    ,
    list(
      IdParamForMapping, IdGlob,
      clusterNodes, mapping,
      model = {
        print(.N)
        samp <- c(.N/2,  1:(.N/2 - 1), .N, (.N/2 + 1):(.N-1))
        #samp <- seq_len(.N)
        print(samp)
        model[samp]
      },
      treeWithRegimes,
      X
    ),
    keyby = list(
      IdTree, IdClusteringForTree,
      IdMappingForClustering, IdSimulationForParam)]

setkey(fits_ShuffledTrueModels_t5,
       IdTree,
       IdClusteringForTree,
       IdMappingForClustering,
       IdParamForMapping,
       IdSimulationForParam,
       IdGlob)

fits_ShuffledTrueModels_t5[
  ,
  logLik:=sapply(.I, function(i) {
    PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]],
           metaI = PCMInfoCpp)
  })]

fits_ShuffledTrueModels_t5[
  ,
  AIC:= sapply(.I, function(i) {
    ll <- logLik[i]
    df <- PCMParamCount(model[[i]],
                        countRegimeChanges = TRUE, countModelTypes = TRUE)
    -2*ll + 2*df
  })]

fits_ShuffledTrueModels_t5[, c("X", "treeWithRegimes"):=NULL]

usethis::use_data(fits_ShuffledTrueModels_t5, overwrite = TRUE)
