library(PCMFit)
library(PCMBase)
library(MGPMSimulations)
library(data.table)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

fits_ModelFWithTrueShiftsOnPRCData_t5 <- rbindlist(
  lapply(
    seq_len(nrow(testData_t5)),
    function(id) {
      if(id %% 10L == 0) {
        cat("id:", id, "\n")
      }

      dfModelF <- NA_integer_

      prefixFiles = paste0("FitModelFWithTrueShiftsOnPRCData_t5_id_", id)
      file <- paste0("Results_FitModelFWithTrueShiftsOnPRCData_t5/", prefixFiles,
                     "/FinalResult_", prefixFiles, ".RData")
      if(file.exists(file)) {
        load(file)
        fit_1 <- fit
        fitPRC_1 <- fitPRC
        dfModelF <- PCMParamCount(fit_1$modelOptim)
      } else {
        fit_1 <- structure(list(0.0), class = "NAFIT")
        fitPRC_1 <- structure(list(0.0), class = "NAFIT")
      }

      file <- paste0("Results_FitModelFWithTrueShiftsOnPRCData_t5_2/", prefixFiles,
                     "/FinalResult_", prefixFiles, ".RData")
      if(file.exists(file)) {
        load(file)
        fit_2 <- fit
        fitPRC_2 <- fitPRC
        dfModelF <- PCMParamCount(fit_2$modelOptim)
      } else {
        fit_2 <- structure(list(0.0), class = "NAFIT")
        fitPRC_2 <- structure(list(0.0), class = "NAFIT")
      }


      file <- paste0("Results_FitModelFWithTrueShiftsOnPRCData_t5_3/", prefixFiles,
                     "/FinalResult_", prefixFiles, ".RData")
      if(file.exists(file)) {
        load(file)
        fit_3 <- fit
        fitPRC_3 <- fitPRC
        dfModelF <- PCMParamCount(fit_3$modelOptim)
      } else {
        fit_3 <- structure(list(0.0), class = "NAFIT")
        fitPRC_3 <- structure(list(0.0), class = "NAFIT")
      }

      logLik.NAFIT <- AIC.NAFIT <- function(object, ...) as.double(NA)

      testData_t5[id,
                  list(
                    IdGlob, IdTree, IdClusteringForTree,
                    IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                    clusterNodes, mapping,
                    #model,
                    logLik = unlist(logLik),
                    nobs,
                    dfModel = sapply(model, PCMParamCount),
                    dfModelF = dfModelF,
                    #model_1 = list(fit_1$modelOptim),
                    #modelPRC_1 = list(fitPRC_1$modelOptim),
                    #model_2 = list(fit_2$modelOptim),
                    #modelPRC_2 = list(fitPRC_2$modelOptim),
                    logLik_fit_1 = logLik(fit_1), AIC_fit_1 = AIC(fit_1),
                    logLik_fitPRC_1 = logLik(fitPRC_1), AIC_fitPRC_1 = AIC(fitPRC_1),
                    logLik_fit_2 = logLik(fit_2), AIC_fit_2 = AIC(fit_2),
                    logLik_fitPRC_2 = logLik(fitPRC_2), AIC_fitPRC_2 = AIC(fitPRC_2),
                    logLik_fit_3 = logLik(fit_3), AIC_fit_3 = AIC(fit_3),
                    logLik_fitPRC_3 = logLik(fitPRC_3), AIC_fitPRC_3 = AIC(fitPRC_3))]
    }))


 usethis::use_data(fits_ModelFWithTrueShiftsOnPRCData_t5, overwrite = TRUE)
