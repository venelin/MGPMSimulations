library(PCMFit)
library(PCMBase)
library(MGPMSimulations)
library(data.table)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

fits_TrueModels_t5 <- rbindlist(
  lapply(
    seq_len(nrow(testData_t5)),
    function(id) {
      if(id %% 10L == 0) {
        cat("id:", id, "\n")
      }

      prefixFiles = paste0("FitTrueModel_t5_id_", id)
      file <- paste0("Results_t5_FitTrueModel_t1/", prefixFiles,
                     "/FinalResult_", prefixFiles, ".RData")
      if(file.exists(file)) {
        load(file)
        fit_1 <- fit
      } else {
        fit_1 <- structure(0.0, class = "NAFIT")
      }

      prefixFiles = paste0("FitTrueModel_t5_2_id_", id)
      file <- paste0("Results_t5_FitTrueModel_t2/", prefixFiles,
                     "/FinalResult_", prefixFiles, ".RData")
      if(file.exists(file)) {
        load(file)
        fit_2 <- fit
      } else {
        fit_2 <- structure(0.0, class = "NAFIT")
      }

      logLik.NAFIT <- AIC.NAFIT <- function(object, ...) as.double(NA)

      fit <- if(!is.na(logLik(fit_1)) && logLik(fit_1) > logLik(fit_2)) fit_1 else fit_2

      fit$modelOptim

      testData_t5[id,
                  list(
                    IdGlob, IdTree, IdClusteringForTree,
                    IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                    clusterNodes, mapping,
                    model = list(fit$modelOptim),
                    logLik = logLik(fit), AIC = AIC(fit),
                    logLik_fit1 = logLik(fit_1), AIC_fit1 = AIC(fit_1),
                    logLik_fit2 = logLik(fit_2), AIC_fit2 = AIC(fit_2))]
    }))


usethis::use_data(fits_TrueModels_t5, overwrite = TRUE)
