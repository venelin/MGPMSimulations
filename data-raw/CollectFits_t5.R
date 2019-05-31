# example calls:

# R --vanilla -f CollectFits_t5.R --args fits_MGPM_A_F_all_AIC_t5 Results_t5_MGPM_A_F_all_NXXX MGPM_A_F_all_id_

# R --vanilla -f CollectFits_t5.R --args fits_MGPM_A_F_all_AIC2_t5 Results_t5_MGPM_A_F_all_AIC2_NXXX MGPM_A_F_all_AIC2_id_

# R --vanilla -f CollectFits_t5.R --args fits_MGPM_A_F_best_clade_2_AIC_t5 Results_t5_MGPM_A_F_best_clade_2_NXXX MGPM_A_F_best_clade_2_id_

# R --vanilla -f CollectFits_t5.R --args fits_MGPM_A_F_best_clade_2_RR_AIC_t5 Results_t5_MGPM_A_F_best_clade_2_RR_NXXX MGPM_A_F_best_clade_2_RR_id_

# R --vanilla -f CollectFits_t5.R --args fits_MGPM_A_F_best_clade_RR_AIC_t5 Results_t5_MGPM_A_F_best_clade_RR_NXXX MGPM_A_F_best_clade_RR_id_


# SURFACE fits
# R --vanilla -f CollectFits_t5.R --args fits_SURFACE_best_clade_2_AICc_t5 Results_t5_NXXX_SURFACE SURFACE_best_clade_2_id_
# 254 of the above had inferred R=1

# R --vanilla -f CollectFits_t5.R --args fits_SURFACE_best_clade_2_AICc_mcs10_t5 Results_t5_NXXX_SURFACE_mcs10 SURFACE_best_clade_2_mcs10_id_
# 63 of the above had R = 1

# SCALAR OU fits
# R --vanilla -f CollectFits_t5.R --args fits_SCALAROU_best_clade_2_AIC_t5 Results_t5_NXXX_SCALAROU SCALAROU_best_clade_2_id_
# 7 of the above had inferred R=1


library(PCMFit)
library(PCMBase)
library(MGPMSimulations)
library(data.table)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
  # e.g. fits_MGPM_A_F_best_clade_2_AIC_t5
  dataName <- as.character(args[[1]])
  # e.g. Results_t5_MGPM_A_F_best_clade_2_
  resultDirPrefix <- as.character(args[[2]])
  # e.g. MGPM_A_F_best_clade_2_id_
  resultFilePrefix <- as.character(args[[3]])
} else {
  stop("data.table name is missing.")
}

data <- rbindlist(
  lapply(
    testData_t5_fittedIds,
    function(id) {
      if(id %% 10L == 0) {
        cat("id:", id, "\n")
      }

      prefixFiles =  paste0(resultFilePrefix, id)

      NSize <- paste0("N", testData_t5[id, nobs])

      file <- paste0(

        gsub("NXXX", NSize, resultDirPrefix), "/", prefixFiles,
        "/FinalResult_", prefixFiles, ".RData")

      if(file.exists(file)) {
        # this should load an object called fitMappings
        load(file)
        fit <- RetrieveBestFitScore(fitMappings)
        ll <- attr(fit$inferredModel, "ll")
        sc <- attr(fit$inferredModel, "score")

        attr(fit$inferredModel, "X") <-
          attr(fit$inferredModel, "SE") <-
          attr(fit$inferredModel, "tree") <-
          attr(fit$inferredModel, "ll") <-
          attr(fit$inferredModel, "score") <- NULL

        testData_t5[id,
                    list(
                      IdGlob, IdTree, IdClusteringForTree,
                      IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                      clusterNodes = list(fit$inferredRegimeNodes),
                      mapping = list(fit$inferredMappingIdx),
                      model = list(fit$inferredModel),
                      logLik = ll,
                      score = sc)]
      } else {
        testData_t5[id,
                    list(
                      IdGlob, IdTree, IdClusteringForTree,
                      IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                      clusterNodes = list(NULL), mapping = list(NULL),
                      model = list(NULL),
                      logLik = as.double(NA), score = as.double(NA))]
      }

    }))

eval(parse(text = paste0(dataName, " <- data")))
eval(parse(text = paste0("usethis::use_data(", dataName, ", overwrite = TRUE)")))
