# The name of the fits_* data.table has to be specified
# as a first argument of the R command executing the script. As a second argument
# a 0 or 1 is specified to indicate whether MGPM model criteria should be calculated,
# e.g. is H a symmetric matrix, are the traits evolving in a correlated fashion, etc.

# examples:
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_TrueModels_t5 1
#
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_MGPM_A_F_all_AIC_t5 1
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_MGPM_A_F_all_AIC2_t5 1
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_MGPM_A_F_best_clade_2_AIC_t5 1
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_MGPM_A_F_best_clade_2_RR_AIC_t5 1
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_MGPM_A_F_best_clade_RR_AIC_t5 1
#
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_SURFACE_best_clade_2_AICc_t5 0
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_SURFACE_best_clade_2_AICc_mcs10_t5 0
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_surface_fwd_t5 0
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_surface_bwd_t5 0
#
# bsub -n 1 R --vanilla -f CalculateBinaryCriteria_t5.R --args fits_SCALAROU_best_clade_2_AIC_t5 0

library(PCMBase)
library(MGPMSimulations)
library(data.table)
library(fpc)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 2) {
  dataName <- as.character(args[1])
  isMGPMFit <- as.logical(as.integer(args[2]))
} else {
  stop("data.table name and/or  is missing.")
}

# short name shortcuts:
tpr <- function(a, b) {
  TruePositiveRate(a, b)
}
fpr <- function(a, b) {
  FalsePositiveRate(a, b)
}

data <- mget(
  dataName, envir=as.environment(-1), ifnotfound=list(NULL), inherits=TRUE)[[1]]

if(is.null(data)) {
  stop(paste0("data.table ", dataName, " was not found."))
}

# testResultData[, crit2:=factor(crit, levels = c("Cluster",
#                                                 "BM process",
#                                                 "OU process",
#                                                 "Uncorrelated traits",
#                                                 "Correlated traits",
#                                                 "NonDiagonal H",
#                                                 "Asymmetric H"))]

data[, c("perf_Cluster_tpr",
         "perf_Cluster_fpr",
         "perf_BM_tpr",
         "perf_BM_fpr",
         "perf_OU_tpr",
         "perf_OU_fpr",
         "perf_Uncorrelated_tpr",
         "perf_Uncorrelated_fpr",
         "perf_Correlated_tpr",
         "perf_Correlated_fpr",
         "perf_NonDiagonalH_tpr",
         "perf_NonDiagonalH_fpr",
         "perf_AsymmetricH_tpr",
         "perf_AsymmetricH_fpr") := {

  cols <- lapply(.I, function(i) {
    if(i %% 10L == 0) {
      cat("i:", i, "\n")
    }

    res <- list(
      perf_Cluster = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_BM = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_OU = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_Uncorrelated = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_Correlated = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_NonDiagonalH = c(tpr = as.double(NA), fpr = as.double(NA)),
      perf_AsymmetricH = c(tpr = as.double(NA), fpr = as.double(NA))
    )

    if( is.PCM(model[[i]]) ) {

      idGlob <- IdGlob[i]

      trueTree <- testData_t5$treeWithRegimes[[idGlob]]
      trueMappedModels <- testData_t5$mapping[[idGlob]][trueTree$edge.regime]

      predictedTree <- PCMTree(testData_t5$tree[[idGlob]])

      if(is.null(names(clusterNodes[[i]]))) {
        # normal MixedGaussian model without lumped models, regimes are named
        # 1,2,3....
        PCMTreeSetPartition(
          predictedTree, nodes = clusterNodes[[i]])
      } else {
        # a MixedGaussian model formed by parsing a fit from the surface
        # package. Regimes are named a,b,c...
        PCMTreeSetPartRegimes(
          predictedTree, part.regime = clusterNodes[[i]], setPartition = TRUE)
      }


      matTRUE_Cluster <- PCMTreeMatrixNodesInSameRegime(trueTree)
      mode(matTRUE_Cluster) <- "double"
      matPredicted_Cluster <- PCMTreeMatrixNodesInSameRegime(predictedTree)
      mode(matPredicted_Cluster) <- "double"


      res$perf_Cluster = c(
          tpr = tpr(matPredicted_Cluster, matTRUE_Cluster),
          fpr = fpr(matPredicted_Cluster, matTRUE_Cluster))

      if(isMGPMFit) {
        predictedMappedModels <- mapping[[i]][
          PCMTreeGetPartRegimes(predictedTree)[
            PCMTreeGetPartsForNodes(predictedTree, predictedTree$edge[,2])]]

        res$perf_BM = c(
          tpr = tpr(is.BM(predictedMappedModels),
                    is.BM(trueMappedModels)),
          fpr = fpr(is.BM(predictedMappedModels),
                    is.BM(trueMappedModels)))

        res$perf_OU = c(
          tpr = tpr(is.OU(predictedMappedModels),
                    is.OU(trueMappedModels)),
          fpr = fpr(is.OU(predictedMappedModels),
                    is.OU(trueMappedModels)))

        res$perf_Uncorrelated = c(
          tpr = tpr(is.Uncorrelated(predictedMappedModels),
                    is.Uncorrelated(trueMappedModels)),
          fpr = fpr(is.Uncorrelated(predictedMappedModels),
                    is.Uncorrelated(trueMappedModels)))

        res$perf_Correlated = c(
          tpr = tpr(is.Correlated(predictedMappedModels),
                    is.Correlated(trueMappedModels)),
          fpr = fpr(is.Correlated(predictedMappedModels),
                    is.Correlated(trueMappedModels)))

        res$perf_NonDiagonalH = c(
          tpr = tpr(is.NonDiagonalH(predictedMappedModels),
                    is.NonDiagonalH(trueMappedModels)),
          fpr = fpr(is.NonDiagonalH(predictedMappedModels),
                    is.NonDiagonalH(trueMappedModels)))

        res$perf_AsymmetricH = c(
          tpr = tpr(is.NonSymmetricH(predictedMappedModels),
                    is.NonSymmetricH(trueMappedModels)),
          fpr = fpr(is.NonSymmetricH(predictedMappedModels),
                    is.NonSymmetricH(trueMappedModels)))
      } else {
        res$perf_BM = c(
          tpr = NA_real_,
          fpr = NA_real_)

        res$perf_OU = c(
          tpr = NA_real_,
          fpr = NA_real_)

        res$perf_Uncorrelated = c(
          tpr = NA_real_,
          fpr = NA_real_)

        res$perf_Correlated = c(
          tpr = NA_real_,
          fpr = NA_real_)

        res$perf_NonDiagonalH = c(
          tpr = NA_real_,
          fpr = NA_real_)

        res$perf_AsymmetricH = c(
          tpr = NA_real_,
          fpr = NA_real_)
      }
    }

    res
  })

  list(sapply(cols, function(.) .[[1]][1]),
       sapply(cols, function(.) .[[1]][2]),
       sapply(cols, function(.) .[[2]][1]),
       sapply(cols, function(.) .[[2]][2]),
       sapply(cols, function(.) .[[3]][1]),
       sapply(cols, function(.) .[[3]][2]),
       sapply(cols, function(.) .[[4]][1]),
       sapply(cols, function(.) .[[4]][2]),
       sapply(cols, function(.) .[[5]][1]),
       sapply(cols, function(.) .[[5]][2]),
       sapply(cols, function(.) .[[6]][1]),
       sapply(cols, function(.) .[[6]][2]),
       sapply(cols, function(.) .[[7]][1]),
       sapply(cols, function(.) .[[7]][2]))}]

eval(parse(text = paste0(dataName, " <- data")))
eval(parse(text = paste0("usethis::use_data(", dataName, ", overwrite = TRUE)")))
