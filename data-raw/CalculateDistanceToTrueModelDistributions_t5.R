# Calculates the mean vector mu and the variance covariance matrix Sigma
# for the (N x k) - variate Gaussian distribution at the tips of each tree, expected
# under a model of evolution described on each row of a given fits_* data.table
# in the MGPMSimulations package. Uses the functions stats::mahalanobis and
# fpc::bhattacharyya.dist to calclate the corresponding distances to the Gaussian
# distribution expected under the true (simulated) model for the corresponding
# data at the tips of the tree. Adds columns dMahalanobis and dBhattacharyya to
# the fits_* data.table. The name of the fits_* data.table has to be specified
# as a first argument of the R command executing the script.
# examples:
#
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_TrueModels_t5
#
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_MGPM_A_F_all_AIC_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_MGPM_A_F_all_AIC2_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_MGPM_A_F_best_clade_2_AIC_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_MGPM_A_F_best_clade_2_RR_AIC_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_MGPM_A_F_best_clade_RR_AIC_t5
#
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_SURFACE_best_clade_2_AICc_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_SURFACE_best_clade_2_AICc_mcs10_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_surface_fwd_t5
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_surface_bwd_t5
#
# bsub -n 1 R --vanilla -f CalculateDistanceToTrueModelDistributions_t5.R --args fits_SCALAROU_best_clade_2_AIC_t5

library(PCMBase)
library(MGPMSimulations)
library(data.table)
library(fpc)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)

load("trueModelDistributions_t5.RData")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  dataName <- as.character(args[1])
} else {
  stop("data.table name is missing.")
}

data <- mget(
  dataName, envir=as.environment(-1), ifnotfound=list(NULL), inherits=TRUE)[[1]]

if(is.null(data)) {
  stop(paste0("data.table ", dataName, " was not found."))
}

dBhattacharyya <- dMahalanobis <- NULL

data[, dBhattacharyya:=NULL]
data[, dMahalanobis:=NULL]

data[, c("dBhattacharyya", "dMahalanobis") := {
  cols12 <- lapply(.I, function(i) {
    if(i %% 10L == 0) {
      cat("i:", i, "\n")
    }

    if((is.null(dBhattacharyya) || is.na(dBhattacharyya[[i]])) &&
       is.PCM(model[[i]])) {

      idGlob <- IdGlob[i]

      trueMu <- trueModelDistributions_t5[
        testData_t5[
          idGlob,
          list(
            IdTree, IdClusteringForTree,
            IdMappingForClustering, IdParamForMapping)], mu[[1L]]]
      trueSigma <- trueModelDistributions_t5[
        testData_t5[
          idGlob,
          list(
            IdTree, IdClusteringForTree,
            IdMappingForClustering, IdParamForMapping)], Sigma[[1L]]]

      tree <- PCMTree(testData_t5$tree[[idGlob]])

      if(is.null(names(clusterNodes[[i]]))) {
        # normal MixedGaussian model without lumped models, regimes are named
        # 1,2,3....
        PCMTreeSetPartition(tree, nodes = clusterNodes[[i]])
      } else {
        # a MixedGaussian model formed by parsing a fit from the surface
        # package. Regimes are named a,b,c...
        PCMTreeSetPartRegimes(tree, part.regime = clusterNodes[[i]], setPartition = TRUE)
      }

      metaI <- PCMInfo(X = NULL, tree, model[[i]])


      mu <- try(as.vector(PCMMean(tree, model[[i]], metaI = metaI)), silent = TRUE)
      Sigma <- try(PCMVar(tree, model[[i]], metaI = metaI), silent = TRUE)

      if(class(mu) == "try-error" || class(Sigma) == "try-error") {
        print(mu)
        print(Sigma)
        cat("error calculating mu or Sigma for ", i, "; returning NA\n")
        c(as.double(NA), as.double(NA))
      } else {
        c(bhattacharyya.dist(mu, trueMu, Sigma, trueSigma),
          mahalanobis(mu, trueMu, trueSigma))
      }
    } else if( !is.null(dBhattacharyya) ) {
      c(as.double(dBhattacharyya[[i]]), as.double(dMahalanobis[[i]]))
    } else {
      c(as.double(NA), as.double(NA))
    }
  })
  list(sapply(cols12, function(.) .[[1]]), sapply(cols12, function(.) .[[2]]))
}]

eval(parse(text = paste0(dataName, " <- data")))
eval(parse(text = paste0("usethis::use_data(", dataName, ", overwrite = TRUE)")))
