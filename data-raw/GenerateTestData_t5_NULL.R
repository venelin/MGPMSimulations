# Generate trait data for k = 2, 4, 6, 8 traits with a single BM_A or BM_B regime
# on the ultrametric and non-ultrametric trees testData_t5.rda.
# This has the goal of assessing how often a false shift and/or OU type model is
# identified when the true data does not contain any shift or OU model type.

library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(MGPMSimulations)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

options(PCMBase.Threshold.EV = 1e-4)
options(PCMBase.Threshold.SV = 1e-5)
options(PCMBase.Threshold.Lambda_ij = 1e-4)

simulatedModels <- MGPMDefaultModelTypes()
argsMixedGaussian_SimulatedModels <- Args_MixedGaussian_MGPMDefaultModelTypes()


# We overwrite the limits for the model parameters. This
# prevents generating data for which some of the transition covariance matrices
# are singular due to too extreme random parameters. Note that these limits have
# been chosent to be adequate with the time-scale of the trees (all trees are of
# depth 166.2)

# these options will affect the H matrix, so that we don't have to specify it
# explicitly in the functions below
options(PCMBase.ParamValue.LowerLimit = -2,
        PCMBase.ParamValue.UpperLimit = 2,
        PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = .01)

PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[matDiag] <- .05
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- .0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[, , r][matDiag] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[, , r][matUppTri] <- .0
      }
    }
  }
  o
}

PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[matDiag] <- .5
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- .2
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[, , r][matDiag] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[, , r][matUppTri] <- .2
      }
    }
  }
  o
}



# number of different modelMappings per per clustering
nRandomModelMappingsPerClustering <- 2L

# For NULL tests, we vary the number of traits: 2, 4, 6, 8
nNumTraits <- 4

# number of generated random models per model-mapping  of a tree
nRandomParamsPerMapping <- 4

# number of simulations per random model
nSimulationsPerRandomParam <- 2

testData_t5_NULL <- testData_t5[
  numClusters == 2,
  list(
    tree = tree[1L],
    numClusters = 1L,
    clusterNodes = list(as.character(PCMTreeNumTips(tree[[1L]]) + 1L)),
    mapping = as.list(
      c(rep(1, nNumTraits * nRandomParamsPerMapping),
        rep(2, nNumTraits * nRandomParamsPerMapping))),
    numTraits = rep(
      c(2L, 4L, 6L, 8L),
      nRandomModelMappingsPerClustering * nRandomParamsPerMapping)),
  by = list(treeType, treeSize, nobs, IdTree)]

testData_t5_NULL[, treeWithRegimes := lapply(1:.N, function(i) {
  PCMTreeSetPartition(PCMTree(tree[[i]]), nodes = clusterNodes[[i]], inplace = FALSE)
})]

set.seed(6, kind = "Mersenne-Twister", normal.kind = "Inversion")

# generate random models
testData_t5_NULL[, model:=lapply(1:.N, function(i) {
  model <- do.call(
    MixedGaussian,
    c(list(k = numTraits[[i]], modelTypes = simulatedModels, mapping = mapping[[i]]),
      argsMixedGaussian_SimulatedModels))

  cat("i=", i, ", k=", numTraits[[i]], "N=", nobs[[i]], ", treeType=", treeType[[i]], ", mapping=", mapping[[i]], ": ")

  while(TRUE) {
    vecParams <- PCMParamRandomVecParams(model, n = 1)
    vecParams <- round(vecParams, digits = 4)


    PCMParamLoadOrStore(model, vecParams, 0, load = TRUE)
    # use a fixed value for X0 in all simulations
    model$X0[] <- rep(c(1.0, -1.0), numTraits[[i]]/2)

    Xtest <- try(PCMSim(treeWithRegimes[[i]], model, X0 = model$X0), silent = TRUE)
    if(class(Xtest) == "try-error") {
      next
    }
    lik <- PCMLik(
      Xtest, treeWithRegimes[[i]], model,
      metaI = PCMInfoCpp(Xtest, treeWithRegimes[[i]], model))
    if(lik > getOption("PCMBase.Value.NA")) {
      cat("\n")
      break
    } else {
      cat(".")
    }
  }

  model
})]

# replicate each row in testData_t5_NULL nSimulationsPerRandomParam times
testData_t5_NULL <- testData_t5_NULL[
  sort(rep_len(1:.N, length.out = .N * nSimulationsPerRandomParam))]

# simulate the trait values
testData_t5_NULL[, X:=lapply(1:.N, function(i) {
  cat(i, ", ")
  res <- try(PCMSim(treeWithRegimes[[i]], model[[i]], X0 = model[[i]]$X0), silent = TRUE)
  if(class(res) == "try-error") {
    print(res)
    NULL
  } else {
    res
  }
})]

# log-likelihood calculated at the model used to generate the data
testData_t5_NULL[, logLik:=lapply(1:.N, function(i) {

  lik <- PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]],
                metaI = PCMInfoCpp(X[[i]], treeWithRegimes[[i]], model[[i]]))
  cat(i, ": ", round(lik, 2), "\n")
  lik
})]

# AIC calculated at the model used to generate the data
testData_t5_NULL[, AIC:=lapply(1:.N, function(i) {
  model <- model[[i]]
  attr(model, "tree") <- treeWithRegimes[[i]]
  attr(model, "X") <- X[[i]][, 1:PCMTreeNumTips(treeWithRegimes[[i]])]
  attr(model, "SE") <- matrix(0.0, numTraits[[i]], PCMTreeNumTips(treeWithRegimes[[i]]))
  attr(model, "PCMInfoFun") <- PCMInfoCpp
  a <- AIC(model)
  cat(i, ": ", round(a, 2), "\n")
  a
})]

testData_t5_NULL[, nobs:=sapply(treeWithRegimes, PCMTreeNumTips)]
testData_t5_NULL[, df:=sapply(model, PCMParamCount, TRUE,TRUE)]

testData_t5_NULL[, IdMappingForClustering := rep(
  rep(1:nRandomModelMappingsPerClustering,
      each=nNumTraits * nRandomParamsPerMapping * nSimulationsPerRandomParam),
  .N / nRandomModelMappingsPerClustering /
    nNumTraits /
    nRandomParamsPerMapping /
    nSimulationsPerRandomParam)]

testData_t5_NULL[, IdParamForMapping:=rep(
  rep(1:nRandomParamsPerMapping,
      each = nNumTraits * nSimulationsPerRandomParam),
  .N / nRandomParamsPerMapping / nSimulationsPerRandomParam / nNumTraits)]

testData_t5_NULL[, IdSimulationForParam:=rep(1:nSimulationsPerRandomParam,
                                         .N/nSimulationsPerRandomParam)]
testData_t5_NULL[, IdClusteringForTree:=1L]

setkey(
  testData_t5_NULL,
  IdTree,
  IdClusteringForTree,
  IdMappingForClustering,
  numTraits,
  IdParamForMapping,
  IdSimulationForParam)

testData_t5_NULL[, IdGlob:=.I]

# store the data within the package MGPMSimulations
usethis::use_data(testData_t5_NULL, overwrite = TRUE)

testData_t5_NULL_fittedIds <- c(
  testData_t5_NULL[nobs %in% c(80) & IdParamForMapping %in% c(1, 3), IdGlob],
  testData_t5_NULL[nobs %in% c(318) & IdParamForMapping %in% c(1, 3), IdGlob])

usethis::use_data(testData_t5_NULL_fittedIds, overwrite = TRUE)
