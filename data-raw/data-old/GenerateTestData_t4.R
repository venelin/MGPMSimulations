library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(phytools)
library(ggtree)

PCMFit::GeneratePCMModelTypes()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

simulatedModels <- PCMFit::MGPMDefaultModelTypes()
argsMixedGaussian_SimulatedModels <- PCMFit::Args_MixedGaussian_MGPMDefaultModelTypes()


# We overwrite the limits for the model parameters. This
# prevents generating data for which some of the transition covariance matrices
# are singular due to too extreme random parameters. Note that these limits have
# been chosent to be adequate with the time-scale of the trees (all trees are of
# depth 166.2)

# these options will affect the H matrix, so that we don't have to specify it
# explicitly in the functions below
options(PCMBase.ParamValue.LowerLimit = -4,
        PCMBase.ParamValue.UpperLimit = 4,
        PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = .1)

PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[1, 1] <- .05
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[1, 1, r] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .0
      }
    }
  }
  o
}

PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .2
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .2
      }
    }
  }
  o
}

PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 3.0
    o$Theta[2] <- 2.0
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 3.0
      o$Theta[2, r] <- 2.0
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .05
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .0
      }
    }
  }
  o
}

PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 3.0 + 3.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
    o$Theta[2] <- 2.0 + 2.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 3.0 + 3.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
      o$Theta[2, r] <- 2.0 + 2.0 * getOption("MGPMSimulations.Theta.factor", 1.0)
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .2
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .2
      }
    }
  }
  o
}



set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# this will generate a non-ultrametric tree with 318 tips
treeFossil318 <- pbtree(n=200, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil318)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil318$edge.length <- treeFossil318$edge.length / max(PCMTreeNodeTimes(treeFossil318)) * 166.2
PCMTreeSetDefaultRegime(treeFossil318, regime = 1)

treeExtant318 <- pbtree(n=318, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant318)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant318$edge.length <- treeExtant318$edge.length / max(PCMTreeNodeTimes(treeExtant318)) * 166.2
PCMTreeSetDefaultRegime(treeExtant318, regime = 1)

# this will generate a non-ultrametric tree with 638 tips
treeFossil638 <- pbtree(n=374, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil638)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil638$edge.length <- treeFossil638$edge.length / max(PCMTreeNodeTimes(treeFossil638)) * 166.2
PCMTreeSetDefaultRegime(treeFossil638, regime = 1)

treeExtant638 <- pbtree(n=638, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant638)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant638$edge.length <- treeExtant638$edge.length / max(PCMTreeNodeTimes(treeExtant638)) * 166.2
PCMTreeSetDefaultRegime(treeExtant638, regime = 1)


# PCMTreePlot(treeFossil318, layout="fan") + geom_nodelab()
# PCMTreePlot(treeExtant318, layout="fan") + geom_nodelab()
# PCMTreePlot(treeFossil638, layout="fan") + geom_nodelab(size = 3)
# PCMTreePlot(treeExtant638, layout="fan") + geom_nodelab(size = 3)

set.seed(2)

# this will generate a non-ultrametric tree with 80 tips
treeFossil80 <- pbtree(n=48, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil80)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil80$edge.length <- treeFossil80$edge.length / max(PCMTreeNodeTimes(treeFossil80)) * 166.2
PCMTreeSetDefaultRegime(treeFossil80, regime = 1)

set.seed(2)
treeExtant80 <- pbtree(n=80, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant80)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant80$edge.length <- treeExtant80$edge.length / max(PCMTreeNodeTimes(treeExtant80)) * 166.2
PCMTreeSetDefaultRegime(treeExtant80, regime = 1)

# this will generate a non-ultrametric tree with 159 tips
treeFossil159 <- pbtree(n=106, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil159)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil159$edge.length <- treeFossil159$edge.length / max(PCMTreeNodeTimes(treeFossil159)) * 166.2
PCMTreeSetDefaultRegime(treeFossil159, regime = 1)

treeExtant159 <- pbtree(n=159, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant159)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant159$edge.length <- treeExtant159$edge.length / max(PCMTreeNodeTimes(treeExtant159)) * 166.2
PCMTreeSetDefaultRegime(treeExtant159, regime = 1)


# PCMTreePlot(treeFossil80, layout="fan") + geom_nodelab()
# PCMTreePlot(treeExtant80, layout="fan") + geom_nodelab()
# PCMTreePlot(treeFossil159, layout="fan") + geom_nodelab(size = 3)
# PCMTreePlot(treeExtant159, layout="fan") + geom_nodelab(size = 3)



# number of different modelMappings per per clustering
nRandomModelMappingsPerClustering <- 4

# number of generated random models per model-mapping  of a tree
nRandomParamsPerMapping <- 4

# number of simulations per random model
nSimulationsPerRandomParam <- 2

# number of traits
k <- 2

set.seed(1)

testData_t4 <- rbindlist(
  list(
    data.table(
      treeType = "ultrametric",
      treeSize = "N=80",
      tree = list(treeExtant80),
      numClusters = 2,
      clusterNodes = list(as.character(c(81, 122))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "ultrametric",
      treeSize= "N=80",
      tree = list(treeExtant80),
      numClusters = 3,
      clusterNodes = list(as.character(c(81, 101, 122))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 3, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=80",
      tree = list(treeFossil80),
      numClusters = 2,
      clusterNodes = list(as.character(c(81, 124))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=80",
      tree = list(treeFossil80),
      numClusters = 3,
      clusterNodes = list(as.character(c(81, 105, 125))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 3, replace = TRUE))
    ),

    data.table(
      treeType = "ultrametric",
      treeSize= "N=159",
      tree = list(treeExtant159),
      numClusters = 2,
      clusterNodes = list(as.character(c(160, 205))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "ultrametric",
      treeSize= "N=159",
      tree = list(treeExtant159),
      numClusters = 5,
      clusterNodes = list(as.character(c(160, 185, 205, 269, 297))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 5, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=159",
      tree = list(treeFossil159),
      numClusters = 2,
      clusterNodes = list(as.character(c(160, 166))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=159",
      tree = list(treeFossil159),
      numClusters = 5,
      clusterNodes = list(as.character(c(160, 161, 168, 208, 239))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 5, replace = TRUE))
    ),


    data.table(
      treeType = "ultrametric",
      treeSize = "N=318",
      tree = list(treeExtant318),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 499))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "ultrametric",
      treeSize= "N=318",
      tree = list(treeExtant318),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 499, 438, 360,
                                         320, 486, 376, 583))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=318",
      tree = list(treeFossil318),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 393))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=318",
      tree = list(treeFossil318),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 393, 430, 484,
                                         517, 486, 600, 343))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),

    data.table(
      treeType = "ultrametric",
      treeSize= "N=638",
      tree = list(treeExtant638),
      numClusters = 2,
      clusterNodes = list(as.character(c(639, 1007))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "ultrametric",
      treeSize= "N=638",
      tree = list(treeExtant638),
      numClusters = 8,
      clusterNodes = list(as.character(c(639, 1105, 1007, 817,
                                         867, 700, 1236, 1177))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=638",
      tree = list(treeFossil638),
      numClusters = 2,
      clusterNodes = list(as.character(c(639, 883))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "non-ultrametric",
      treeSize= "N=638",
      tree = list(treeFossil638),
      numClusters = 8,
      clusterNodes = list(as.character(c(639, 685, 644, 641,
                                         914, 971, 1136, 975))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    )

    ))


testData_t4[, treeWithRegimes:=lapply(1:.N, function(i) {
  PCMTreeSetRegimes(tree[[i]], nodes = clusterNodes[[i]], inplace = FALSE)
})]

# this will reorder the nodes in clusterNodes column according to the preorder traversal of
# the tree
testData_t4[, clusterNodes:=lapply(1:.N, function(i) {
  PCMTreeGetLabels(treeWithRegimes[[i]])[PCMTreeGetStartingNodesRegimes(treeWithRegimes[[i]])]
})]

# see how the treees with regimes look like
#PCMTreePlot(testData_t4$treeWithRegimes[[5]], layout="fan") + geom_nodelab()

# replicate each row in testData_t4 nRandomParamsPerMapping times
testData_t4 <- testData_t4[rep(1:.N, each = nRandomParamsPerMapping)]


# generate random models
testData_t4[, model:=lapply(1:.N, function(i) {
  model <- do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = simulatedModels, mapping = mapping[[i]]),
      argsMixedGaussian_SimulatedModels))

  if((i %% nRandomParamsPerMapping) %in% seq_len(nRandomParamsPerMapping / 2)) {
    options(MGPMSimulations.Theta.factor = 1.0)
  } else {
    options(MGPMSimulations.Theta.factor = 0.1)
  }

  vecParams <- PCMParamRandomVecParams(model, n = 1)
  vecParams <- round(vecParams, digits = 2)

  PCMParamLoadOrStore(model, vecParams, 0, load = TRUE)
  # use a fixed value for X0 in all simulations
  model$X0[] <- c(1.0, -1.0)
  model
})]

# replicate each row in testData_t4 nSimulationsPerRandomParam times
testData_t4 <- testData_t4[sort(rep_len(1:.N, length.out = .N * nSimulationsPerRandomParam))]

# simulate the trait values
testData_t4[, X:=lapply(1:.N, function(i) {
  cat(i, ", ")
  PCMSim(treeWithRegimes[[i]], model[[i]], X0 = model[[i]]$X0)
})]

# for a few simulations there was an eigenvalue below the threshold of 1e-5 for the
# variance covariance transition matrix on 1 node. To calculate the likelihood
# and AIC for these 10 simulations we set the threshold to a more tolerant value.
# In the inference, the threshold will be left as per defaul (i.e. 1e-5).
options(PCMBase.Threshold.EV = 1e-7)

# log-likelihood calculated at the model used to generate the data
testData_t4[, logLik:=lapply(1:.N, function(i) {

  lik <- PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]],
                metaI = PCMInfoCpp(X[[i]], treeWithRegimes[[i]], model[[i]]))
  cat(i, ": ", round(lik, 2), "\n")
  lik
})]

# AIC calculated at the model used to generate the data
testData_t4[, AIC:=lapply(1:.N, function(i) {
  model <- model[[i]]
  attr(model, "tree") <- treeWithRegimes[[i]]
  attr(model, "X") <- X[[i]][, 1:PCMTreeNumTips(treeWithRegimes[[i]])]
  attr(model, "SE") <- matrix(0.0, 2, PCMTreeNumTips(treeWithRegimes[[i]]))
  attr(model, "PCMInfoFun") <- PCMInfoCpp
  a <- AIC(model)
  cat(i, ": ", round(a, 2), "\n")
  a
})]

testData_t4[, nobs:=sapply(treeWithRegimes, PCMTreeNumTips)]
testData_t4[, df:=sapply(model, PCMParamCount, TRUE,TRUE)]

testData_t4[, IdMappingForClustering:=rep(
  rep(1:nRandomModelMappingsPerClustering,
      each=nSimulationsPerRandomParam*nRandomParamsPerMapping),
  .N/nRandomModelMappingsPerClustering/nRandomParamsPerMapping/nSimulationsPerRandomParam)]

testData_t4[, IdParamForMapping:=rep(
  rep(1:nRandomParamsPerMapping,
      each=nSimulationsPerRandomParam),
  .N/nRandomParamsPerMapping/nSimulationsPerRandomParam)]

testData_t4[, IdSimulationForParam:=rep(1:nSimulationsPerRandomParam,
                                         .N/nSimulationsPerRandomParam)]

# store the data within the package TestPCMFit
usethis::use_data(testData_t4, overwrite = TRUE)
