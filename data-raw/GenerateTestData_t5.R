library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(phytools)
library(ggtree)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

options(PCMBase.Threshold.EV = 1e-7)

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



set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")

# number of regimes
R <- 2
# number of traits
k <- 2

# this will generate a non-ultrametric tree with 318 tips
treeFossil318 <- pbtree(n=200, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil318)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil318$edge.length <- treeFossil318$edge.length / max(PCMTreeNodeTimes(treeFossil318)) * 166.2
treeFossil318 <- PCMTree(treeFossil318)
PCMTreeSetPartition(treeFossil318)

treeExtant318 <- pbtree(n=318, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant318)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant318$edge.length <- treeExtant318$edge.length / max(PCMTreeNodeTimes(treeExtant318)) * 166.2
treeExtant318 <- PCMTree(treeExtant318)
PCMTreeSetPartition(treeExtant318)

# this will generate a non-ultrametric tree with 638 tips
treeFossil638 <- pbtree(n=374, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil638)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil638$edge.length <- treeFossil638$edge.length / max(PCMTreeNodeTimes(treeFossil638)) * 166.2
treeFossil638 <- PCMTree(treeFossil638)
PCMTreeSetPartition(treeFossil638)


treeExtant638 <- pbtree(n=638, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant638)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant638$edge.length <- treeExtant638$edge.length / max(PCMTreeNodeTimes(treeExtant638)) * 166.2
treeExtant638 <- PCMTree(treeExtant638)
PCMTreeSetPartition(treeExtant638)


# PCMTreePlot(treeFossil318, layout="fan") + geom_nodelab()
# PCMTreePlot(treeExtant318, layout="fan") + geom_nodelab()
# PCMTreePlot(treeFossil638, layout="fan") + geom_nodelab(size = 3)
# PCMTreePlot(treeExtant638, layout="fan") + geom_nodelab(size = 3)

set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")

# this will generate a non-ultrametric tree with 80 tips
treeFossil80 <- pbtree(n=48, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil80)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil80$edge.length <- treeFossil80$edge.length / max(PCMTreeNodeTimes(treeFossil80)) * 166.2
treeFossil80 <- PCMTree(treeFossil80)
PCMTreeSetPartition(treeFossil80)

set.seed(2, kind = "Mersenne-Twister", normal.kind = "Inversion")
treeExtant80 <- pbtree(n=80, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant80)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant80$edge.length <- treeExtant80$edge.length / max(PCMTreeNodeTimes(treeExtant80)) * 166.2
treeExtant80 <- PCMTree(treeExtant80)
PCMTreeSetPartition(treeExtant80)

# this will generate a non-ultrametric tree with 159 tips
treeFossil159 <- pbtree(n=106, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil159)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossil159$edge.length <- treeFossil159$edge.length / max(PCMTreeNodeTimes(treeFossil159)) * 166.2
treeFossil159 <- PCMTree(treeFossil159)
PCMTreeSetPartition(treeFossil159)

treeExtant159 <- pbtree(n=159, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant159)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtant159$edge.length <- treeExtant159$edge.length / max(PCMTreeNodeTimes(treeExtant159)) * 166.2
treeExtant159 <- PCMTree(treeExtant159)
PCMTreeSetPartition(treeExtant159)


# PCMTreePlot(treeFossil80, layout="fan") + geom_nodelab()
# PCMTreePlot(treeExtant80, layout="fan") + geom_nodelab()
# PCMTreePlot(treeFossil159, layout="fan") + geom_nodelab(size = 3)
# PCMTreePlot(treeExtant159, layout="fan") + geom_nodelab(size = 3)



# number of different modelMappings per per clustering
nRandomModelMappingsPerClustering <- 4

# number of generated random models per model-mapping  of a tree
nRandomParamsPerMapping <- 8

# number of simulations per random model
nSimulationsPerRandomParam <- 4

# number of traits
k <- 2

set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")

testData_t5 <- rbindlist(
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


testData_t5[, treeWithRegimes:=lapply(1:.N, function(i) {
  PCMTreeSetPartition(tree[[i]], nodes = clusterNodes[[i]], inplace = FALSE)
})]

# this will reorder the nodes in clusterNodes column according to the preorder traversal of
# the tree
testData_t5[, clusterNodes:=lapply(1:.N, function(i) {
  PCMTreeGetLabels(treeWithRegimes[[i]])[PCMTreeGetPartition(treeWithRegimes[[i]])]
})]

# see how the treees with regimes look like
#PCMTreePlot(testData_t5$treeWithRegimes[[5]], layout="fan") + geom_nodelab()

# replicate each row in testData_t5 nRandomParamsPerMapping times
testData_t5 <- testData_t5[rep(1:.N, each = nRandomParamsPerMapping)]


# generate random models
testData_t5[, model:=lapply(1:.N, function(i) {
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

# replicate each row in testData_t5 nSimulationsPerRandomParam times
testData_t5 <- testData_t5[sort(rep_len(1:.N, length.out = .N * nSimulationsPerRandomParam))]

# simulate the trait values
testData_t5[, X:=lapply(1:.N, function(i) {
  cat(i, ", ")
  res <- try(PCMSim(treeWithRegimes[[i]], model[[i]], X0 = model[[i]]$X0), silent = TRUE)
  if(class(res) == "try-error") {
    print(res)
    NULL
  } else {
    res
  }
})]

# The TWO models for rows 1997:2000 and 2025:2028 had defective H matrices for
# some of the OU regimes. Both of these models were with MGPMSimulations.Theta.factor = 1.0
# Here, we regenerate these models and the corresponding data

for(i in c(1997, 2025)) {
  testData_t5[i:(i+nSimulationsPerRandomParam-1), model:=rep({
    model <- do.call(
      MixedGaussian,
      c(list(k = k, modelTypes = simulatedModels, mapping = mapping[[1L]]),
        argsMixedGaussian_SimulatedModels))

    options(MGPMSimulations.Theta.factor = 1.0)

    vecParams <- PCMParamRandomVecParams(model, n = 1)
    vecParams <- round(vecParams, digits = 2)

    PCMParamLoadOrStore(model, vecParams, 0, load = TRUE)
    # use a fixed value for X0 in all simulations
    model$X0[] <- c(1.0, -1.0)
    list(model)
  }, nSimulationsPerRandomParam)]
}

# simulate the trait values
testData_t5[sapply(X, is.null), X:=lapply(1:.N, function(i) {
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
testData_t5[, logLik:=lapply(1:.N, function(i) {

  lik <- PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]],
                metaI = PCMInfoCpp(X[[i]], treeWithRegimes[[i]], model[[i]]))
  cat(i, ": ", round(lik, 2), "\n")
  lik
})]

# AIC calculated at the model used to generate the data
testData_t5[, AIC:=lapply(1:.N, function(i) {
  model <- model[[i]]
  attr(model, "tree") <- treeWithRegimes[[i]]
  attr(model, "X") <- X[[i]][, 1:PCMTreeNumTips(treeWithRegimes[[i]])]
  attr(model, "SE") <- matrix(0.0, 2, PCMTreeNumTips(treeWithRegimes[[i]]))
  attr(model, "PCMInfoFun") <- PCMInfoCpp
  a <- AIC(model)
  cat(i, ": ", round(a, 2), "\n")
  a
})]

testData_t5[, nobs:=sapply(treeWithRegimes, PCMTreeNumTips)]
testData_t5[, df:=sapply(model, PCMParamCount, TRUE,TRUE)]

testData_t5[, IdMappingForClustering:=rep(
  rep(1:nRandomModelMappingsPerClustering,
      each=nSimulationsPerRandomParam*nRandomParamsPerMapping),
  .N/nRandomModelMappingsPerClustering/nRandomParamsPerMapping/nSimulationsPerRandomParam)]

testData_t5[, IdParamForMapping:=rep(
  rep(1:nRandomParamsPerMapping,
      each=nSimulationsPerRandomParam),
  .N/nRandomParamsPerMapping/nSimulationsPerRandomParam)]

testData_t5[, IdSimulationForParam:=rep(1:nSimulationsPerRandomParam,
                                         .N/nSimulationsPerRandomParam)]

testData_t5[, IdTree:=as.integer(max(.I)/256), by=list(nobs, treeType)]
testData_t5[, IdClusteringForTree:=rep(1:2, each=128), by=IdTree]

testData_t5[, IdGlob := .I]

setkey(
  testData_t5,
  IdTree,
  IdClusteringForTree,
  IdMappingForClustering,
  IdParamForMapping,
  IdSimulationForParam,
  IdGlob)

# store the data within the package TestPCMFit
usethis::use_data(testData_t5, overwrite = TRUE)


# Selected simulations for fitting best-clade-2 algorithm

# N = 8     0
# Tree−type ultrametric / N=80 / 2 regimes  / Mapping 1. BC             1:32
# Tree−type ultrametric / N=80 / 2 regimes  / Mapping 2. DF             33:64
# Tree−type ultrametric / N=80 / 3 regimes  / Mapping 1. DAB            129:160
# Tree−type ultrametric / N=80 / 3 regimes  / Mapping 2. BEC            161:192

# Tree−type non−ultrametric / N=80 / 2 regimes  / Mapping 1. FB         257:288
# Tree−type non−ultrametric / N=80 / 2 regimes  / Mapping 2. DA         289:320
# Tree−type non−ultrametric / N=80 / 3 regimes  / Mapping 1. FCC        385:416
# Tree−type non−ultrametric / N=80 / 3 regimes  / Mapping 2. DCB        417:448

# N = 159
# Tree−type ultrametric / N=159 / 2 regimes  / Mapping 1. ED            513:544
# Tree−type ultrametric / N=159 / 2 regimes  / Mapping 4. AC            609:640
# Tree−type non−ultrametric / N=159 / 2 regimes  / Mapping 1. AF        769:800
# Tree−type non−ultrametric / N=159 / 2 regimes  / Mapping 2. CF        801:832

# Tree−type ultrametric / N=159 / 5 regimes  / Mapping 1. EECFC         641:672
# Tree−type ultrametric / N=159 / 5 regimes  / Mapping 3. DCFBC         705:736
# Tree−type non−ultrametric / N=159 / 5 regimes  / Mapping 1. FCEFC     897:928
# Tree−type non−ultrametric / N=159 / 5 regimes  / Mapping 4. ADFEE     993:1024

# N = 318
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 3. DC            1089:1120
# Tree−type ultrametric / N=318 / 2 regimes  / Mapping 4. BF            1121:1152
# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 1. DD        1281:1312
# Tree−type non−ultrametric / N=318 / 2 regimes  / Mapping 3. ED        1345:1376

# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 1. DBACFDFE      1153:1184
# Tree−type ultrametric / N=318 / 8 regimes  / Mapping 2. CCAAEACD      1185:1216
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 1. ECBEAFDD  1409:1440
# Tree−type non−ultrametric / N=318 / 8 regimes  / Mapping 3. BFCEFCAC  1473:1504

# N = 638
# Tree−type ultrametric / N=638 / 2 regimes  / Mapping 1. DE            1537:1568
# Tree−type ultrametric / N=638 / 2 regimes  / Mapping 2. DF            1569:1600
# Tree−type ultrametric / N=638 / 8 regimes  / Mapping 1. FBEEFDEC      1665:1696
# Tree−type ultrametric / N=638 / 8 regimes  / Mapping 2. AFBDAFBE      1697:1728

# Tree−type non−ultrametric / N=638 / 2 regimes  / Mapping 2. FC        1825:1856
# Tree−type non−ultrametric / N=638 / 2 regimes  / Mapping 4. BD        1889:1920
# Tree−type non−ultrametric / N=638 / 8 regimes  / Mapping 1. FDBACFCA  1921:1952
# Tree−type non−ultrametric / N=638 / 8 regimes  / Mapping 4. CFEFCDCA  2017:2048


testData_t5_fittedIds <- c(
  1:32,
  33:64,
  129:160,
  161:192,
  257:288,
  289:320,
  385:416,
  417:448,
  513:544,
  609:640,
  769:800,
  801:832,
  641:672,
  705:736,
  897:928,
  993:1024,
  1089:1120,
  1121:1152,
  1281:1312,
  1345:1376,
  1153:1184,
  1185:1216,
  1409:1440,
  1473:1504,
  1537:1568,
  1569:1600,
  1665:1696,
  1697:1728,
  1825:1856,
  1889:1920,
  1921:1952,
  2017:2048)
usethis::use_data(testData_t5_fittedIds, overwrite = TRUE)
