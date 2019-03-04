# This script can be run locally or on a cluster using a command like:
#
# bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_MGPM_A_F_best_clade_RR.R --args 96 MGPM_A_F_best_clade_2_id_ ../../Results_t5_MGPM_A_F_best_clade_2_N638_24h
#
# bsub -M 30000 -n 20 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_MGPM_A_F_best_clade_RR.R --args $id MGPM_A_F_best_clade_2_id_ ../../Results_t5_MGPM_A_F_best_clade_2_N318
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package MGPMSimulations The blank XXX stays for the type of model to be inferred, e.g. ScalarOU or MixedGaussian
#
library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(MGPMSimulations)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
  prefixFilesPrev <- as.character(args[2])
  resDirPrev <- as.character(args[3])
} else {
  id <- 1
}

prefixFilesPrev = paste0(prefixFilesPrev, id)

prefixFiles = paste0("MGPM_A_F_best_clade_RR_id_", id)

if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a cluster)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE),
                                     outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  }
}

tableFitsPrev <- NULL
# try using a previously stored tableFits from a previous run that was interupted
filePrevResults <- paste0(
  resDirPrev, "/", prefixFilesPrev, "/FinalResult_", prefixFilesPrev, ".RData")

if(file.exists(filePrevResults)) {
  cat("Loading previously stored tableFits from file", filePrevResults, "...\n")

  # this should load a list object named listResults containg a named entry "tableFits":
  load(filePrevResults)

  tableFitsPrev <- fitMappings$tableFits

  setkey(tableFitsPrev, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)
}

# try using a previously stored tableFits from a previous run that was interupted
tableFitsCurrent <- NULL
fileCurrentResults <- paste0("CurrentResults_", prefixFiles, ".RData")
if(file.exists(fileCurrentResults)) {
  cat("Loading previously stored tableFits from file", fileCurrentResults, "...\n")

  # this should load a list object named listResults containg a named entry "tableFits":
  load(fileCurrentResults)

  tableFitsCurrent <- listResults$tableFits

  tempFiles <- list.files(pattern = paste0("^", prefixFiles, ".*.RData"))
  if(length(tempFiles) > 0) {
    cat("Loading previously stored tableFits from temporary files (", toString(tempFiles), ")...\n")
    tableFitsTempFiles <- rbindlist(
      lapply(tempFiles, function(file) {
        load(file)
        fits
      }))
    tableFitsCurrent <- PCMFit:::UpdateTableFits(tableFitsCurrent, tableFitsTempFiles)
  }

  setkey(tableFitsCurrent, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)
}

tableFits <- PCMFit:::UpdateTableFits(tableFitsPrev, tableFitsCurrent)

print(tableFits)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-7)

print(PCMOptions())

tree <- testData_t5$tree[[id]]
values <- testData_t5$X[[id]][, seq_len(PCMTreeNumTips(tree))]

# Specify F for each node as the only possible model during step 2. - not done,
# we prefer specifying best-clade instead.
# listAllowedModelTypesIndices = {
#   nodes <- PCMTreeListCladePartitions(tree, nNodes = 1L, minCladeSize = 20)
#   lstAllowedModelTypesInd <- as.list(rep(6L, length(nodes)))
#   names(lstAllowedModelTypesInd) <- as.character(unlist(nodes))
#   lstAllowedModelTypesInd[[as.character(PCMTreeNumTips(tree)+1)]] <- 6L
#   lstAllowedModelTypesInd
# }
#

fitMappings <- PCMFitMixed(

  values, tree,

  metaIFun = PCMInfoCpp, positiveValueGuard = 10000,

  tableFits = tableFits,

  listPartitions = NULL,

  minCladeSizes = 20,

  maxCladePartitionLevel = 100L, maxNumNodesPerCladePartition = 1L,

  listAllowedModelTypesIndices = "best-clade",

  argsConfigOptim1 = DefaultArgsConfigOptim(
    numCallsOptim = 400,
    numRunifInitVecParams = 100000,
    numGuessInitVecParams = 50000),
  argsConfigOptim2 = DefaultArgsConfigOptim(
    numCallsOptim = 4,
    numRunifInitVecParams = 1000,
    numGuessInitVecParams = 10000,
    numJitterRootRegimeFit = 1000,
    sdJitterRootRegimeFit = 0.05,
    numJitterAllRegimeFits = 1000,
    sdJitterAllRegimeFits = 0.05),
  argsConfigOptim3 = DefaultArgsConfigOptim(
    numCallsOptim = 4,
    numRunifInitVecParams = 1000,
    numGuessInitVecParams = 10000,
    numJitterRootRegimeFit = 1000,
    sdJitterRootRegimeFit = 0.05,
    numJitterAllRegimeFits = 1000,
    sdJitterAllRegimeFits = 0.05),

  maxNumRoundRobins = 5,
  maxNumPartitionsInRoundRobins = 3,

  doParallel = TRUE,

  prefixFiles = prefixFiles,
  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = FALSE,
  verbose = TRUE,
  debug = FALSE)

save(fitMappings, file = paste0("FinalResult_", prefixFiles, ".RData"))

if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

