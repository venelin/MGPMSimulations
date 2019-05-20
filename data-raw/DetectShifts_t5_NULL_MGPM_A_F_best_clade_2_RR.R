# This script can be run locally or on a cluster using a command like:
#
# bsub -M 200000 -n 200 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args 96
#
# bsub -M 30000 -n 20 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_NULL_MGPM_A_F_best_clade_2_RR.R --args $id
# In the above command the argument 96 after --args denotes the row-number in
# testData_t5_NULL data.table from the package MGPMSimulations.
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
} else {
  id <- 1
}

prefixFiles = paste0("t5_NULL_MGPM_A_F_best_clade_2_RR_id_", id)

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

# try using a previously stored tableFits from a previous run that was interupted
tableFitsCurrent <- NULL
fileCurrentResults <- paste0("CurrentResults_", prefixFiles, ".RData")
if(file.exists(fileCurrentResults)) {
  cat("Loading previously stored tableFits from file", fileCurrentResults, "...\n")

  # this should load a list object named listResults containg a named entry "tableFits":
  load(fileCurrentResults)

  tableFitsCurrent <- listResults$tableFits

  setkey(tableFitsCurrent, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)
}

tempFiles <- list.files(pattern = paste0("^", prefixFiles, ".*.RData"))
if(length(tempFiles) > 0) {
  cat("Loading previously stored tableFits from temporary files (", toString(tempFiles), ")...\n")
  tableFitsTempFiles <- rbindlist(
    lapply(tempFiles, function(file) {
      load(file)
      fits
    }))
  tableFitsCurrent <- PCMFit:::UpdateTableFits(tableFitsCurrent, tableFitsTempFiles)
  setkey(tableFitsCurrent, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)
}

tableFits <- tableFitsCurrent

print(tableFits)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-4)
options(PCMBase.Threshold.SV = 1e-5)
options(PCMBase.Threshold.Lambda_ij = 1e-4)

print(PCMOptions())

tree <- testData_t5_NULL$tree[[id]]
values <- testData_t5_NULL$X[[id]][, seq_len(PCMTreeNumTips(tree))]

options(MGPMSimulations.LowerLimitTheta = apply(values, 1, min))
options(MGPMSimulations.UpperLimitTheta = apply(values, 1, max))

numTraits <- testData_t5_NULL$numTraits[[id]]
logLik <- testData_t5_NULL$logLik[[id]]

fitMappings <- PCMFitMixed(

  values, tree,

  metaIFun = PCMInfoCpp, positiveValueGuard = 100,

  tableFits = tableFits,

  listPartitions = NULL,

  minCladeSizes = 20,

  maxCladePartitionLevel = 100L, maxNumNodesPerCladePartition = 1L,

  listAllowedModelTypesIndices = "best-clade-2",

  argsConfigOptim1 = DefaultArgsConfigOptim(
    numCallsOptim = 100,
    numRunifInitVecParams = 100000 * numTraits / 2,
    numGuessInitVecParams = 100000 * numTraits / 2),
  argsConfigOptim2 = DefaultArgsConfigOptim(
    numCallsOptim = 4 * numTraits / 2,
    numRunifInitVecParams = 1000 * numTraits / 2,
    numGuessInitVecParams = 10000 * numTraits / 2,
    numJitterRootRegimeFit = 1000 * numTraits / 2,
    sdJitterRootRegimeFit = 0.05,
    numJitterAllRegimeFits = 1000 * numTraits / 2,
    sdJitterAllRegimeFits = 0.05),
  argsConfigOptim3 = DefaultArgsConfigOptim(
    numCallsOptim = 4,
    numRunifInitVecParams = 1000 * numTraits / 2,
    numGuessInitVecParams = 10000 * numTraits / 2,
    numJitterRootRegimeFit = 1000 * numTraits / 2,
    sdJitterRootRegimeFit = 0.05,
    numJitterAllRegimeFits = 1000 * numTraits / 2,
    sdJitterAllRegimeFits = 0.05),

  maxNumRoundRobins = 5,
  maxNumPartitionsInRoundRobins = 3,

  listPCMOptions = c(
    PCMOptions(),
    options()[
      c('MGPMSimulations.LowerLimitTheta', 'MGPMSimulations.UpperLimitTheta')]),

  doParallel = TRUE,

  prefixFiles = prefixFiles,
  saveTempWorkerResults = TRUE,
  printFitVectorsToConsole = TRUE,
  verbose = TRUE,
  debug = FALSE)

save(fitMappings, file = paste0("FinalResult_", prefixFiles, ".RData"))

if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

