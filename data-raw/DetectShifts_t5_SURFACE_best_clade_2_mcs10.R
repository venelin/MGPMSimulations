# This script can be run locally or on a cluster using a command like:
#
# bsub -n 10 -W 23:59 -R ib sh R --vanilla --slave -f ../../DetectShifts_t5_SURFACE_best_clade_2.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package MGPMSimulations The blank XXX stays for the type of model to be inferred, e.g. SURFACE or MixedGaussian
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

prefixFiles = paste0("SURFACE_best_clade_2_mcs10_id_", id)

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

tableFits <- NULL
# try using a previously stored tableFits from a previous run that was interupted
fileCurrentResults <- paste0("CurrentResults_", prefixFiles, ".RData")
if(file.exists(fileCurrentResults)) {
  cat("Loading previously stored tableFits from file", fileCurrentResults, "...\n")

  # this should load a list object named listResults containg a named entry "tableFits":
  load(fileCurrentResults)

  tableFits <- listResults$tableFits

  tempFiles <- list.files(pattern = paste0("^", prefixFiles, ".*.RData"))
  if(length(tempFiles) > 0) {
    cat("Loading previously stored tableFits from temporary files (", toString(tempFiles), ")...\n")
    tableFitsTempFiles <- rbindlist(
      lapply(tempFiles, function(file) {
        load(file)
        fits
      }))
    tableFits <- rbindlist(list(tableFits, tableFitsTempFiles))
  }

  setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

  tableFits <- unique(tableFits, by = key(tableFits))
  tableFits[, duplicated:=FALSE]

  print(tableFits)
}

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-7)

print(PCMOptions())

AICc_SURFACE <- function(object, ..., k = 2) {
  llValue <- logLik(object)

  # number of parameters. The df attribute is calculated by the PCMFit package
  # using a call:
  # PCMParamCount(
  #    object$model,
  #    countRegimeChanges = TRUE,
  #    countModelTypes = TRUE)
  # The countModelTypes = TRUE parameter does not matter since modelTypes is
  # of length 1, hence there is no additional degree of freedom.
  # According to Ingam and Mahler 2013, "counting the placement of the ancestral
  # regime as a shift", hence 1 is added.
  p <- attr(llValue, "df") + 1

  # total number of observations: number of traits * number of tips
  N <- 2 * attr(llValue, "nobs")

  # final AICc value:
  -2*llValue + 2*p + (2*p*(p+1)/(N-p-1))
}

tree <- testData_t5$tree[[id]]
values <- testData_t5$X[[id]][, 1:PCMTreeNumTips(tree)]

fitMappings <- PCMFitMixed(

  values, tree, modelTypes = MGPMSurfaceOUType(),

  argsMixedGaussian = Args_MixedGaussian_MGPMSurfaceOUType(),

  metaIFun = PCMInfoCpp, positiveValueGuard = 10000,

  tableFits = tableFits,

  listPartitions = NULL,

  minCladeSizes = 10,

  maxCladePartitionLevel = 100L, maxNumNodesPerCladePartition = 1L,

  listAllowedModelTypesIndices = "best-clade-2",

  scoreFun = AICc_SURFACE,

  argsConfigOptim1 = DefaultArgsConfigOptim(
    numCallsOptim = 400,
    numRunifInitVecParams = 100000,
    numGuessInitVecParams = 50000),
  argsConfigOptim2 = DefaultArgsConfigOptim(
    numCallsOptim = 20,
    numRunifInitVecParams = 1000,
    numGuessInitVecParams = 10000,
    numJitterRootRegimeFit = 1000,
    sdJitterRootRegimeFit = 0.05,
    numJitterAllRegimeFits = 1000,
    sdJitterAllRegimeFits = 0.05),

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

