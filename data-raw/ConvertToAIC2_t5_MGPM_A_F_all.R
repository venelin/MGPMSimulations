# This script modifies the fits for trees of 80 tips by subtracting the number of
# regimes in each fit from the AIC score. This is in response to a reviewer's
# comment about adding 2 instead of 1 to the AIC score for each new regime in
# the tree. We can do this for the trees of 80 tips only since these are the
# trees for which we could conduct a search through all MGPMs (minCladeSize = 20).
# Hence, it is not needed to re-run the mode fits on a computing cluster. Instead,
# we only adjust the score based on the previously calculated AIC score.
#
#
# sh R --vanilla --slave -f ../../DetectShifts_t5_MGPM_A_F_all.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package MGPMSimulations The blank XXX stays for the type of model to be inferred, e.g. ScalarOU or MixedGaussian
#
library(ape)
library(PCMBase)
library(PCMFit)
library(data.table)
library(MGPMSimulations)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
} else {
  id <- 1
}

prefixFiles = paste0("MGPM_A_F_all_id_", id)

prefixFilesAIC2 <- paste0("MGPM_A_F_all_AIC2_id_", id)

fileFinalResultsAIC <- paste0("../../Results_t5_MGPM_A_F_all_N80/", prefixFiles, "/FinalResult_", prefixFiles, ".RData")
if(file.exists(fileFinalResultsAIC)) {
  cat("Loading ", fileFinalResultsAIC, "...\n")

  # this should load a list object named listResults containg a named entry "tableFits":
  load(fileFinalResultsAIC)

  fitMappings$tableFits[, score:=score - sapply(startingNodesRegimesLabels, length)]
  fitMappings$tableFits[, df:=df - sapply(startingNodesRegimesLabels, length)]
  fitMappings$arguments$scoreFun <-  function(object, ..., k=2) {
    s <- AIC(object)
    s - PCMNumRegimes(object$modelOptim)
  }

  save(fitMappings, file = paste0("FinalResult_", prefixFilesAIC2, ".RData"))
} else {
  stop(paste0("File not found ", fileFinalResultsAIC, "."))
}


