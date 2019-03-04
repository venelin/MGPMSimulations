# This script can be run locally or on a cluster using a command like:
#
# bsub -n 1 -W 3:59 sh R --vanilla --slave -f ../../DetectShifts_t5_surfaceFwdBwd.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package MGPMSimulations

library(PCMBase)
library(ape)
library(MGPMSimulations)
library(surface)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
} else {
  id <- 1
}

prefixFiles = paste0("surfaceFwdBwd_id_", id)

tree <- testData_t5$tree[[id]]
dat <- as.data.frame(t(testData_t5$X[[id]][, seq_len(PCMTreeNumTips(tree))]))

olist <- convertTreeData(tree,dat)
otree <- olist[[1]]; odata <- olist[[2]]

fwd <- surfaceForward(otree, odata, aic_threshold = 0, exclude = 0,
                      verbose = TRUE, plotaic = FALSE)

save(fwd, file = paste0("FinalResult_fwd_", prefixFiles, ".RData"))

k <- length(fwd)

bwd <- surfaceBackward(otree, odata, fwd[[k]])
save(bwd, file = paste0("FinalResult_bwd_", prefixFiles, ".RData"))

