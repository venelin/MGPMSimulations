# example calls:
# R --vanilla -f CollectFits_surfaceFwd_t5.R

library(PCMFit)
library(PCMBase)
library(MGPMSimulations)
library(data.table)
# tested with surface v 0.4-1
library(surface)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)


dataNameFwd <- "fits_surface_fwd_t5"

resultDirPrefix <- "Results_t5_N80_surfaceFwdBwd"
# e.g. MGPM_A_F_best_clade_2_id_
resultDirPrefix2 <- "surfaceFwdBwd_id_"
resultFilePrefixFwd <- "FinalResult_fwd_surfaceFwdBwd_id_"

Args_MixedGaussian_MGPMSurfaceOUType2 <- function() {
  list(
    X0 = structure(
      0.0,
      class = c("VectorParameter", "_Fixed", "_Global"),
      description = "root value fixed to the parameter theta of regime 1."),
    H = structure(
      0.0,
      class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      description = "adaptation rate matrix"),
    Sigma_x = structure(
      0.0,
      class = c("MatrixParameter", "_Diagonal", "_WithNonNegativeDiagonal", "_Global"),
      description = "unit-time variance parameter of the OU-process"),
    Sigmae_x = structure(
      0.0, class = c("MatrixParameter", "_Omitted"),
      description = "upper triangular Choleski factor of the non-phylogenetic variance-covariance")
  )
}

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

getTreeWithRegimesFromHansenFit <- function(
  tree, hansenfit, cols = NULL, convcol = TRUE, labelshifts = FALSE,
  ...) {

  # copied and slightly modified from package surface, function surfaceTreePlot:

  fit <- hansenfit$fit[[1]]
  otree <- as(fit, "data.frame")
  otree <- data.frame(otree, shifts = rep(NA, length(otree$nodes)))
  otree$shifts[match(names(hansenfit$savedshifts), otree$nodes)] <- 1:length(hansenfit$savedshifts)
  ntip <- (dim(otree)[1] + 1)/2
  nnode <- ntip - 1
  otree2 <- otree[match(c(tree$tip.label, tree$node.label),
                        otree$labels), ]
  otree2 <- otree2[tree$edge[, 2], ]
  if (length(cols) == 1)
    cols <- rep(cols, length(unique(hansenfit$savedshifts)))
  if (is.null(cols)) {
    xx <- summary(factor(hansenfit$savedshifts))
    if (convcol) {
      cols <- character(length(xx))
      cols[xx > 1] <- rainbow(sum(xx > 1))
      if (any(xx == 1))
        cols[xx == 1] <- c("black", grey(seq(0.7, 0.3,
                                             length.out = sum(xx == 1) - 1)))
    }
    else {
      cols <- c("black", rainbow(length(xx) - 1))
    }
  }
  edgecols <- cols[as.numeric(factor(otree2[, 5]))]
  tree$edge.regime <- as.character(factor(otree2[, 5]))
  tree
  # plot(tree, edge.color = edgecols, ...)
  # if (labelshifts) {
  #   nodelabels(node = tree$edge[, 2][which(!is.na(otree2$shifts))],
  #              bg = "white", text = otree2$shifts[!is.na(otree2$shifts)],
  #              cex = 0.6, frame = "circle")
  # }
}

id <- testData_t5_fittedIds[[1]]
data <- rbindlist(
  lapply(
    testData_t5_fittedIds,
    function(id) {
      if(id %% 10L == 0) {
        cat("id:", id, "\n")
      }

      NSize <- paste0("N", testData_t5[id, nobs])

      fileFwd <- paste0(resultDirPrefix, "/", resultDirPrefix2, id , "/", resultFilePrefixFwd, id, ".RData")

      if(file.exists(fileFwd)) {
        # loads an object called fwd
        load(fileFwd)

        sfwd <- surfaceSummary(fwd)

        tree <- testData_t5[IdGlob == id]$tree[[1]]
        X <- testData_t5[IdGlob == id]$X[[1]]


        modelMGPMSurfaceFwd <- do.call(
          MixedGaussian,
          c(list(
            k = 2,
            modelTypes = MGPMSurfaceOUType(),
            mapping = structure( rep(1L, nrow(sfwd$theta)),
                                 names = rownames(sfwd$theta)) ),
            Args_MixedGaussian_MGPMSurfaceOUType2()) )


        modelMGPMSurfaceFwd$X0[] <- sfwd$theta[1, ]
        diag(modelMGPMSurfaceFwd$H) <- sfwd$alpha
        diag(modelMGPMSurfaceFwd$Sigma_x) <- sqrt(sfwd$sigma_squared)

        for(regimeName in rownames(sfwd$theta)) {
          modelMGPMSurfaceFwd[[regimeName]]$Theta[, 1] <- sfwd$theta[regimeName, ]
        }

        vecParams <- PCMParamGetShortVector(modelMGPMSurfaceFwd)

        PCMParamLoadOrStore(
          modelMGPMSurfaceFwd, vecParams,
          offset = 0, k = 2, R = PCMNumRegimes(modelMGPMSurfaceFwd), load = TRUE)

        treeWithRegimesFwd <- getTreeWithRegimesFromHansenFit(tree, fwd[[length(fwd)]])

        ll <- PCMLik(X, treeWithRegimesFwd, modelMGPMSurfaceFwd)

        attr(modelMGPMSurfaceFwd, "tree") <- treeWithRegimesFwd
        attr(modelMGPMSurfaceFwd, "X") <- X[, seq_len(PCMTreeNumTips(treeWithRegimesFwd))]
        attr(modelMGPMSurfaceFwd, "SE") <- attr(modelMGPMSurfaceFwd, "X") * 0.0
        sc <- AICc_SURFACE(modelMGPMSurfaceFwd)

        testData_t5[id,
                    list(
                      IdGlob, IdTree, IdClusteringForTree,
                      IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                      clusterNodes = list(PCMTreeGetPartRegimes(PCMTree(treeWithRegimesFwd))),
                      mapping = list(rep(1L, PCMNumRegimes(modelMGPMSurfaceFwd))),
                      model = list(modelMGPMSurfaceFwd),
                      logLik = ll,
                      score = sc,
                      logLikSurfacePackage = colSums(sbwd$lnls)[which.min(sbwd$aics)],
                      scoreSurfacePackage = min(sfwd$aics))]
      } else {
        testData_t5[id,
                    list(
                      IdGlob, IdTree, IdClusteringForTree,
                      IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                      clusterNodes = list(NULL), mapping = list(NULL),
                      model = list(NULL),
                      logLik = NA_real_, score = NA_real_,
                      logLikSurfacePackage = NA_real_,
                      scoreSurfacePackage = NA_real_)]
      }

    }))

eval(parse(text = paste0(dataNameFwd, " <- data")))
eval(parse(text = paste0("usethis::use_data(", dataNameFwd, ", overwrite = TRUE)")))
