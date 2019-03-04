# example calls:
# R --vanilla -f CollectFits_surfaceFwdBwd_t5.R

library(PCMFit)
library(PCMBase)
library(MGPMSimulations)
library(data.table)
# tested with surface v 0.4-1
library(surface)

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)
options(PCMBase.Threshold.EV = 1e-9)


dataNameBwd <- "fits_surface_bwd_t5"

resultDirPrefix <- "Results_t5_N80_surfaceFwdBwd"
# e.g. MGPM_A_F_best_clade_2_id_
resultDirPrefix2 <- "surfaceFwdBwd_id_"
resultFilePrefixBwd <- "FinalResult_bwd_surfaceFwdBwd_id_"

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

AICc_SURFACE <- function(object, ..., tree = NULL) {
  k = PCMNumTraits(object)
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
  if(!is.null(tree)) {
    r <- PCMTreeNumParts(tree)
    # I think that this should be with +1 to stick with the definition in
    # Ingram et. al. 2013. However, commenting out the +1 results in
    # perfect match with the surface resulting AIC.
    p <- r + PCMParamCount(object) # +1
  } else {
    p <- attr(llValue, "df") + 1
  }


  # total number of observations: number of traits * number of tips
  N <- k * attr(llValue, "nobs")

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

      fileBwd <- paste0(resultDirPrefix, "/", resultDirPrefix2, id , "/", resultFilePrefixBwd, id, ".RData")

      if(file.exists(fileBwd)) {
        # loads an object called bwd
        load(fileBwd)

        sbwd <- surfaceSummary(bwd)

        tree <- testData_t5[IdGlob == id]$tree[[1]]
        X <- testData_t5[IdGlob == id]$X[[1]]


        modelMGPMSurfaceBwd <- do.call(
          MixedGaussian,
          c(list(
            k = 2,
            modelTypes = MGPMSurfaceOUType(),
            mapping = structure( rep(1L, nrow(sbwd$theta)),
                                 names = rownames(sbwd$theta)) ),
            Args_MixedGaussian_MGPMSurfaceOUType2()) )


        modelMGPMSurfaceBwd$X0[] <- sbwd$theta[1, ]
        diag(modelMGPMSurfaceBwd$H) <- sbwd$alpha
        diag(modelMGPMSurfaceBwd$Sigma_x) <- sqrt(sbwd$sigma_squared)

        for(regimeName in rownames(sbwd$theta)) {
          modelMGPMSurfaceBwd[[regimeName]]$Theta[, 1] <- sbwd$theta[regimeName, ]
        }

        vecParams <- PCMParamGetShortVector(modelMGPMSurfaceBwd)

        PCMParamLoadOrStore(
          modelMGPMSurfaceBwd, vecParams,
          offset = 0, k = 2, R = PCMNumRegimes(modelMGPMSurfaceBwd), load = TRUE)

        treeWithRegimesBwd <- PCMTree(getTreeWithRegimesFromHansenFit(tree, bwd[[length(bwd)]]))

        ll <- PCMLik(X, treeWithRegimesBwd, modelMGPMSurfaceBwd)

        attr(modelMGPMSurfaceBwd, "tree") <- treeWithRegimesBwd
        attr(modelMGPMSurfaceBwd, "X") <- X[, seq_len(PCMTreeNumTips(treeWithRegimesBwd))]
        attr(modelMGPMSurfaceBwd, "SE") <- attr(modelMGPMSurfaceBwd, "X") * 0.0
        sc <- AICc_SURFACE(modelMGPMSurfaceBwd, tree = treeWithRegimesBwd)

        testData_t5[id,
                    list(
                      IdGlob, IdTree, IdClusteringForTree,
                      IdMappingForClustering, IdParamForMapping, IdSimulationForParam,
                      clusterNodes = list(PCMTreeGetPartRegimes(PCMTree(treeWithRegimesBwd))),
                      mapping = list(rep(1L, PCMNumRegimes(modelMGPMSurfaceBwd))),
                      model = list(modelMGPMSurfaceBwd),
                      logLik = ll,
                      score = sc,
                      logLikSurfacePackage = colSums(sbwd$lnls)[which.min(sbwd$aics)],
                      scoreSurfacePackage = min(sbwd$aics))]
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

eval(parse(text = paste0(dataNameBwd, " <- data")))
eval(parse(text = paste0("usethis::use_data(", dataNameBwd, ", overwrite = TRUE)")))
