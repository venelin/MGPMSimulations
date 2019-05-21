#' @export
PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- getOption(
        'MGPMSimulations.LowerSigma_xUppTri', -.0)
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[, , r][matUppTri] <- getOption(
          'MGPMSimulations.LowerSigma_xUppTri', -.0)
      }
    }
  }
  o
}

#' @export
PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[matDiag] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[, , r][matDiag] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[, , r][matUppTri] <- 1.0
      }
    }
  }
  o
}

#' @export
PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Theta)) {
    o$Theta[] <- rep(getOption(
      'MGPMSimulations.LowerLimitTheta', c(0.0, -1.2)), length.out = k)
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- rep(getOption(
        'MGPMSimulations.LowerLimitTheta', c(0.0, -1.2)), length.out = k)
    }
  }
  if(is.Global(o$Sigma_x)) {
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- getOption(
        'MGPMSimulations.LowerSigma_xUppTri', -.0)
    }
  } else {
    if(!is.Diagonal(o$Sigma_x)) {
      for(r in seq_len(R)) {
        o$Sigma_x[, , r][matUppTri] <- getOption(
          'MGPMSimulations.LowerSigma_xUppTri', -.0)
      }
    }
  }
  o
}

#' @export
PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  mat <- matrix(FALSE, k, k)
  matUppTri <- upper.tri(mat)
  matDiag <- mat
  diag(matDiag) <- TRUE

  if(is.Global(o$Theta)) {
    o$Theta[] <- rep(getOption(
      'MGPMSimulations.UpperLimitTheta', c(7.8, 4.2)), length.out = k)
  } else {
    for(r in seq_len(R)) {
      o$Theta[, r] <- rep(getOption(
        'MGPMSimulations.UpperLimitTheta', c(7.8, 4.2)), length.out = k)
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[matDiag] <- 1.0
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[matUppTri] <- 1.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[, , r][matDiag] <- 1.0
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[, , r][matUppTri] <- 1.0
      }
    }
  }
  o
}
