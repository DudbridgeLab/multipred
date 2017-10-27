#' Analytic joint measures of predictive accuracy
#'
#' Analytic calculation of joint sensitivity, specificity, positive and negative predictive value, concordance and net benefit, under a multivariate liability threshold model.
#'
#' Details here

#' @param VL Variance-covariance matrix of liability.  Must have 1 on diagonal.
#' @param VX Variance-covariance matrix of predictors.  Diagonal entries are the liability variances explained for each trait ("heritabilities")
#' @param VLX Cross-covariance matrix between liabilities and predictors.  Entry on row i, column j, is covariance between liability i and predictor j.
#' @param thresh Vector of risk thresholds for predicting an event.
#' @param prev Vector of prevalences, ie population risks, for each trait
#' @param target Target trait vector for calculating joint measures of accuracy.
#' @param target2 Second target vector for calculating joint concordance.
#' @param targetProb Probability of all the events occurring in \code{target}, used for calculating net benefit.
#'
#' @export
analyticJoint = function(VL,VX,VLX,thresh,prev,target,target2=NULL,targetProb) {

  # coerce VL, VX and VLX to matrices
  VL = as.matrix(VL)
  VX = as.matrix(VX)
  VLX = as.matrix(VLX)

  # coerce prev, thresh and target to vectors
  prev = as.vector(prev)
  thresh = as.vector(thresh)
  target = as.vector(target)

  check=checkMatrixDimensions(VL,VX)
  if (!is.null(check)) {
    print(paste("VL and VX",check))
    return(NULL)
  }

  check=checkMatrixDimensions(VL,VLX)
  if (!is.null(check)) {
    print(paste("VL and VLX",check))
    return(NULL)
  }

  check=checkVectorMatrixDimensions(thresh,VL)
  if (!is.null(check)) {
    print(paste("thresh",check))
    return(NULL)
  }

  check=checkVectorMatrixDimensions(prev,VL)
  if (!is.null(check)) {
    print(paste("prev",check))
    return(NULL)
  }

    check=checkVectorMatrixDimensions(target,VL)
  if (!is.null(check)) {
    print(paste("target",check))
    return(NULL)
  }

  # liability thresholds
  liabThresh = qnorm(prev,lower=F)

  # scores corresponding to the risk thresholds
  scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VLX))

  # truncation points for multivariate normal liability
  lowerTrunc = rep(-Inf,length(target))
  lowerTrunc[target==1] = liabThresh[target==1]
  upperTrunc = rep(Inf,length(target))
  upperTrunc[target==0] = liabThresh[target==0]

  # mean liability in subjects with target traits
  liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=lowerTrunc, upper=upperTrunc)

  # mean score in selected subjects
  invVL = solve(VL)
  meanScoreTrunc = as.vector(VLX %*% invVL %*% liabTrunc$tmean)
  # covariance matrix of scores in selected subjects
  varScoreTrunc = VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
  # make varScoreTrunc symmetrical
  varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2

  # sensitivity
  lowerSens = scoreThresh
  lowerSens[target==0] = -Inf
  upperSens = rep(Inf,length(target))
  sens = pmvnorm(lower=lowerSens, upper=upperSens, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(sens) = NULL

  # specificity
  lowerSpec = rep(-Inf,length(target))
  upperSpec = scoreThresh
  upperSpec[target==1] = Inf
  spec = pmvnorm(lower=lowerSpec, upper=upperSpec, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(spec) = NULL

  # truncation points for multivariate normal scores
  lowerTrunc = rep(-Inf,length(target))
  lowerTrunc[target==1] = scoreThresh[target==1]
  upperTrunc = rep(Inf,length(target))
  upperTrunc[target==0] = scoreThresh[target==0]

  # mean score in subjects with target predictions
  scoreTrunc = tmvtnorm::mtmvnorm(sigma=VX, lower=lowerTrunc, upper=upperTrunc)

  # mean liability in selected subjects
  invVX = solve(VX)
  meanLiabTrunc = as.vector(t(VLX) %*% invVX %*% scoreTrunc$tmean)
  # covariance matrix of liability in selected subjects
  varLiabTrunc = VL - t(VLX) %*% (invVX - invVX %*% scoreTrunc$tvar %*% invVX) %*% VLX
  # make varLiabTrunc symmetrical
  varLiabTrunc = (varLiabTrunc + t(varLiabTrunc)) / 2

  # positive predictive value
  lowerPPV = liabThresh
  lowerPPV[target==0] = -Inf
  upperPPV = rep(Inf,length(target))
  PPV = pmvnorm(lower=lowerPPV, upper=upperPPV, mean=meanLiabTrunc, sigma=varLiabTrunc)
  attributes(PPV) = NULL

  # negative predictive value
  lowerNPV = rep(-Inf,length(target))
  upperNPV = liabThresh
  upperNPV[target==1] = Inf
  NPV = pmvnorm(lower=lowerNPV, upper=upperNPV, mean=meanLiabTrunc, sigma=varLiabTrunc)
  attributes(NPV) = NULL

  # concordance
  if (is.null(target2)) {
    C = 0
  }
  else {
    # distribution of scores for subjects with target2
    # truncation points for multivariate normal liability
    lowerTrunc = rep(-Inf,length(target2))
    lowerTrunc[target2==1] = liabThresh[target2==1]
    upperTrunc = rep(Inf,length(target2))
    upperTrunc[target2==0] = liabThresh[target2==0]

    # mean liability in subjects with target2 traits
    liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=lowerTrunc, upper=upperTrunc)

    # distribution of liability differeces between subjects with target and target2
    # mean score in selected subjects
    invVL = solve(VL)
    meanScoreTrunc = meanScoreTrunc - as.vector(VLX %*% invVL %*% liabTrunc$tmean)
    # covariance matrix of scores in selected subjects
    varScoreTrunc = varScoreTrunc + VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
    # make varScoreTrunc symmetrical
    varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2

    # concordance
    lowerTrunc = rep(-Inf,length(thresh))
    lowerTrunc[target>target2] = 0
    upperTrunc = rep(Inf,length(thresh))
    upperTrunc[target<target2] = 0
    C = pmvnorm(lower=lowerTrunc, upper=upperTrunc, mean=meanScoreTrunc, sigma=varScoreTrunc)
    attributes(C) = NULL

  }

  # net benefit calculations
  # lower bound for cost/benefit ratio
  costBenefit = prod(target*thresh) / (1-prod(target*thresh))

  # family-wise specificity for complementary target
  # truncation points for multivariate normal liability
  lowerTrunc = rep(-Inf,length(target))
  lowerTrunc[target==0] = liabThresh[target==0]
  upperTrunc = rep(Inf,length(target))
  upperTrunc[target==1] = liabThresh[target==1]

  # mean liability in subjects with target traits
  liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=lowerTrunc, upper=upperTrunc)

  # mean score in selected subjects
  invVL = solve(VL)
  meanScoreTrunc = as.vector(VLX %*% invVL %*% liabTrunc$tmean)
  # covariance matrix of scores in selected subjects
  varScoreTrunc = VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
  # make varScoreTrunc symmetrical
  varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2

  # specificity
  lowerSpec = scoreThresh
  lowerSpec[target==0] = -Inf
  upperSpec = rep(Inf,length(target))
  familyWiseSpec = 1-pmvnorm(lower=lowerSpec, upper=upperSpec, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(familyWiseSpec) = NULL

  # finally, the net benefit
  NB = sens - (1-targetProb)/targetProb*costBenefit * (1-familyWiseSpec)

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              NB=NB))
}
