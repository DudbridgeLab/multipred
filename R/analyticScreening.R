#' Analytic screening measures of predictive accuracy
#'
#' Analytic calculation of screening sensitivity, specificity, positive and negative predictive value, concordance and net benefit, under a multivariate liability threshold model.
#'
#' Details here

#' @param VL Variance-covariance matrix of liability.  Must have 1 on diagonal.
#' @param VX Variance-covariance matrix of predictors.  Diagonal entries are the liability variances explained for each trait ("heritabilities")
#' @param VLX Cross-covariance matrix between liabilities and predictors.  Entry on row i, column j, is covariance between liability i and predictor j.
#' @param thresh Vector of risk thresholds for predicting an event.
#' @param prev Vector of prevalences, ie population risks, for each trait
#'
#' @export
analyticScreening = function(VL,VX,VLX,thresh,prev) {

  # coerce VL, VX and VLX to matrices
  VL = as.matrix(VL)
  VX = as.matrix(VX)
  VLX = as.matrix(VLX)

  # coerce prev, thresh and target to vectors
  prev = as.vector(prev)
  thresh = as.vector(thresh)

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

  check=checkVectorMatrixDimensions(prev,VL)
  if (!is.null(check)) {
    print(paste("prev",check))
    return(NULL)
  }

  check=checkVectorMatrixDimensions(thresh,VL)
  if (!is.null(check)) {
    print(paste("thresh",check))
    return(NULL)
  }

  # liability thresholds
  liabThresh = qnorm(prev,lower=F)

  # scores corresponding to the risk thresholds
  scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VLX))

  # truncation points for multivariate normal liability
  lowerTrunc = rep(-Inf,length(thresh))
  upperTrunc = liabThresh

  # prevalence of subjects with no traits
  prev0 = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VL)
  attributes(prev0) = NULL

  # mean liability in subjects with at least one trait
  liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=lowerTrunc, upper=upperTrunc)
  liabTruncMean = -prev0/(1-prev0) * liabTrunc$tmean

  # variance of liability in subjects with at least one trait
  liabTruncVar = (VL - prev0 * (liabTrunc$tvar + liabTrunc$tmean %*% t(liabTrunc$tmean))) / (1-prev0) -
    prev0^2/(1-prev0)^2 * liabTrunc$tmean %*% t(liabTrunc$tmean)

  # mean score in selected subjects
  invVL = solve(VL)
  meanScoreTrunc = as.vector(VLX %*% invVL %*% liabTruncMean)
  # covariance matrix of scores in selected subjects
  varScoreTrunc = VX - VLX %*% (invVL - invVL %*% liabTruncVar %*% invVL) %*% t(VLX)
  # make varScoreTrunc symmetrical
  varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2

  # sensitivity
  lowerSens = rep(-Inf,length(thresh))
  upperSens = scoreThresh
  sens = 1 - pmvnorm(lower=lowerSens, upper=upperSens, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(sens) = NULL

  # specificity
  meanScoreTrunc = as.vector(VLX %*% invVL %*% liabTrunc$tmean)
  varScoreTrunc = VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
  varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2
  lowerSpec = rep(-Inf,length(thresh))
  upperSpec = scoreThresh
  spec = pmvnorm(lower=lowerSpec, upper=upperSpec, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(spec) = NULL

  ### NOW FOR PPV AND NPV

  # prevalence of subjects with no predicted traits
  lowerTrunc=rep(-Inf,length(thresh))
  upperTrunc=scoreThresh
  prev0 = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VX)
  attributes(prev0) = NULL

  # mean score in subjects with no predicted traits
  scoreTrunc = tmvtnorm::mtmvnorm(sigma=VX, lower=lowerTrunc, upper=upperTrunc)
  scoreTruncMean = -prev0/(1-prev0) * scoreTrunc$tmean

  # variance of score in subjects with no predicted traits
  scoreTruncVar = (VX - prev0 * (scoreTrunc$tvar + scoreTrunc$tmean %*% t(scoreTrunc$tmean))) / (1-prev0) -
    prev0^2/(1-prev0)^2 * scoreTrunc$tmean %*% t(scoreTrunc$tmean)

  # mean liability in selected subjects
  invVX = solve(VX)
  meanLiabTrunc = as.vector(t(VLX) %*% invVX %*% scoreTruncMean)
  # covariance matrix of liabilities in selected subjects
  varLiabTrunc = VL - t(VLX) %*% (invVX - invVX %*% scoreTruncVar %*% invVX) %*% VLX
  # make varLiabTrunc symmetrical
  varLiabTrunc = (varLiabTrunc + t(varLiabTrunc)) / 2

  # PPV
  lowerPPV = rep(-Inf,length(thresh))
  upperPPV = liabThresh
  PPV = 1 - pmvnorm(lower=lowerPPV, upper=upperPPV, mean=meanLiabTrunc, sigma=varLiabTrunc)
  attributes(PPV) = NULL

  # NPV
  meanLiabTrunc = as.vector(t(VLX) %*% invVX %*% scoreTrunc$tmean)
  varLiabTrunc = VL - t(VLX) %*% (invVX - invVX %*% scoreTrunc$tvar %*% invVX) %*% VLX
  varLiabTrunc = (varLiabTrunc + t(varLiabTrunc)) / 2
  lowerNPV = rep(-Inf,length(thresh))
  upperNPV = liabThresh
  NPV = pmvnorm(lower=lowerNPV, upper=upperNPV, mean=meanLiabTrunc, sigma=varLiabTrunc)
  attributes(NPV) = NULL

  # concordance
  # truncation points for multivariate normal liability
  lowerTrunc = rep(-Inf,length(thresh))
  upperTrunc = liabThresh

  # prevalence of subjects with no traits
  prev0 = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VL)
  attributes(prev0) = NULL

  # mean liability in subjects with no traits
  liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=lowerTrunc, upper=upperTrunc)
  # in subjects with at least one trait
  liabTruncMean = -prev0/(1-prev0) * liabTrunc$tmean

  # variance of liability in subjects with at least one trait
  liabTruncVar = (VL - prev0 * (liabTrunc$tvar + liabTrunc$tmean %*% t(liabTrunc$tmean))) / (1-prev0) -
    prev0^2/(1-prev0)^2 * liabTrunc$tmean %*% t(liabTrunc$tmean)

  # mean score in subjects with at least one trait
  invVL = solve(VL)
  meanScoreTrunc = as.vector(VLX %*% invVL %*% liabTruncMean)
  # covariance matrix of scores in selected subjects
  varScoreTrunc = VX - VLX %*% (invVL - invVL %*% liabTruncVar %*% invVL) %*% t(VLX)

  # now subtract mean score in subjects with no traits
  meanScoreTrunc = meanScoreTrunc - as.vector(VLX %*% invVL %*% liabTrunc$tmean)
  varScoreTrunc = varScoreTrunc + VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
  # make varScoreTrunc symmetrical
  varScoreTrunc = (varScoreTrunc + t(varScoreTrunc)) / 2

  lowerTrunc=rep(-Inf,length(thresh))
  upperTrunc=rep(0,length(thresh))
  C = 1-pmvnorm(lower=lowerTrunc, upper=upperTrunc, mean=meanScoreTrunc, sigma=varScoreTrunc)
  attributes(C) = NULL

  # net benefit
  # truncation points for multivariate normal liability
  lowerTrunc = rep(-Inf,length(thresh))
  upperTrunc = liabThresh

  # prevalence of subjects with at least one trait
  targetProb = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VL)
  attributes(targetProb) = NULL

  NB = sens - (1-spec) * min(thresh/(1-thresh)) * (1-targetProb)/targetProb

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              NB=NB))
}
