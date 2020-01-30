#' Analytic outcome-wise measures of predictive accuracy
#'
#' Analytic calculation of outcome-wise sensitivity, specificity, positive and negative predictive value, concordance and relative utility, under a multivariate liability threshold model.
#'
#' Outcome-wise measures consider the prediction of individual outcomes summed over individuals.
#' When \code{weight} is a vector of 1's (default), outcome-wise measures correspond to classical univariate measures with the \code{x} matrix vectorised into a column vector.
#' More generally, \code{weight} allows different outcomes to contribute more or less to the calculations.
#'
#' Outcome-wise sensitivity, specificity and concordance are weighted sums of the univariate measures,
#' where the weights depend on \code{prev}.
#'
#' @template analyticParams
#' @param weight Vector of weights.
#'
#' @examples
#' attach(PRSdata)
#' analyticOutcomeWise(VL,VX,VX,thresh=prevalence,prev=prevalence)
#'
#' # $sens
#' # [1] 0.6243863
#'
#' # $spec
#' # [1] 0.6132883
#'
#' # $PPV
#' # [1] 0.04641913
#'
#' # $NPV
#' # [1] 0.9818697
#'
#' # $C
#' # [1] 0.6533142
#'
#' # $RU
#' # [1]  0.2376747
#'
#' @export
analyticOutcomeWise = function(VL,VX,VLX=NULL,thresh=NULL,weight=NULL,prev) {

  # coerce VL, VX and VLX to matrices
  if (is.null(VLX)) VLX = VX
  VL = as.matrix(VL)
  VX = as.matrix(VX)
  VLX = as.matrix(VLX)

  # coerce prev, thresh and target to vectors
  prev = as.vector(prev)
  thresh = as.vector(thresh)
  weight = as.vector(weight)

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

  if (!is.null(thresh)) {
    check=checkVectorMatrixDimensions(thresh,VL)
    if (!is.null(check)) {
      print(paste("thresh",check))
      return(NULL)
    }
  }

  if (!is.null(weight)) {
    check=checkVectorMatrixDimensions(weight,VL)
    if (!is.null(check)) {
      print(paste("weight",check))
      return(NULL)
    }
  }

  check=checkVectorMatrixDimensions(prev,VL)
  if (!is.null(check)) {
    print(paste("prev",check))
    return(NULL)
  }

  ntrait = ncol(VL)
  if (is.null(weight)) weight=rep(1,ntrait)

  sens = NULL
  spec = NULL
  PPV = NULL
  NPV = NULL
  C = NULL
  RU = NULL

  # liability thresholds
  liabThresh = qnorm(prev,lower=F)

  if (!is.null(thresh)) {

    # scores corresponding to the risk thresholds
    scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VX))

    probPredicted = pnorm(scoreThresh, sd=sqrt(diag(VX)), lower=F)
    probOutcomePredicted = NULL
    for(i in 1:ntrait)
      probOutcomePredicted[i] = pmvnorm(lower=c(liabThresh[i],scoreThresh[i]), upper=rep(Inf,2), sigma=rbind(c(VL[i,i],VLX[i,i]),c(VLX[i,i],VX[i,i])))

    # sensitivity
    sens = probOutcomePredicted / prev
    sens = t(prev * weight) %*% sens / sum(prev * weight)
    sens = as.vector(sens)

    # specificity
    spec = 1 - (probPredicted-probOutcomePredicted) / (1-prev)
    spec = t((1-prev) * weight) %*% spec / sum((1-prev) * weight)
    spec = as.vector(spec)

    # positive predictive value
    PPV = probOutcomePredicted / probPredicted
    PPV = t(probPredicted * weight) %*% PPV / sum(probPredicted * weight)
    PPV = as.vector(PPV)

    # negative predictive value
    NPV = 1 - (prev-probOutcomePredicted) / (1-probPredicted)
    NPV = t ((1-probPredicted) * weight) %*% NPV / sum((1-probPredicted) * weight)
    NPV = as.vector(NPV)

    # relative utility
    RU = sens - (thresh %*% weight) / ((1-thresh) %*% weight) * ((1-prev) %*% weight) / (prev %*% weight) * (1-spec)
    RU = as.vector(RU)
  }

  # concordance
  threshDensity = dnorm(liabThresh) / prev
  caseMean = diag(VLX) * threshDensity
  caseVar = diag(VLX) * (1 - diag(VLX) * threshDensity * (threshDensity - liabThresh))
  threshDensity = -dnorm(liabThresh) / (1-prev)
  controlMean = diag(VLX) * threshDensity
  controlVar = diag(VLX) * (1 - diag(VLX) * threshDensity * (threshDensity - liabThresh))
  C = pnorm((caseMean - controlMean) / sqrt(caseVar + controlVar))
  C = (prev * (1-prev) * weight) %*% C / sum(prev * (1-prev) * weight)
  C = as.vector(C)

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              RU=RU))
}
