#' Analytic event-wise measures of predictive accuracy
#'
#' Analytic calculation of event-wise sensitivity, specificity, positive and negative predictive value, concordance and net benefit, under a multivariate liability threshold model.
#'
#' Details here

#' @param VL Variance-covariance matrix of liability.  Must have 1 on diagonal.
#' @param VX Variance-covariance matrix of predictors.  Diagonal entries are the liability variances explained for each trait ("heritabilities")
#' @param VLX Cross-covariance matrix between liabilities and predictors.  Entry on row i, column j, is covariance between liability i and predictor j.
#' @param prev Vector of prevalences, ie population risks, for each trait
#' @param thresh Vector of risk thresholds for predicting an event.
#'
#' @export
analyticEventWise = function(VL,VX,VLX,thresh,prev) {

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

  # liability thresholds
  liabThresh = qnorm(prev,lower=F)

  # scores corresponding to the risk thresholds
  scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VLX))

  # sensitivity
  sens = NULL
  for(i in 1:length(thresh)) {
    sens[i] = pmvnorm(lower=c(liabThresh[i],scoreThresh[i]),upper=c(Inf,Inf),
                      sigma=rbind(c(VL[i,i],VLX[i,i]),c(VLX[i,i],VX[i,i]))) / prev[i]
  }
  sens = t(prev) %*% sens / sum(prev)

  # specificity
  spec = NULL
  for(i in 1:length(thresh)) {
    spec[i] = pmvnorm(lower=c(-Inf,-Inf),upper=c(liabThresh[i],scoreThresh[i]),
                      sigma=rbind(c(VL[i,i],VLX[i,i]),c(VLX[i,i],VX[i,i]))) / (1-prev[i])
  }
  spec = t(1-prev) %*% spec / sum(1-prev)

  # probabilities of positive predictions
  scoreRate = pnorm(scoreThresh,sd=sqrt(VX[i,i]),lower=F)

  # positive predictive value
  PPV = NULL
  for(i in 1:length(thresh)) {
    PPV[i] = pmvnorm(lower=c(liabThresh[i],scoreThresh[i]),upper=c(Inf,Inf),
                      sigma=rbind(c(VL[i,i],VLX[i,i]),c(VLX[i,i],VX[i,i]))) / scoreRate[i]
  }
  PPV = t(scoreRate) %*% PPV / sum(scoreRate)

  # negative predictive value
  NPV = NULL
  for(i in 1:length(thresh)) {
    NPV[i] = pmvnorm(lower=c(-Inf,-Inf),upper=c(liabThresh[i],scoreThresh[i]),
                     sigma=rbind(c(VL[i,i],VLX[i,i]),c(VLX[i,i],VX[i,i]))) / (1-scoreRate[i])
  }
  NPV = t(1-scoreRate) %*% NPV / sum(1-scoreRate)

  # concordance
  # for each trait
  C = NULL
  for(i in 1:length(thresh)) {
    threshDensity = dnorm(liabThresh[i]) / prev[i]
    caseMean = VLX[i,i] * threshDensity
    caseVar = VLX[i,i] * (1 - VLX[i,i] * threshDensity * (threshDensity - liabThresh[i]))
    threshDensity = -dnorm(liabThresh[i]) / (1-prev[i])
    controlMean = VLX[i,i] * threshDensity
    controlVar = VLX[i,i] * (1 - VLX[i,i] * threshDensity * (threshDensity - liabThresh[i]))
    C[i] = pnorm((caseMean - controlMean) / sqrt(caseVar + controlVar))
  }
  C = (prev * (1-prev)) %*% C / sum(prev * (1-prev))

  # net benefit
  NB = sens - (1-spec) * min(thresh/(1-thresh)) * sum(1-prev) / sum(prev)

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              NB=NB))
}
