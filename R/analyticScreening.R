#' Analytic screening measures of predictive accuracy
#'
#' Analytic calculation of screening sensitivity, specificity, positive and negative predictive value, concordance and relative utility, under a multivariate liability threshold model.
#'
#' Screening measures consider the prediction of at least one outcome to occur, without regard to whether the correct outcomes are predicted.
#' For example, screening sensitivity is the probability that, for an individual in which at least one outcome did occur, the predicted risk
#' exceeds the threshold for at least one outcome (but not necessary the ones that did occur).
#' Screening specificity is the probability that, for an individual in which no outcomes did occur, the predicted risk is lower than the
#' threshold for all outcomes.
#'
#' Screening concordance is the probability that given one individual in which at least one outcome did occur, and another in which no outcomes did occur,
#' the maximum predicted risk over all outcomes is higher in the former individual.
#' It is calculated by randomly simulating \code{nsample} such pairs of individuals from the specified model.

#' @template analyticParams
#' @template analyticnsample
#'
#' @examples
#' # results will vary due to random sampling in computing multvariate integrals
#' attach(PRSdata)
#' analyticScreening(VL,VX,VX,thresh=prevalence,prev=prevalence,nsample=1e5)
#'
#' # $sens
#' # [1] 0.9591925
#'
#' # $spec
#' # [1] 0.06055228
#'
#' # $PPV
#' # [1] 0.1604819
#'
#' # $NPV
#' # [1] 0.8879618
#'
#' # $C
#' # [1] 0.59799
#'
#' # $RU
#' # [1] -0.04479532
#'
#' @export
analyticScreening = function(VL,VX,VLX=NULL,thresh=NULL,prev,nsample=NULL) {

  # coerce VL and VX to matrices
  if (is.null(VLX)) VLX = VX
  VL = as.matrix(VL)
  VX = as.matrix(VX)
  VLX = as.matrix(VLX)

  # coerce prev and thresh to vectors
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

  if (!is.null(thresh)) {
    check=checkVectorMatrixDimensions(thresh,VL)
    if (!is.null(check)) {
      print(paste("thresh",check))
      return(NULL)
    }
  }

  check=checkVectorMatrixDimensions(prev,VL)
  if (!is.null(check)) {
    print(paste("prev",check))
    return(NULL)
  }

  ntrait = ncol(VL)

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

    # probability of no events
    probNoEvents = pmvnorm(lower=rep(-Inf,ntrait), upper=liabThresh, sigma=VL)

        # probability of no predicted events
    probNoPredictions = pmvnorm(lower=rep(-Inf,ntrait), upper=scoreThresh, sigma=VX)

    # probability of no events and no predicted events
    probNoEventsNoPredictions = pmvnorm(lower=rep(-Inf,2*ntrait), upper=c(liabThresh,scoreThresh), sigma=rbind(cbind(VL,VLX),cbind(t(VLX),VX)))

    # sensitivity
    sens = 1 - (probNoPredictions-probNoEventsNoPredictions) / (1-probNoEvents)
    attributes(sens) = NULL

    # specificity
    spec = probNoEventsNoPredictions / probNoEvents
    attributes(spec) = NULL

    # PPV
    PPV = 1 - (probNoEvents-probNoEventsNoPredictions) / (1-probNoPredictions)
    attributes(PPV) = NULL

    # NPV
    NPV = probNoEventsNoPredictions / probNoPredictions
    attributes(NPV) = NULL

    # relative utility
    # probability of no events given predictions at the threshold
    probNoEventsConditional = pmvnorm(lower=rep(-Inf,ntrait), upper=liabThresh, mean=as.vector((VLX %*% solve(VX) %*% as.matrix(thresh))), sigma=VL-VLX %*% solve(VX) %*% t(VLX))
    #RU = sens - (1-spec) * (1-prod(1-thresh))/prod(1-thresh) * probNoEvents/(1-probNoEvents)
    RU = sens - (1-spec) * (1-probNoEventsConditional)/probNoEventsConditional * probNoEvents/(1-probNoEvents)
    attributes(RU) = NULL

  }

  # concordance
  if (!is.null(nsample)) {

    # probability of no events
    probNoEvents = pmvnorm(lower=rep(-Inf,ntrait), upper=liabThresh, sigma=VL)
    attributes(probNoEvents) = NULL

    # mean liability in subjects with at least one event
    liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=rep(-Inf,ntrait), upper=liabThresh)
    liabTruncMean = -probNoEvents/(1-probNoEvents) * liabTrunc$tmean

    # variance of liability in subjects with at least one event
    liabTruncVar = (VL - probNoEvents * (liabTrunc$tvar + liabTrunc$tmean %*% t(liabTrunc$tmean))) / (1-probNoEvents) -
      probNoEvents^2/(1-probNoEvents)^2 * liabTrunc$tmean %*% t(liabTrunc$tmean)

    # mean score in subjects with at least one event
    invVL = solve(VL)
    meanScoreSens = as.vector(VLX %*% invVL %*% liabTruncMean)
    # covariance matrix of scores in subjects with at least one event
    varScoreSens = VX - VLX %*% (invVL - invVL %*% liabTruncVar %*% invVL) %*% t(VLX)
    # make varScoreSens symmetrical
    varScoreSens = (varScoreSens + t(varScoreSens)) / 2

    # mean score in subjects with no evemts
    meanScoreSpec = as.vector(VLX %*% invVL %*% liabTrunc$tmean)
    # covariance matrix of scores in subejcts with no events
    varScoreSpec = VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
    # make varScoreSpec symmetrical
    varScoreSpec = (varScoreSpec + t(varScoreSpec)) / 2

    # sample of scores in subjects with at least one event
    sensSample = rmvnorm(nsample,mean=meanScoreSens,sigma=varScoreSens)
    # convert to risks
    risk=matrix(nrow=nsample,ncol=ntrait)
    for(i in 1:ntrait) {
      risk[,i]=pnorm(liabThresh[i],mean=sensSample[,i],sd=sqrt(1-VX[i,i]),lower=F)
    }
    # take maximum for each subject
    sensSample = apply(risk,1,max)

    # sample of scores in subjects with no events
    specSample = rmvnorm(nsample,mean=meanScoreSpec,sigma=varScoreSpec)
    # convert to risks
    risk=matrix(nrow=nsample,ncol=ntrait)
    for(i in 1:ntrait) {
      risk[,i]=pnorm(liabThresh[i],mean=specSample[,i],sd=sqrt(1-VX[i,i]),lower=F)
    }
    # take maximum for each subject
    specSample = apply(risk,1,max)

    # concordance
    C = mean(sensSample > specSample)
  }

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              RU=RU))
}
