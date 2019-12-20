#' Analytic joint measures of predictive accuracy
#'
#' Analytic calculation of joint sensitivity, specificity, positive and negative predictive value, concordance and relative utility, under a multivariate liability threshold model.
#'
#' Joint measures consider the prediction of all events occuring in an individual.
#' For example, joint sensitivity is the probability that, given an individual with events for all traits,
#' the predicted risk exceeds the threshold for all traits.
#' Joint specificity is the probability that, given an individual with at least one non-event,
#' the predicted risk is lower than the threshold for at least one trait.
#'
#' Joint concordance is the probability that, given one individual with events for all traits
#' and another with at least one non-event, the minimum predicted risk over all traits is greater in the former
#' individual.  It is calculated by randomly simulating \code{nsample} such pairs of individuals from the specified model.
#'
#' @template analyticParams
#' @template nsample
#'
#' @export
analyticJoint = function(VL,VX,VLX=NULL,thresh=NULL,prev,nsample=NULL) {

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

  # truncation points for multivariate normal liability
  lowerTrunc = liabThresh
  upperTrunc = rep(Inf,ntrait)

  if (!is.null(thresh)) {

    # scores corresponding to the risk thresholds
    scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VX))

    # probability of all events
    probAllEvents = pmvnorm(lower=liabThresh, upper=rep(Inf,ntrait), sigma=VL)

    # probability of all predicted events
    probAllPredicted = pmvnorm(lower=scoreThresh, upper=rep(Inf,ntrait), sigma=VX)

    # probability of all events and all predicted events
    probAllEventsAllPredicted = pmvnorm(lower=c(liabThresh,scoreThresh), upper=rep(Inf,2*ntrait), sigma=rbind(cbind(VL,VLX),cbind(t(VLX),VX)))

    # sensitivity
    sens = probAllEventsAllPredicted / probAllEvents
    attributes(sens) = NULL

    # specificity
    spec = 1 - (probAllPredicted-probAllEventsAllPredicted) / (1-probAllEvents)
    attributes(spec) = NULL

    # relative utility
    # probability of all events given predictions at the threshold
    probAllEventsConditional = pmvnorm(lower=liabThresh, upper=rep(Inf,ntrait), mean=as.vector(VLX %*% solve(VX) %*% as.matrix(thresh)), sigma=VL-VLX %*% solve(VX) %*% t(VLX))
    #RU = sens - (1-spec) * prod(thresh)/(1-prod(thresh)) * (1-probAllEvents)/probAllEvents
    RU = sens - (1-spec) * probAllEventsConditional/(1-probAllEventsConditional) * (1-probAllEvents)/probAllEvents
    attributes(RU) = NULL

    # PPV
    PPV = probAllEventsAllPredicted / probAllPredicted
    attributes(PPV) = NULL

    # NPV
    NPV = 1 - (probAllEvents-probAllEventsAllPredicted) / (1-probAllEvents)
    attributes(NPV) = NULL
  }

  #concordance
  if (!is.null(nsample)) {

    # probability of all events
    probAllEvents = pmvnorm(lower=liabThresh,upper=rep(Inf,ntrait),sigma=VL)
    attributes(probAllEvents) = NULL

    # mean liability in subjects with all events
    liabTrunc = tmvtnorm::mtmvnorm(sigma=VL, lower=liabThresh, upper=rep(Inf,ntrait))
    # in subjects with at least one non-event
    liabTruncMean = -probAllEvents/(1-probAllEvents) * liabTrunc$tmean

    # variance of liability in subjects with at least one non-event
    liabTruncVar = (VL - probAllEvents * (liabTrunc$tvar + liabTrunc$tmean %*% t(liabTrunc$tmean))) / (1-probAllEvents) -
      probAllEvents^2/(1-probAllEvents)^2 * liabTrunc$tmean %*% t(liabTrunc$tmean)

    # mean score in subjects with all events
    invVL = solve(VL)
    meanScoreSens = as.vector(VLX %*% invVL %*% liabTrunc$tmean)
    # covariance matrix of scores in subjects with all events
    varScoreSens = VX - VLX %*% (invVL - invVL %*% liabTrunc$tvar %*% invVL) %*% t(VLX)
    # make varScoreSens symmetrical
    varScoreSens = (varScoreSens + t(varScoreSens)) / 2

    # mean score in subjects with at least one non-event
    meanScoreSpec = as.vector(VX %*% invVL %*% liabTruncMean)
    varScoreSpec = VX - VLX %*% (invVL - invVL %*% liabTruncVar %*% invVL) %*% t(VLX)
    # make varScoreSpec symmetrical
    varScoreSpec = (varScoreSpec + t(varScoreSpec)) / 2

    # sample of scores in subjects with all events
    # convert to risks
    sensSample = rmvnorm(nsample,mean=meanScoreSens,sigma=varScoreSens)
    risk=matrix(nrow=nsample,ncol=ntrait)
    for(i in 1:ntrait) {
      risk[,i]=pnorm(liabThresh[i],mean=sensSample[,i],sd=sqrt(1-VX[i,i]),lower=F)
    }
    # take minimum over all subjects
    sensSample = apply(risk,1,min)

    # sample of scores in subjects with at least one non-event
    specSample = rmvnorm(nsample,mean=meanScoreSpec,sigma=varScoreSpec)
    risk=matrix(nrow=nsample,ncol=ntrait)
    for(i in 1:ntrait) {
      risk[,i]=pnorm(liabThresh[i],mean=specSample[,i],sd=sqrt(1-VX[i,i]),lower=F)
    }
    # take minimum over all subjects
    specSample = apply(risk,1,min)

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
