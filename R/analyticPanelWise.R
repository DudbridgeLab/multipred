#' Analytic panel-wise measures of predictive accuracy
#'
#' Analytic calculation of panel-wise sensitivity, specificity, positive and negative predictive value, concordance and relative utility, under a multivariate liability threshold model.
#'
#' Panel-wise measures consider the prediction of at least one outcome to occur.
#' At least one outcome that did occur must be predicted to occur.
#' For example, panel-wise sensitivity is the probability that, for an individual in which at least one outcome did occur, the predicted risk
#' exceeds the threshold for at least one of the outcomes that did occur.
#' Panel-wise specificity is the probability that, for an individual in which at least one outcome did not occur, the predicted risk is lower than the
#' threshold for all the outcomes that did not occur.
#'
#' Panel-wise concordance is the probability that given one individual in which at least one outcome did occur, and another in which at least one did not occur,
#' the maximum predicted risk over all outcomes that occurred in the former is higher than the maximum over all outcomes that did not occur in the latter.
#' Note that under this definition an individual can be either concordant or discordant with itself.
#' Concordance is calculated by randomly simulating \code{nsample} such pairs of individuals from the specified model.
#'
#' @template analyticParams
#' @template nsample
#'
#' @examples
#' # results will vary due to random sampling in computing multivariate integrals
#' attach(PRSdata)
#' analyticPanelWise(VL,VX,VX,thresh=prevalence,prev=prevalence,nsample=1e5)
#'
#' # $sens
#' # [1] 0.6463497
#'
#' # $spec
#' # [1] 0.0708455
#'
#' # $PPV
#' # [1] 0.1081343
#'
#' # $NPV
#' # [1] 0.9371735
#'
#' # $C
#' # [1] 0.49142
#'
#' # $RU
#' # [1] -0.31006
#'
#' @export
analyticPanelWise = function(VL,VX,VLX=NULL,thresh=NULL,prev,nsample=NULL) {

  # coerce VL, VX and VX to matrices
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

  if (!is.null(nsample)) {
    # concordance
    # jointly simulate liabilities and scores
    LX=rmvnorm(nsample,sigma=rbind(cbind(VL,VLX),cbind(t(VLX),VX)))
    risk=matrix(nrow=nsample,ncol=ntrait)
    for(i in 1:ntrait) {
      risk[,i]=pnorm(liabThresh[i],mean=LX[,i+ntrait],sd=sqrt(1-VX[i,i]),lower=F)
    }
    w = which(apply(t(LX[,1:ntrait]) > liabThresh, 2, max)==1)
    sensSample = apply( t(t(LX[w,1:ntrait]) > liabThresh) * risk[w,], 1,function(x) max(x[x!=0]))
    w = which(apply(t(LX[,1:ntrait]) < liabThresh, 2, max)==1)
    specSample = apply( t(t(LX[w,1:ntrait]) < liabThresh) * risk[w,], 1,function(x) max(x[x!=0]))

    C = mean(sample(sensSample,nsample,replace=TRUE) > sample(specSample,nsample,replace=TRUE))
  }

  if (!is.null(thresh)) {

    # scores corresponding to the risk thresholds
    scoreThresh = liabThresh - qnorm(1-thresh) * sqrt(1-diag(VX))

    # sensitivity

    # probability of no events
    lowerTrunc = rep(-Inf,ntrait)
    upperTrunc = liabThresh
    probNoEvents = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VL)
    attributes(probNoEvents) = NULL

    # loop through event vectors with at least one event
    sens = 0
    for(i in 1:(2^ntrait-1)) {
      event = as.numeric(intToBits(i)[1:ntrait])
      lowerSens = rep(-Inf,ntrait)
      lowerSens[event==1] = liabThresh[event==1]
      upperSens = rep(Inf,ntrait)
      upperSens[event==0] = liabThresh[event==0]
      lowerSens =  c(lowerSens,rep(-Inf,ntrait))
      eventScore = scoreThresh
      eventScore[event==0] = Inf
      upperSens = c(upperSens,eventScore)
      sens = sens + pmvnorm(lower=lowerSens, upper=upperSens, sigma=rbind(cbind(VL,VLX),cbind(t(VLX),VX)))
    }
    sens = 1-sens/(1-probNoEvents)
    attributes(sens) = NULL

    # specificity

    # probability of all events
    lowerTrunc = liabThresh
    upperTrunc = rep(Inf,ntrait)
    probAllEvents = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VL)
    attributes(probAllEvents) = NULL

    # loop through event vectors with at least one non-event
    spec = 0
    for(i in 0:(2^ntrait-2)) {
      event = as.numeric(intToBits(i)[1:ntrait])
      lowerSpec = rep(-Inf,ntrait)
      lowerSpec[event==1] = liabThresh[event==1]
      upperSpec = rep(Inf,ntrait)
      upperSpec[event==0] = liabThresh[event==0]
      lowerSpec =  c(lowerSpec,rep(-Inf,ntrait))
      eventScore = scoreThresh
      eventScore[event==1] = Inf
      upperSpec = c(upperSpec,eventScore)
      spec = spec + pmvnorm(lower=lowerSpec, upper=upperSpec, sigma=rbind(cbind(VL,VLX),cbind(t(VLX),VX)))
    }
    spec = spec/(1-probAllEvents)
    attributes(spec) = NULL

    # relative utility
    # probability of all events given predictions at the threshold
    probAllEventsConditional = pmvnorm(lower=liabThresh, upper=rep(Inf,ntrait), mean=as.vector(VLX %*% solve(VX) %*% as.matrix(thresh)), sigma=VL-VLX %*% solve(VX) %*% t(VLX))
    # probability of no events given predictions at the threshold
    probNoEventsConditional = pmvnorm(lower=rep(-Inf,ntrait), upper=liabThresh, mean=as.vector(VLX %*% solve(VX) %*% as.matrix(thresh)), sigma=VL-VLX %*% solve(VX) %*% t(VLX))

    RU = sens - (1-spec) * (1-probNoEventsConditional)/(1-probAllEventsConditional) * (1-probAllEvents)/(1-probNoEvents)
    attributes(RU) = NULL

    # positive predictive value

    # prevalence of subjects with no predicted events
    lowerTrunc = rep(-Inf,ntrait)
    upperTrunc = scoreThresh
    probNoEvents = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VX)
    attributes(probNoEvents) = NULL

    # loop through event vectors with at least one predicted event
    PPV = 0
    for(i in 1:(2^ntrait-1)) {
      event = as.numeric(intToBits(i)[1:ntrait])
      lowerPPV = rep(-Inf,ntrait)
      lowerPPV[event==1] = scoreThresh[event==1]
      upperPPV = rep(Inf,ntrait)
      upperPPV[event==0] = scoreThresh[event==0]
      lowerPPV =  c(lowerPPV,rep(-Inf,ntrait))
      eventLiab = liabThresh
      eventLiab[event==0] = Inf
      upperPPV = c(upperPPV,eventLiab)
      PPV = PPV + pmvnorm(lower=lowerPPV, upper=upperPPV, sigma=rbind(cbind(VX,t(VLX)),cbind(VLX,VL)))
    }
    PPV = 1-PPV/(1-probNoEvents)
    attributes(PPV) = NULL

    # negative predictive value

    # prevalence of subjects with all predicted events
    lowerTrunc = scoreThresh
    upperTrunc = rep(Inf,ntrait)
    probAllEvents = pmvnorm(lower=lowerTrunc,upper=upperTrunc,sigma=VX)
    attributes(probAllEvents) = NULL

    # loop through event vectors with at least one predicted non-event
    NPV = 0
    for(i in 0:(2^ntrait-2)) {
      event = as.numeric(intToBits(i)[1:ntrait])
      lowerNPV = rep(-Inf,ntrait)
      lowerNPV[event==1] = scoreThresh[event==1]
      upperNPV = rep(Inf,ntrait)
      upperNPV[event==0] = scoreThresh[event==0]
      lowerNPV =  c(lowerNPV,rep(-Inf,ntrait))
      eventLiab = liabThresh
      eventLiab[event==1] = Inf
      upperNPV = c(upperNPV,eventLiab)
      NPV = NPV + pmvnorm(lower=lowerNPV, upper=upperNPV, sigma=rbind(cbind(VX,t(VLX)),cbind(VLX,VL)))
    }
    NPV = NPV/(1-probAllEvents)
    attributes(NPV) = NULL

  }

  return(list(sens=sens,
              spec=spec,
              PPV=PPV,
              NPV=NPV,
              C=C,
              RU=RU))
}

