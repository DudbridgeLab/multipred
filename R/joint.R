#' Joint measures of predictive accuracy
#'
#' Calculates joint sensitivity, specificity, positive and negative predictive value, concordance and relative utility for a vector of predictors.
#'
#' Joint measures consider the prediction of all outcomes occuring in an individual.
#' For example, joint sensitivity is the probability that, given an individual in which all outcomes did occur,
#' the predicted risk exceeds the threshold for all outcomes.
#' Joint specificity is the probability that, given an individual in which at least one outcome did not occur,
#' the predicted risk is lower than the threshold for at least one outcome.
#'
#' Joint concordance is the probability that, given one individual in which all outcomes did occur,
#' and another in which at least one outcome did not occur, the minimum predicted risk over all outcomes is greater in the former
#' individual.  It is calculated by randomly drawing such pairs of individuals from \code{y}.
#' If \code{nsample} is zero, all such pairs are drawn from \code{y}; this might be time-consuming.
#' Therefore the default is not to calculate concordance.
#' However, a good estimate of concordance can be obtained from a limited number of random samples \code{nsample}.
#'
#' \code{prev} and \code{condprev} are only required to calculate relative utility, and can be omitted otherwise.

#'
#' @template sharedParams
#' @template nsample
#' @param prev Probability of all events occurring, required for calculating relative utility.
#' If NULL, which is the default, then \code{prev} is estimated from the \code{y} matrix, ignoring ascertainment.
#' @param condprev Probability of all events occuring, conditional on the risk predictions being equal to \code{thresh}.
#' If NULL, which is the default, then \code{prev} is set to the product of the elements of \code{thresh}.
#' This working definition is exact when predictions and outcomes both are jointly independent.
#'
#' @examples
#'
#' attach(PRSdata)
#' joint(risk[,1:2],disease[,1:2],thresh=prevalence[1:2],nsample=1e5)
#'
#' # $sens
#' # [1] 0.4041096
#'
#' # $spec
#' # [1] 0.7653745
#'
#' # $PPV
#' # [1] 0.02488402
#'
#' # $NPV
#' # [1] 0.9885961
#'
#' # $C
#' # [1] 0.63282
#'
#' # $RU
#' # [1] 0.3292956
#'
#' @export
joint=function(x,y,thresh=NULL,prev=NULL,condprev=NULL,nsample=NULL) {

  # coerce x and y to matrices
  x = as.matrix(x)
  y = as.matrix(y)

  check=checkMatrixDimensions(x,y)
  if (!is.null(check)) {
    print(paste("x and y",check))
    return(NULL)
  }

  if (!is.null(thresh)) {
    check=checkVectorMatrixDimensions(thresh,y)
    if (!is.null(check)) {
      print(paste("thresh",check))
      return(NULL)
    }
  }

  check=checkBinary(y)
  if (!is.null(check)) {
    print(paste("y",check))
    return(NULL)
  }

  sens = NULL
  spec = NULL
  PPV = NULL
  NPV = NULL
  C = NULL
  RU = NULL

  nsubject = nrow(y)
  ntrait = ncol(y)
  predictedTrait = matrix(0, nrow=nsubject, ncol=ntrait)

  if (!is.null(thresh)) {

    # predicted binary traits
    for(i in 1:nsubject) predictedTrait[i,]=x[i,]>=thresh

    # sensitivity
    sens = mean(apply(as.matrix(predictedTrait[apply(y,1,min)==1,]),1,min)==1)

    # specificity
    spec = mean(apply(as.matrix(predictedTrait[apply(y,1,min)==0,]),1,min)==0)

    # positive predictive value
    PPV= mean(apply(as.matrix(y[apply(predictedTrait,1,min)==1,]),1,min)==1)

    # negative predictive value
    NPV = mean(apply(as.matrix(y[apply(predictedTrait,1,min)==0,]),1,min)==0)

    # relative utility
    # probability of all events
    if (is.null(prev)) prev = mean(apply(y,1,min)==1)
    # probability of all events, give predictions equal to the threshold
    if (is.null(condprev)) condprev = prod(thresh)

    RU = sens - (1-spec) * condprev/(1-condprev) * (1-prev)/prev

  }

  # concordance
  if (!is.null(nsample)) {
    Cnumer = 0
    Cdenom = 0
    # complete enumeration
    if (nsample==0) {
      # individuals with all events
      for(i in apply(y,1,min)==1) {
        # individuals with at least one non-event
        for(j in apply(y,1,min)==0) {
          Cnumer = Cnumer + (min(x[i,]) > min(x[j,]))
          Cdenom = Cdenom + 1
        }
      }
    }
    # random sampling
    else {
      isample=sample(which(apply(y,1,min)==1),nsample,replace=TRUE)
      jsample=sample(which(apply(y,1,min)==0),nsample,replace=TRUE)
      for(i in 1:nsample) {
        Cnumer = Cnumer + (min(x[isample[i],]) > min(x[jsample[i],]))
      }
      Cdenom = nsample
    }
    C = as.numeric(Cnumer / Cdenom)
  }


  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       RU=RU
  )
}
