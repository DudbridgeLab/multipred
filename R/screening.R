#' Screening measures of predictive accuracy
#'
#' Calculates screening sensitivity, specificity, positive and negative predictive value, concordance and relative utility for a vector of predictors.
#'
#' Screening measures consider the prediction of at least one outcome to occur, without regard to whether the correct outcomes are predicted.
#' For example, screening sensitivity is the probability that, for an individual in which at least one outcome did occur, the predicted risk
#' exceeds the threshold for at least one outcome (but not necessary the ones that did occur).
#' Screening specificity is the probability that, for an individual in which no outcomes did occur, the predicted risk is lower than the
#' threshold for all outcomes.
#'
#' Screening concordance is the probability that given one individual in which at least one outcome did occur, and another in which no outcomes did occur,
#' the maximum predicted risk over all outcomes is higher in the former individual.
#' pairs of individuals from \code{y}.  If \code{nsample} is zero, all such pairs are drawn from \code{y};
#' this might be time-consuming.  Therefore the default is not to calculate concordance.
#' However, a good estimate can be obtained from a limited number of random samples \code{nsample}.
#'
#' \code{prev} and \code{condprev} are only required to calculate relative utility, and can be omitted otherwise.

#' @template sharedParams
#' @template nsample
#' @param prev Probability of at least one event, required for calculating relative utility.
#' If NULL, which is the default, then \code{prev} is estimated from the \code{y} matrix, ignoring ascertainment.
#' @param condprev Probability of at least one events, conditional on the risk predictions being equal to \code{thresh}.
#' If NULL, which is the default, then \code{prev} is set to 1- the product of the elements of (1-\code{thresh}).
#' This working definition is exact when predictions and outcomes both are jointly independent.
#'
#'
#' @examples
#'
#' attach(PRSdata)
#' screening(risk,disease,thresh=prevalence,nsample=1e5)
#'
#' # $sens
#' # [1] 0.9653894
#'
#' # $spec
#' # [1] 0.05989024
#'
#' # $PPV
#' # [1] 0.1654311
#'
#' # $NPV
#' # [1] 0.8996416
#'
#' # $C
#' # [1] 0.59156
#'
#' # $RU
#' # [1] -0.009099365
#'
#' @export
screening=function(x,y,thresh=NULL,prev=NULL,condprev=NULL,nsample=NULL) {

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
    for(i in 1:nsubject) for(j in 1:ntrait) predictedTrait[i,j]=x[i,j]>=thresh[j]

    # sensitivity
    sens = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==1,]),1,max)==1)

    # specificity
    spec = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==0,]),1,max)==0)

    # positive predictive value
    PPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==1,]),1,max)==1)

    # negative predictive value
    NPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==0,]),1,max)==0)

    # relative utility
    # probability of at least one event
    if (is.null(prev)) prev = mean(apply(y,1,max))
    # probability of at least one event, given predictions equal to the threshold
    if (is.null(condprev)) condprev = (1-prod(1-thresh))

    RU = sens - (1-spec) * condprev/(1-condprev) * (1-prev)/prev

  }

  # concordance
  if (!is.null(nsample)) {
    Cnumer = 0
    Cdenom = 0
    # complete enumeration
    if (nsample==0) {
      # individuals with at least one event
      for(i in which(apply(y,1,max)==1)) {
        # individuals with no events
        for(j in which(apply(y,1,max)==0)) {
          Cnumer = Cnumer + (max(x[i,]) > max(x[j,]))
          Cdenom = Cdenom + 1
        }
      }
    }
    # random sampling
    if (nsample>0) {
      sensSample=NULL
      specSample=NULL
      predictedTrait1=predictedTrait
      isample=sample(which(apply(y,1,max)==1),nsample,replace=TRUE)
      jsample=sample(which(apply(y,1,max)==0),nsample,replace=TRUE)
      for(i in 1:nsample) {
        Cnumer = Cnumer + (max(x[isample[i],]) > max(x[jsample[i],]))
      }
      Cdenom = nsample
    }
    C = as.numeric(Cnumer/Cdenom)
  }

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       RU=RU
  )
}
