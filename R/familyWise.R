#' Family-wise measures of predictive accuracy
#'
#' Calculates family-wise sensitivity, specificity, positive and negative predictive value, concordance and relative utility for a vector of predictors.
#'
#' Family-wise measures consider the prediction of at least one outcome to occur.
#' An outcome that did occur must be predicted to occur in at least one case.
#' For example, family-wise sensitivity is the probability that, for an individual in which at least one outcome did occur, the predicted risk
#' exceeds the threshold for at least one of the outcomes that did occur.
#' Family-wise specificity is the probability that, for an individual in which at least one outcome did not occur, the predicted risk is lower than the
#' threshold for all the outcomes that did not occur.
#'
#' Family-wise concordance is the probability that given one individual in which at least one outcome did occur, and another in which at least one did not occur,
#' the maximum predicted risk over all outcomes that occurred in the former is higher than the maximum over all outcomes that did not occur in the latter.
#' Note that under this definition an individual can be either concordant or discordant with itself.
#' Concordance is calculated by randomly drawing such pairs of individuals from \code{y}.
#' If \code{nsample} is zero, all such pairs are drawn from \code{y}; this might be time-consuming.
#' Therefore the default is not to calculate condcordance.
#' However, a good estimate of concordance can be obtained from a limited number of random samples \code{nsample}.
#'
#' \code{prev0}, \code{prev1}, \code{condprev0} and \code{condprev1} are only required to calculate relative utility, and can be omitted otherwise.

#' @template sharedParams
#' @template nsample
#' @param prev0 Probability of at least one non-event, required for calculating relative utility.
#' If NULL, which is the default, then \code{prev0} is estimated from the \code{y} matrix, ignoring ascertainment.
#' @param prev1 Probability of at least one event, required for calculating relative utility.
#' If NULL, which is the default, then \code{prev1} is estimated from the \code{y} matrix, ignoring ascertainment.
#' @param condprev0 Probability of at least one non-event, conditional on the risk predictions being equal to \code{thresh}.
#' If NULL, which is the default, then \code{prev} is set to 1 - the product of the elements of \code{thresh}.
#' This working definition is exact when predictions and outcomes both are jointly independent.
#' @param condprev1 Probability of at least one event, conditional on the risk predictions being equal to \code{thresh}.
#' If NULL, which is the default, then \code{prev} is set to 1 - the product of the elements of (1-\code{thresh}).
#' This working definition is exact when predictions and outcomes both are jointly independent.
#'
#' @examples
#'
#' attach(PRSdata)
#' familyWise(risk,disease,thresh=prevalence,nsample=1e5)
#'
#' # $sens
#' # [1] 0.6266996
#'
#' # $spec
#' # [1] 0.0701
#'
#' # $PPV
#' # [1] 0.1073925
#'
#' # $NPV
#' # [1] 0.9311696
#'
#' # $C
#' # [1] 0.47999
#'
#' # $RU
#' # [1] -5.120519
#'
#' @export
familyWise=function(x,y,thresh=NULL,prev0=NULL,prev1=NULL,condprev0=NULL,condprev1=NULL,nsample=NULL) {

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
    for(i in 1:nsubject) predictedTrait[i,] = x[i,]>=thresh

    # sensitivity
    sens = sum(apply(y*predictedTrait,1,max)) / sum(apply(y,1,max))

    # specificity
    spec = 1- sum(apply((1-y)*predictedTrait,1,max)) / sum(apply(1-y,1,max))

    # positive predictive value
    PPV = sum(apply(predictedTrait*y,1,max)) / sum(apply(predictedTrait,1,max))

    # negative predictive value
    NPV = 1- sum(apply((1-predictedTrait)*y,1,max)) / sum(apply(1-predictedTrait,1,max))

    # relative utility
    # probability of at least one non-event
    if (is.null(prev0)) prev0 = 1-mean(apply(y,1,min))
    # probability of at least one event
    if (is.null(prev1)) prev1 = mean(apply(y,1,max))
    # probability of at least one non-event, given predictions at the threshold
    if (is.null(condprev0)) condprev0 = 1-prod(thresh)
    # probability of at least one event, given predictions at the threshold
    if (is.null(condprev1)) condprev1 = 1-prod(1-thresh)
    RU = sens - (1-spec) * condprev1/condprev1 * prev0/prev1

  }

  # concordance
  if (!is.null(nsample)) {
    Cnumer = 0
    Cdenom = 0
    # complete enumeration
    if (nsample==0) {
      for(i in which(apply(y,1,max)==1)) {
        for(j in which(apply(y,1,min)==0)) {
          Cnumer = Cnumer + ( max(x[i,y[i,]==1]) > max(x[j,y[j,]==0]) )
          Cdenom = Cdenom + 1
        }
      }
    }
    # random sampling
    else {
      isample=sample(which(apply(y,1,max)==1),nsample,replace=TRUE)
      jsample=sample(which(apply(y,1,min)==0),nsample,replace=TRUE)

      Cnumer = 0
      Cdenom = 0
      for(i in 1:nsample) {
        Cnumer = Cnumer + ( max(x[isample[i],y[isample[i],]==1]) > max(x[jsample[i],y[jsample[i],]==0]) )
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
