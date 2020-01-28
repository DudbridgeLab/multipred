#' Outcome-wise measures of predictive accuracy
#'
#' Calculates outcome-wise sensitivity, specificity, positive and negative predictive value, concordance and relative utility for a vector of predictors.
#'
#' Outcome-wise measures consider the prediction of individual outcomes summed over individuals.
#' When \code{weight} is a vector of 1's (default), outcome-wise measures correspond to classical univariate measures with the \code{x} matrix vectorised into a column vector.
#' More generally, \code{weight} allows different outcomes to contribute more or less to the calculations.
#'
#' Outcome-wise sensitivity, specificity and concordance are weighted sums of the univariate measures,
#' where the weights depend on \code{prev}.  By default, \code{prev} is estimated from the outcome
#' rates in \code{y}, but external estimates of population risk may be used instead.
#'
#' @param prev Vector of prevalences, ie population risks, for each trait.  Defaults to NULL, in which case prevalences are estimated in the data, ignoring ascertainment.
#' @param weight Vector of weights.  Defaults to a vector of 1's.
#'
#' @template sharedParams
#'
#' @examples
#'
#' attach(PRSdata)
#' outcomeWise(risk,disease,thresh=prevalence)
#'
#' # $sens
#' # [1] 0.6017748
#'
#' # $spec
#' # [1] 0.6129354
#'
#' # $PPV
#' # [1] 0.04595316
#'
#' # $NPV
#' # [1] 0.9802688
#'
#' # $C
#' # [1] 0.6442582
#'
#' # $RU
#' # [1] 0.2251043
#'
#' @export
outcomeWise=function(x,y,thresh=NULL,weight=NULL,prev=NULL) {

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

  if (!is.null(weight)) {
    check=checkVectorMatrixDimensions(weight,y)
    if (!is.null(check)) {
      print(paste("weight",check))
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

  if (is.null(weight)) weight = rep(1,ntrait)

  # if no prevalences given, estimate from the data
  if (is.null(prev)) prev = apply(y,2,mean)

  if (!is.null(thresh)) {

    # predicted binary traits
    for(i in 1:nsubject) predictedTrait[i,] = x[i,]>=thresh

    # sensitivity
    sens = apply(y * predictedTrait,2,sum) / apply(y,2,sum)
    sens = t(prev * weight) %*% sens / sum(prev * weight)
    sens = as.vector(sens)

    # specificity
    spec = apply((1-y) * (1-predictedTrait),2,sum) / apply(1-y,2,sum)
    spec = t((1-prev) * weight) %*% spec / sum((1-prev) * weight)
    spec = as.vector(spec)

    # positive predictive value
    probPredicted = apply(predictedTrait,2,mean)
    PPV = apply(y*predictedTrait,2,sum) / apply(predictedTrait,2,sum)
    PPV = t(probPredicted * weight) %*% PPV / sum(probPredicted * weight)
    PPV = as.vector(PPV)

    # negative predictive value
    NPV = apply((1-y) * (1-predictedTrait),2,sum) / apply(1-predictedTrait,2,sum)
    NPV = t((1-probPredicted) * weight) %*% NPV / sum((1-probPredicted) * weight)
    NPV = as.vector(NPV)

    # relative utility
    RU = sens - (thresh %*% weight) / ((1-thresh) %*% weight) * ((1-prev) %*% weight) / (prev %*% weight) * (1-spec)
    RU = as.vector(RU)

  }

  # concordance
  # for each trait
  C = NULL
  for(j in 1:ntrait) {
    C[j] = as.numeric(wilcox.test(x[y[,j]==1,j],x[y[,j]==0,j],alternative="greater")$statistic/sum(y[,j]==0)/sum(y[,j]==1))
  }
  # combined
  C = (prev * (1-prev) * weight) %*% C / sum(prev * (1-prev) * weight)
  C = as.vector(C)

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       RU=RU
  )
}
