#' Event-wise measures of predictive accuracy
#'
#' Calculates event-wise sensitivity, specificity, positive and negative predictive value, concordance and net benefit for a vector of predictors.
#'
#' Event-wise measures consider the prediction of individual events summed over individuals.
#' When \code{weight} is the identity matrix, event-wise measures correspond to classical univariate measures with the x matrix vectorised into a column vector.
#' More generally, \code{weight} matrices allow different events to contribute more or less to the calculations, and to allow for co-occurence of events within individuals.

#' @param x Matrix of predicted risks.  Each row corresponds to an individual, each column to a trait.  Each entry should be a risk between 0 and 1.
#' @param y Matrix of traits.  Each row corresponds to an individual, each column to a trait.  Must contain binary events coded as 0 and 1.
#' @param thresh Vector of risk thresholds.  For each row of x, an event is predicted for each trait that exceeds the corresponding element of thresh.  These predictions are then compared to the elements of y.
#' @param prev Vector of prevalences, ie population risks, for each trait.  Defaults to NULL, in which case prevalences are estimated in the data, ignoring ascertainment.
#' all pairs of individuals in x are considered.
#'
#' @export
eventWise=function(x,y,thresh,prev=NULL) {

  # coerce x and y to matrices
  x = as.matrix(x)
  y = as.matrix(y)

  check=checkMatrixDimensions(x,y)
  if (!is.null(check)) {
    print(paste("x and y",check))
    return(NULL)
  }

  check=checkVectorMatrixDimensions(thresh,y)
  if (!is.null(check)) {
    print(paste("thresh",check))
    return(NULL)
  }

  check=checkBinary(y)
  if (!is.null(check)) {
    print(paste("y",check))
    return(NULL)
  }

  # predicted binary traits
  predictedTrait = x
  for(i in 1:dim(x)[1]) predictedTrait[i,] = x[i,]>=thresh

  # if no prevalences given, estimate from the data
  if (is.null(prev)) {
    prev=apply(y,2,mean)
  }

  # sensitivity
  sens = apply(y * predictedTrait,2,sum) / apply(y,2,sum)
  sens = t(prev) %*% sens / sum(prev)

  # specificity
  spec = apply((1-y) * (1-predictedTrait),2,sum) / apply(1-y,2,sum)
  spec = t(1-prev) %*% spec / sum(1-prev)

  # positive predictive value
  PPV = sum(predictedTrait * y) / sum(predictedTrait)

  # negative predictive value
  NPV = sum((1-predictedTrait) * (1-y)) / sum(1-predictedTrait)

  # concordance
  # for each trait
  C = NULL
  for(j in 1:dim(y)[2]) {
    C[j] = as.numeric(wilcox.test(x[y[,j]==1,j],x[y[,j]==0,j],alternative="greater")$statistic/sum(y[,j]==0)/sum(y[,j]==1))
  }
  C = (prev * (1-prev)) %*% C / sum(prev * (1-prev))

  # net benefit
  NB = sens - (1-spec) * min(thresh/(1-thresh)) * sum(1-y) / sum(y)

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       NB=NB
  )
}
