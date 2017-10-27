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
#' @param weight Weighting matrix. Defaults to identity.
#' @param sample Number of random samples to draw when estimating concordance.  Defaults to 0, in which case
#' all pairs of individuals in x are considered.
#'
#' @export
eventWise=function(x,y,thresh,weight=NULL,sample=0) {

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

  check=checkVectorMatrixDimensions(weight[1,],y)
  if (!is.null(check)) {
    print(paste("weight",check))
    return(NULL)
  }

  check=checkBinary(y)
  if (!is.null(check)) {
    print(paste("y",check))
    return(NULL)
  }

  # predicted binary traits
  predictedTrait=x
  for(i in 1:dim(x)[1]) predictedTrait[i,]=x[i,]>=thresh

  # weight matrix is identity if not specified as input
  if (is.null(weight)) weight=diag(length(thresh))

  # sensitivity
  sens = sum((y %*% weight) * predictedTrait) / sum((y %*% weight) * y)

  # specificity
  spec = sum(((1-y) %*% weight) * (1-predictedTrait)) / sum(((1-y) %*% weight) * (1-y))

  # positive predictive value
  PPV = sum((predictedTrait %*% weight) * y) / sum((predictedTrait %*% weight) * predictedTrait)

  # negative predictive value
  NPV = sum(((1-predictedTrait) %*% weight) * (1-y)) / sum(((1-predictedTrait) %*% weight) * (1-predictedTrait))

  # concordance
  # identity weights
  if (sample==0 & identical(weight,diag(1,length(thresh)))) {
    C=as.numeric(wilcox.test(x[y==1],x[y==0],alternative="greater")$statistic/sum(y==0)/sum(y==1))
  }
  else {
    Cnumer=0
    Cdenom=0
    # complete enumeration
    if (sample==0) {
      for(i in 1:dim(y)[1]) {
        xvector=NULL
        yvector=NULL
        for(j in 1:dim(y)[1]) {
          xvector=x[i,]>x[j,]
          yvector=y[i,]>y[j,]
          Cnumer=Cnumer+t(xvector) %*% weight %*% yvector
          Cdenom=Cdenom+t(yvector) %*% weight %*% yvector
          xvector=x[i,]<x[j,]
          yvector=y[i,]<y[j,]
          Cnumer=Cnumer+t(xvector) %*% weight %*% yvector
          Cdenom=Cdenom+t(yvector) %*% weight %*% yvector
        }
      }
    }
    #random sampling
    else {
      isample=sample(dim(y)[1],sample,replace=TRUE)
      jsample=sample(dim(y)[1],sample,replace=TRUE)
      for(i in 1:sample) {
        xvector=x[isample[i],]>x[jsample[i],]
        yvector=y[isample[i],]>y[jsample[i],]
        Cnumer=Cnumer+t(xvector) %*% weight %*% yvector
        Cdenom=Cdenom+t(yvector) %*% weight %*% yvector
        xvector=x[isample[i],]<x[jsample[i],]
        yvector=y[isample[i],]<y[jsample[i],]
        Cnumer=Cnumer+t(xvector) %*% weight %*% yvector
        Cdenom=Cdenom+t(yvector) %*% weight %*% yvector
      }
    }
    C=as.numeric(Cnumer/Cdenom)
  }

  # net benefit
  NB = sens - (1-spec) * min(diag(weight)*thresh/(1-diag(weight)*thresh)) * sum(((1-y) %*% weight) * (1-y)) / sum((y %*% weight) * y)

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       NB=NB
  )
}
