#' Family-wise measures of predictive accuracy
#'
#' Calculates family-wise sensitivity, specificity, positive and negative predictive value, concordance and net benefit for a vector of predictors.
#'
#' Family-wise measures consider the prediction of at least one event within a fixed target vector.
#' If \code{target} is NULL, which is the default, then the marginal measures are calculated across the distribution of
#' trait vectors.  Unlike classical sensitivity, specificity and concordance, the measures then depend upon the
#' distribution of traits in the input data, which may differ from that in the population.
#'
#' If \code{targetProb} is not specified, it is estimated from the \code{y} matrix, ignoring ascertainment.
#'
#' @param target Target trait vector for evaluating family-wise measures of accuracy.
#' @param target2 Second target vector for calculating family-wise concordance.
#' @param targetProb Probability of all the events occurring in \code{target}, used for calculating net benefit.
#' @inheritParams eventWise

#'
#' @export
familyWise=function(x,y,thresh=NULL,target=NULL,target2=NULL,nsample=0,targetProb=NULL) {

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

  if (!is.null(target)) {
    check=checkVectorMatrixDimensions(target,y)
    if (!is.null(check)) {
      print(paste("target",check))
      return(NULL)
    }
  }

  if (!is.null(target2)) {
    check=checkVectorMatrixDimensions(target2,y)
    if (!is.null(check)) {
      print(paste("target2",check))
      return(NULL)
    }
  }

  check=checkBinary(y)
  if (!is.null(check)) {
    print(paste("y",check))
    return(NULL)
  }

  # predicted binary traits
  predictedTrait=x
  for(i in 1:dim(x)[1]) predictedTrait[i,]=x[i,]>=thresh

  # rows of y matching target
  targetIndex = apply(y,1,function(x) sum(x!=target)==0)

  # sensitivity
  sens = mean(y[targetIndex,] %*% t(predictedTrait[targetIndex,])>0)

  # specificity
  spec = mean((1-y[targetIndex,]) %*% t(1-predictedTrait[targetIndex,])>0)

  # rows of predictedTrait matching target
  targetPredictedTrait = apply(predictedTrait,1,function(x) sum(x!=target)==0)

  # positive predictive value
  PPV = mean(predictedTrait[targetPredictedTrait,] %*% t(y[targetPredictedTrait,])>0)

  # negative predictive value
  NPV = mean((1-predictedTrait[targetPredictedTrait,]) %*% t(1-y[targetPredictedTrait,])>0)

  # rows of y matching target2
  target2Index = apply(y,1,function(x) sum(x!=target2)==0)

  # concordance
  # complete enumeration
  if (nsample==0) {
    C = 0
    Ccount = 0
    for(i in which(targetIndex)) {
      for(j in which(target2Index)) {
        if (!identical(y[i,], y[j,])) {
          C = C + ( t(y[i,]>y[j,]) %*% (x[i,]>x[j,]) + t(y[i,]<y[j,]) %*% (x[i,]<x[j,]) >0)
          Ccount = Ccount + 1
        }
      }
    }
    C = as.numeric(C/Ccount)
  }
  # random sampling
  else {
    isample=sample(which(targetIndex),nsample,replace=TRUE)
    jsample=sample(which(target2Index),nsample,replace=TRUE)
    Cnumer = 0
    Cdenom = 0
    for(i in 1:nsample) {
      if (!identical(y[isample[i],], y[jsample[i],])) {
        # select pairs in which one score is systematically higher than the other
        if ( t(y[isample[i],]!=y[jsample[i],]) %*% (x[isample[i],]<x[jsample[i],]) == 0 |
             t(y[isample[i],]!=y[jsample[i],]) %*% (x[isample[i],]>x[jsample[i],]) == 0 ) {
          Cnumer = Cnumer + ( t(y[isample[i],]>y[jsample[i],]) %*% (x[isample[i],]>x[jsample[i],]) +
                    t(y[isample[i],]<y[jsample[i],]) %*% (x[isample[i],]<x[jsample[i],]) >0)
          Cdenom = Cdenom + 1
        }
      }
    }
    C = as.numeric(Cnumer/Cdenom)
    print(Cdenom)

  }

  # net benefit
  # probability of at least one event in the target vector
  if (is.null(targetProb)) {
    yMatch = 0
    for(i in 1:dim(y)[1]) {
      # changed this to allow for NULL target, even though it gives a net benefit of NaN
      #yMatch = yMatch + (target %*% (y[i,]==target) > 0)
      yMatch = yMatch + (y[i,] %*% target > 0)
    }
    targetProb = yMatch/dim(y)[1]
  }

  # cost/benefit ratio
  costBenefit = min(thresh[which(target==1)])
  costBenefit = costBenefit/(1-costBenefit)

  # joint specificity for complementary target
  cTargetIndex = apply(y,1,function(x) sum(x!=(1-target))==0)
  jointSpec = mean((1-y[cTargetIndex,]) %*% t(predictedTrait[cTargetIndex,])==0)
  NB = sens - (1-targetProb)/targetProb*costBenefit * (1-jointSpec)

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       NB=NB
  )
}
